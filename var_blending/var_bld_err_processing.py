import os
from datetime import timedelta
from var_blending import fft_process
from var_blending import MODELTYPES, BLD_VAR_ERR_RATIO, WET_CRITERIAS
from wrf import getvar
from netCDF4 import Dataset
import numpy
from scipy.sparse import csc_matrix
import gzip
import cPickle
from scipy.ndimage.filters import gaussian_filter1d
from scipy import signal

try:
    from AMPSAws import resources
except:
    pass


def wet_criteria(datapath, max_wavenumber,
                 lowest_land_height=5.0,
                 lowest_dbz_consider_dry=0.25,
                 max_allowed_large_scale_power_from_hist=None,
                 max_allowed_wet_ratio_from_hist=None):
    """whether to calculate model error or run blending based on
    (1) the ratio of the large scale information in dbz
    (2) the ratio of rainfall areas over lands
    """
    ncfile = Dataset(datapath)
    mdbz = getvar(ncfile, 'mdbz').values
    topo = getvar(ncfile, "HGT").values
    ncfile.close()

    # the ratio of the large scale information in dbz
    data_shape = mdbz.shape
    model_power_along_lat = []
    for i in range(0, data_shape[0]):
        cur_dbz = mdbz[i, :]
        model_power_along_lat.append(signal.welch(
            cur_dbz, nperseg=data_shape[1])[1])

    model_power_along_lat = numpy.asarray(model_power_along_lat)
    model_power_along_lat = numpy.mean(model_power_along_lat, 0)

    total_power = numpy.sum(model_power_along_lat)
    large_scale_power = numpy.sum(model_power_along_lat[0:max_wavenumber])
    large_scale_power_ratio = large_scale_power/total_power

    # the ratio of rainfall areas over lands
    mdbz[topo < lowest_land_height] = numpy.NaN
    wet_points = numpy.nansum(mdbz > lowest_dbz_consider_dry)
    dry_points = numpy.nansum(mdbz < lowest_dbz_consider_dry)
    wet_ratio = wet_points/float(dry_points + wet_points)

    # if (1) LAM does not have strong large scale forcing
    #        (e.g., by determined the dbz at spectrual space) and
    #    (2) the rainfalll area over lands are smaller than a certain value,
    # => blending/error calculation will be carried out
    if (not max_allowed_large_scale_power_from_hist) or \
       (not max_allowed_wet_ratio_from_hist):
        # used in the step of checking the data for creating errors
        if large_scale_power_ratio < \
            WET_CRITERIAS['max_allowed_large_scale_power'] \
                and wet_ratio < WET_CRITERIAS['max_allowed_wet_ratio']:
            return True, large_scale_power_ratio, wet_ratio
        else:
            return False, large_scale_power_ratio, wet_ratio
    else:
        # used in the step of checking the data for blending
        max_allowed_large_scale_power = min(
            WET_CRITERIAS['max_allowed_large_scale_power'],
            max_allowed_large_scale_power_from_hist * 1.3)

        max_allowed_wet_ratio = min(
            WET_CRITERIAS['max_allowed_wet_ratio'],
            max_allowed_wet_ratio_from_hist * 1.3)

        if (large_scale_power_ratio <= max_allowed_large_scale_power) and \
           (wet_ratio <= max_allowed_wet_ratio):
            return True
        else:
            return False


def obtain_griddata(datapath, model_variable):
    # Open the NetCDF file
    ncfile = Dataset(datapath)

    model_height = getvar(ncfile, "height").values
    topo = getvar(ncfile, "HGT").values
    if 'UV' == model_variable:
        var = getvar(ncfile, 'wspd_wdir')[0].values
    elif 'T' == model_variable:
        var = getvar(ncfile, 'temp').values
    else:
        var = getvar(ncfile, model_variable).values

    ncfile.close()
    return var, topo, model_height


def download_data_from_s3_to_local(
        work_dir, cur_val_datetime,
        global_fcst_on_s3, global_ana_on_s3,
        lam_fcst_on_s3, lam_model_name, domain_id,
        fcst_length):

    def _copy_from_s3_to_local(data_on_s3, data_on_local):
        if not os.path.exists(data_on_local):
            try:
                print 'copy glb from {} to {}'.format(
                    data_on_s3, data_on_local)
                resources.copy(data_on_s3, data_on_local)
            except:
                print ' --- failed coping {} to {}'.format(
                    data_on_s3, data_on_local)
                pass

    cur_hist_dir = os.path.join(work_dir, 'hist_data',
                                cur_val_datetime.strftime('%Y%m%d%H'))

    if not os.path.exists(cur_hist_dir):
        os.makedirs(cur_hist_dir)

    # 0.0: copy global fcst
    cur_ana_datetime = cur_val_datetime - timedelta(seconds=3600*fcst_length)
    cur_hist_glb_fcst_path_on_s3 = os.path.join(
        global_fcst_on_s3,
        'global_a{}_v{}.nc'.format(cur_ana_datetime.strftime('%Y%m%d%H'),
                                   cur_val_datetime.strftime('%Y%m%d%H')))
    cur_hist_glb_fcst_path_on_local = os.path.join(
        cur_hist_dir, 'glb_l{}_d0{}.nc'.format(fcst_length, domain_id))

    _copy_from_s3_to_local(cur_hist_glb_fcst_path_on_s3,
                           cur_hist_glb_fcst_path_on_local)

    # 0.1: copy global analysis
    cur_hist_glb_ana_path_on_s3 = os.path.join(
        global_ana_on_s3,
        'global_a{}_v{}.nc'.format(cur_val_datetime.strftime('%Y%m%d%H'),
                                   cur_val_datetime.strftime('%Y%m%d%H')))
    cur_hist_glb_ana_path_on_local = os.path.join(
        cur_hist_dir, 'ana_l0_d0{}.nc'.format(domain_id))
    _copy_from_s3_to_local(cur_hist_glb_ana_path_on_s3,
                           cur_hist_glb_ana_path_on_local)

    # 0.2: copy lam fcst
    cur_hist_lam_fcst_path_on_s3 = os.path.join(
        lam_fcst_on_s3, lam_model_name,
        cur_ana_datetime.strftime('%y'),
        cur_ana_datetime.strftime('%m'), cur_ana_datetime.strftime('%d'),
        cur_ana_datetime.strftime('%H'),
        'wrf_hourly_{}_d0{}_{}'.format(
             lam_model_name,
             domain_id,
             cur_val_datetime.strftime('%Y-%m-%d_%H:00:00')))
    cur_hist_lam_fcst_path_on_local = os.path.join(
        cur_hist_dir, 'lam_l{}_d0{}.nc'.format(fcst_length, domain_id))
    _copy_from_s3_to_local(cur_hist_lam_fcst_path_on_s3,
                           cur_hist_lam_fcst_path_on_local)


def init_err_vars(use_inv_fft_as_err):
    gridata = {}
    fft_coef_useful = {}
    power_spectrum_useful = {}
    freq_rows_useful = {}
    freq_cols_useful = {}
    data_number = 0

    lam_wet_ratio = {}
    lam_large_scale_power_ratio = {}

    if use_inv_fft_as_err:
        inv_fft_err = {}
    else:
        inv_fft_err = None

    griddata_shape = {}

    return (gridata, fft_coef_useful, power_spectrum_useful,
            freq_rows_useful, freq_cols_useful, data_number,
            inv_fft_err, griddata_shape,
            lam_wet_ratio, lam_large_scale_power_ratio)


def init_all_data_available(fcst_length):
    cur_data_available = {}
    for modeltype in MODELTYPES:
        cur_data_available[modeltype] = {}
        cur_data_available[modeltype][fcst_length] = True
    cur_data_available['ana'] = {}
    cur_data_available['ana'][0] = True
    return cur_data_available


def extract_model_data_and_err(
                       start_datetime, end_datetime,
                       work_dir, model_variable,
                       all_data_available,
                       fcst_length, analysis_interval,
                       gridata, data_number,
                       max_wavenumber,
                       fft_coef_useful, power_spectrum_useful,
                       freq_rows_useful, freq_cols_useful, inv_fft_err,
                       total_model_levels,
                       griddata_shape, domain_id,
                       lam_wet_ratio, lam_large_scale_power_ratio,
                       use_wet_criteria=False,
                       adjust_lam_err_use_topo=False,
                       use_inv_fft_as_err=False):

    def _check_all_data_available(
            cur_data_available, fcst_length):
        for modeltype in MODELTYPES:
            if not cur_data_available[modeltype][fcst_length]:
                return False

        if not cur_data_available['ana'][0]:
            return False
        return True

    def _adjust_lam_err_over_topo(model_height, topo,
                                  model_field, fcst_length):
        for i in range(0, model_field['lam'][fcst_length].shape[0]):
            cur_diff = model_height[i, :, :] - topo
            model_field['lam'][fcst_length][i, :, :][cur_diff < topo] = \
                model_field['ana'][0][i, :, :][cur_diff < topo]

        return model_field

    """-----------------------------
    # main function starts
    -----------------------------"""
    start_val_datetime = start_datetime
    end_val_datetime = end_datetime
    cur_val_datetime = start_val_datetime

    if model_variable in ['U', 'V']:
        gridata_wind_spd = {}
    else:
        gridata_wind_spd = None
        fft_coef_wind_spd = None

    while cur_val_datetime <= end_val_datetime:
        print cur_val_datetime
        if not _check_all_data_available(all_data_available[
                    cur_val_datetime.strftime('%Y%m%d%H')], fcst_length):
            cur_val_datetime = cur_val_datetime + timedelta(
                seconds=3600*analysis_interval)
            continue

        print '{}: found full dataset for {}'.format(
            data_number, cur_val_datetime)
        if use_wet_criteria:
            print '    step 0: wet criteria checks'
            datapath = os.path.join(
                work_dir,  'hist_data', cur_val_datetime.strftime('%Y%m%d%H'),
                'lam_l{}_d0{}.nc'.format(fcst_length, domain_id))
            pass_wet_check, large_scale_power_ratio, wet_ratio = wet_criteria(
                datapath, max_wavenumber)
            if not pass_wet_check:
                print ('     wet criteria check failed: '
                       'LAM large scale power ratio :{}; '
                       'wet ratio: {}').format(
                        large_scale_power_ratio, wet_ratio)
                cur_val_datetime = cur_val_datetime + timedelta(
                    seconds=3600*analysis_interval)
                continue

        print '    step 1: extracting gridded dataset'
        if all_data_available:
            gridata[cur_val_datetime.strftime('%Y%m%d%H')] = {}
            # read fcst
            for modeltype in MODELTYPES:
                gridata[cur_val_datetime.strftime('%Y%m%d%H')][modeltype] = {}
                gridata[cur_val_datetime.strftime('%Y%m%d%H')][
                    modeltype][fcst_length] = {}
                model_path = os.path.join(
                    work_dir,  'hist_data',
                    cur_val_datetime.strftime('%Y%m%d%H'),
                    '{}_l{}_d0{}.nc'.format(modeltype, fcst_length, domain_id))

                (gridata[cur_val_datetime.strftime('%Y%m%d%H')][
                    modeltype][fcst_length],
                 topo, model_height) = obtain_griddata(
                    model_path, model_variable)

                if model_variable in ['U', 'V']:
                    gridata_wind_spd[modeltype] = {}
                    (gridata_wind_spd[modeltype][fcst_length],
                     _, _) = obtain_griddata(model_path, 'UV')

            # read analysis
            gridata[cur_val_datetime.strftime('%Y%m%d%H')]['ana'] = {}
            gridata[cur_val_datetime.strftime('%Y%m%d%H')]['ana'][0] = {}
            model_path = os.path.join(
                work_dir,  'hist_data',
                cur_val_datetime.strftime('%Y%m%d%H'),
                'ana_l0_d0{}.nc'.format(domain_id))
            (gridata[cur_val_datetime.strftime('%Y%m%d%H')]['ana'][0],
             topo, model_height) = obtain_griddata(model_path, model_variable)
            if model_variable in ['U', 'V']:
                gridata_wind_spd['ana'] = {}
                gridata_wind_spd['ana'][0], _, _ = obtain_griddata(
                    model_path, 'UV')

            # read model shape
            griddata_shape = gridata[cur_val_datetime.strftime('%Y%m%d%H')][
                modeltype][fcst_length].shape

        if adjust_lam_err_use_topo:
            gridata[cur_val_datetime.strftime('%Y%m%d%H')] = \
                _adjust_lam_err_over_topo(
                    model_height, topo,
                    gridata[cur_val_datetime.strftime('%Y%m%d%H')],
                    fcst_length)

        print '    step 2: extracting FFT dataset'
        (cur_fft_coef_useful, cur_power_spectrum_useful,
         cur_freq_rows_useful, cur_freq_cols_useful) = \
            fft_process.obtain_fftdata(
                gridata[cur_val_datetime.strftime('%Y%m%d%H')],
                MODELTYPES, fcst_length, 0, total_model_levels)
        if model_variable in ['U', 'V']:
            (fft_coef_wind_spd, _, _, _) = fft_process.obtain_fftdata(
                gridata_wind_spd,
                MODELTYPES, fcst_length, 0, total_model_levels)

        print '    step 3: extract power spectral using welch method'
        power_spectrum_welch = fft_process.obtain_power_spectrum_using_welch(
            gridata[cur_val_datetime.strftime('%Y%m%d%H')],
            fcst_length, MODELTYPES,
            total_model_levels, max_wavenumber)

        if use_inv_fft_as_err:
            print '    step 4: extract RMS for each freqencies'
            inv_fft_err[cur_val_datetime.strftime('%Y%m%d%H')] = {}
            inv_fft_err[cur_val_datetime.strftime('%Y%m%d%H')] = \
                fft_process.obtain_inv_fft_rms(
                    cur_fft_coef_useful, fcst_length,
                    model_variable, total_model_levels,
                    max_wavenumber_threshold=max_wavenumber,
                    fft_coef_wind_spd=fft_coef_wind_spd)

        fft_coef_useful[
            cur_val_datetime.strftime('%Y%m%d%H')] = cur_fft_coef_useful
        power_spectrum_useful[
            cur_val_datetime.strftime('%Y%m%d%H')] = power_spectrum_welch
        freq_rows_useful[
            cur_val_datetime.strftime('%Y%m%d%H')] = cur_freq_rows_useful
        freq_cols_useful[
            cur_val_datetime.strftime('%Y%m%d%H')] = cur_freq_cols_useful

        if use_wet_criteria:
            lam_large_scale_power_ratio[
              cur_val_datetime.strftime('%Y%m%d%H')] = large_scale_power_ratio
            lam_wet_ratio[cur_val_datetime.strftime('%Y%m%d%H')] = wet_ratio

        data_number = data_number + 1

        cur_val_datetime = cur_val_datetime + timedelta(
            seconds=3600*analysis_interval)

    gridata = None
    return (fft_coef_useful, power_spectrum_useful,
            freq_rows_useful, freq_cols_useful, inv_fft_err,
            griddata_shape,
            lam_large_scale_power_ratio, lam_wet_ratio)


def obtain_bld_var_err_ratio(
        cur_model_var, total_model_levels, modeltype):
    err_ratio_bottom = BLD_VAR_ERR_RATIO[cur_model_var][modeltype][0]
    err_ratio_top = BLD_VAR_ERR_RATIO[cur_model_var][modeltype][1]

    err_ratio = numpy.linspace(
        err_ratio_bottom, err_ratio_top, total_model_levels)

    return err_ratio


def obtain_model_errors(model_error_matrix_over_frequencies,
                        lvl, fcst_length, griddata_shape,
                        cur_model_var,
                        glb_model_err_ratio, lam_model_err_ratio,
                        use_fix_wavenumber, fix_wavenumber,
                        global_model_err_ratio0,
                        global_model_err_ratio1):
    def _adjust_err(cur_lvl_fcst_power_spectrum,
                    cur_lvl_ana_power_spectrum,
                    res_thres,
                    err_ratio=None, err_add=None):
        """below resolution_threshold: large scale info
           above resolution_threshold: small scale info
          if use_hard_cutoff_err = False:
              anything lower than resolution_threshold would
              be very similar to global analysis,
              therefore later they will have very small errors
           """
        mm, nn = cur_lvl_fcst_power_spectrum.shape

        if err_ratio is not None:
            if res_thres > min(cur_lvl_fcst_power_spectrum.shape):
                cur_lvl_fcst_power_spectrum = \
                    cur_lvl_ana_power_spectrum * err_ratio
            else:
                cur_lvl_fcst_power_spectrum[0:res_thres, 0:res_thres] = \
                    cur_lvl_ana_power_spectrum[
                        0:res_thres,  0:res_thres]*err_ratio
                cur_lvl_fcst_power_spectrum[0:res_thres, nn-res_thres:nn] = \
                    cur_lvl_ana_power_spectrum[
                        0:res_thres, nn-res_thres:nn]*err_ratio
                cur_lvl_fcst_power_spectrum[mm-res_thres:mm,
                                            nn-res_thres:nn] = \
                    cur_lvl_ana_power_spectrum[
                        mm-res_thres:mm, nn-res_thres:nn]*err_ratio
                cur_lvl_fcst_power_spectrum[mm-res_thres:mm, 0:res_thres] = \
                    cur_lvl_ana_power_spectrum[
                        mm-res_thres:mm, 0:res_thres]*err_ratio

        if err_add is not None:
            if res_thres > min(griddata_shape):
                cur_lvl_fcst_power_spectrum = \
                    cur_lvl_ana_power_spectrum + err_add
            cur_lvl_fcst_power_spectrum[0:res_thres, 0:res_thres] = \
                cur_lvl_ana_power_spectrum[
                    0:res_thres, 0:res_thres] + err_add
            cur_lvl_fcst_power_spectrum[0:res_thres, nn-res_thres:nn] = \
                cur_lvl_ana_power_spectrum[
                    0:res_thres, nn-res_thres:nn] + err_add
            cur_lvl_fcst_power_spectrum[mm-res_thres:mm, nn-res_thres:nn] = \
                cur_lvl_ana_power_spectrum[
                    mm-res_thres:mm, nn-res_thres:nn] + err_add
            cur_lvl_fcst_power_spectrum[mm-res_thres:mm, 0:res_thres] = \
                cur_lvl_ana_power_spectrum[
                    mm-res_thres:mm, 0:res_thres] + err_add

        return cur_lvl_fcst_power_spectrum

    def _obtain_fcst_lengths_in_model_err(
            model_error_matrix_over_frequencies, fcst_length):
        """get the forecast length for the model error:
           * if it is provided by user: use the user's one
           * if not: use the first one in the model error matrix"""
        if not fcst_length:
            first_model_err_datetime = \
                model_error_matrix_over_frequencies.keys()[0]
            model_err_fcst_length = model_error_matrix_over_frequencies[
                first_model_err_datetime]['glb'].keys()[0]
        else:
            model_err_fcst_length = int(fcst_length)

        return model_err_fcst_length

    def _update_model_error(
            cur_glb_error, cur_lam_error,
            glb_model_err_ratio, lam_model_err_ratio):
        cur_glb_error = cur_glb_error * glb_model_err_ratio
        cur_lam_error = cur_lam_error * lam_model_err_ratio

        return cur_glb_error, cur_lam_error

    """-----------------------------
    # main function starts
    -----------------------------"""
    lvl_glb_fcst_err_list = []
    lvl_lam_fcst_err_list = []

    if use_fix_wavenumber:
        fix_wavenumber = int(fix_wavenumber)
        # if use the user defined error, reinitialize the power matrix
        cur_lvl_glb_ana_power_spectrum = numpy.zeros(griddata_shape)
        cur_lvl_glb_fcst_power_spectrum = numpy.ones(griddata_shape) * 999.0
        cur_lvl_lam_fcst_power_spectrum = numpy.ones(griddata_shape) * 999.0

        # update the lam power (for large scale information)
        cur_lvl_lam_fcst_power_spectrum = _adjust_err(
            cur_lvl_lam_fcst_power_spectrum,
            cur_lvl_glb_ana_power_spectrum,
            err_add=(1.0 - global_model_err_ratio0) * lam_model_err_ratio[lvl],
            res_thres=fix_wavenumber)

        # update the lam power err (the small scale information)
        cur_lvl_lam_fcst_power_spectrum[
            cur_lvl_lam_fcst_power_spectrum == 999.0] = \
            1.0 - global_model_err_ratio1

        # update the glb error (for large scale information)
        cur_lvl_glb_fcst_power_spectrum = _adjust_err(
            cur_lvl_glb_fcst_power_spectrum,
            cur_lvl_glb_ana_power_spectrum,
            err_add=global_model_err_ratio0 * glb_model_err_ratio[lvl],
            res_thres=fix_wavenumber)

        # update the global fcst err (for small scale information)
        cur_lvl_glb_fcst_power_spectrum[
            cur_lvl_glb_fcst_power_spectrum == 999.0] = global_model_err_ratio1

        cur_lvl_lam_fcst_power_spectrum = numpy.fft.fftshift(
            cur_lvl_lam_fcst_power_spectrum)
        cur_lvl_glb_fcst_power_spectrum = numpy.fft.fftshift(
            cur_lvl_glb_fcst_power_spectrum)

        cur_lvl_glb_ana_power_spectrum_vector = matrix2vector(
            cur_lvl_glb_ana_power_spectrum)
        cur_lvl_glb_fcst_power_spectrum_vector = matrix2vector(
            cur_lvl_glb_fcst_power_spectrum)
        cur_lvl_lam_fcst_power_spectrum_vector = matrix2vector(
            cur_lvl_lam_fcst_power_spectrum)

        cur_lvl_glb_fcst_err = numpy.asarray(
            cur_lvl_glb_fcst_power_spectrum_vector) - \
            numpy.asarray(cur_lvl_glb_ana_power_spectrum_vector)
        cur_lvl_lam_fcst_err = numpy.asarray(
            cur_lvl_lam_fcst_power_spectrum_vector) - \
            numpy.asarray(cur_lvl_glb_ana_power_spectrum_vector)

        lvl_glb_fcst_err_list.append(cur_lvl_glb_fcst_err)
        lvl_lam_fcst_err_list.append(cur_lvl_lam_fcst_err)

    else:
        for cur_datetime in model_error_matrix_over_frequencies:
            # get the forecast length for the model error:
            # * if it is provided by user: use the user's one
            # * if not: use the first one in the model error matrix
            model_err_fcst_length = _obtain_fcst_lengths_in_model_err(
                        model_error_matrix_over_frequencies, fcst_length)

            # convert model error matrix to vectors
            cur_glb_error = model_error_matrix_over_frequencies[
                    cur_datetime]['glb'][model_err_fcst_length][lvl]
            cur_lam_error = model_error_matrix_over_frequencies[
                    cur_datetime]['lam'][model_err_fcst_length][lvl]

            cur_glb_error, cur_lam_error = _update_model_error(
                cur_glb_error, cur_lam_error,
                glb_model_err_ratio[lvl], lam_model_err_ratio[lvl])

            cur_inv_fft_err_glb = numpy.asarray(matrix2vector(cur_glb_error))
            cur_inv_fft_err_lam = numpy.asarray(matrix2vector(cur_lam_error))

            cur_lvl_glb_fcst_err = cur_inv_fft_err_glb*cur_inv_fft_err_glb
            cur_lvl_lam_fcst_err = cur_inv_fft_err_lam*cur_inv_fft_err_lam

            lvl_glb_fcst_err_list.append(cur_lvl_glb_fcst_err)
            lvl_lam_fcst_err_list.append(cur_lvl_lam_fcst_err)

    lvl_glb_fcst_err = numpy.nanmean(
                numpy.asarray(lvl_glb_fcst_err_list), 0)
    lvl_lam_fcst_err = numpy.nanmean(
                numpy.asarray(lvl_lam_fcst_err_list), 0)

    lvl_glb_fcst_err = numpy.reshape(
        lvl_glb_fcst_err, (lvl_glb_fcst_err.shape[0], 1))
    lvl_lam_fcst_err = numpy.reshape(
        lvl_lam_fcst_err, (lvl_lam_fcst_err.shape[0], 1))

    return lvl_glb_fcst_err, lvl_lam_fcst_err


def setup_blue(all_err_coef, lvl_glb_fcst_err, lvl_lam_fcst_err):
    err_componenet0 = (lvl_glb_fcst_err + lvl_lam_fcst_err)**(-1)

    matrix_length = lvl_lam_fcst_err.shape[0]
    matrix_range = numpy.array(range(0, matrix_length))

    lvl_lam_fcst_err = csc_matrix(
        (lvl_lam_fcst_err.reshape(matrix_length),
         (matrix_range, matrix_range)),
        shape=(matrix_length, matrix_length))

    err_componenet0 = csc_matrix(
        (err_componenet0.reshape(matrix_length), (matrix_range, matrix_range)),
        shape=(matrix_length, matrix_length))

    err_coef = lvl_lam_fcst_err.dot(err_componenet0)

    all_err_coef.append(err_coef)

    return all_err_coef


def smooth_model_weights(model_weights,
                         total_model_levels, gaussian_sigma=1.0):
    """ using a gaussian filter to smooth the model weights"""
    # 1: extract the model weight diagnoal
    model_weights_diag = []
    for lvl in range(0, int(total_model_levels)):
        model_weights_diag.append(model_weights[lvl].diagonal())
    model_weights_diag_array = numpy.asarray(model_weights_diag)

    # 2: filter the model weight diagnoals vertically
    model_weights_diag_filtered = gaussian_filter1d(
        model_weights_diag_array, gaussian_sigma, axis=0)

    # 3: revert the diagnoals back to a sparse matrix
    matrix_length = model_weights_diag_array.shape[1]
    matrix_range = numpy.array(range(0, matrix_length))

    smoothed_model_weights = []
    for lvl in range(0, int(total_model_levels)):
        cur_model_weights = csc_matrix(
            (model_weights_diag_filtered[lvl, :],
             (matrix_range, matrix_range)),
            shape=(matrix_length, matrix_length))
        smoothed_model_weights.append(cur_model_weights)

    return smoothed_model_weights


def matrix2vector(matrix_in):
    vector_out = []
    for i in range(0, matrix_in.shape[0]):
        for j in range(0, matrix_in.shape[1]):
            vector_out.append(matrix_in[i, j])

    return vector_out


def save_error(error_data_path, data_err):
    fp = gzip.open(error_data_path, 'wb')
    cPickle.dump(data_err, fp)
    fp.close()


def save_wet_ratio(wet_ratio_path, lam_large_scale_power_ratio, lam_wet_ratio):
    fp = gzip.open(wet_ratio_path, 'wb')
    cPickle.dump(lam_large_scale_power_ratio, fp)
    cPickle.dump(lam_wet_ratio, fp)
    fp.close()
