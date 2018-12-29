import numpy
from scipy import fftpack, signal
from copy import deepcopy
from var_blending import BLD_VAR_ERR_RATIO
from wrf import smooth2d


def obtain_power_spectrum_using_welch(
        cur_gridata, fcst_length, modeltypes,
        total_model_levels,
        max_wavenumber):

    # init power_spectrum_welch dict
    power_spectrum_welch = {}
    for modeltype in modeltypes:
        power_spectrum_welch[modeltype] = {}
        power_spectrum_welch[modeltype][fcst_length] = []

    # obtain power spectrum along lat
    data_shape = cur_gridata['lam'][12].shape
    for cur_lvl in range(0, total_model_levels):
        lam_power_along_lat = []
        glb_power_along_lat = []
        lam_power_along_lon = []
        glb_power_along_lon = []
        for i in range(0, data_shape[1]):
            lam = cur_gridata['lam'][fcst_length][cur_lvl, i, :]
            glb = cur_gridata['glb'][fcst_length][cur_lvl, i, :]
            lam_power_along_lat.append(
                signal.welch(lam, nperseg=data_shape[2])[1])
            glb_power_along_lat.append(
                signal.welch(glb, nperseg=data_shape[2])[1])

        for i in range(0, data_shape[2]):
            lam = cur_gridata['lam'][fcst_length][cur_lvl, :, i]
            glb = cur_gridata['glb'][fcst_length][cur_lvl, :, i]
            lam_power_along_lon.append(
                signal.welch(lam, nperseg=data_shape[1])[1])
            glb_power_along_lon.append(
                signal.welch(glb, nperseg=data_shape[1])[1])

        lam_power_along_lat = numpy.asarray(lam_power_along_lat)
        glb_power_along_lat = numpy.asarray(glb_power_along_lat)

        lam_power_along_lat_mean = numpy.mean(
            lam_power_along_lat, 0)[0: max_wavenumber]
        glb_power_along_lat_mean = numpy.mean(
            glb_power_along_lat, 0)[0: max_wavenumber]

        lam_power_along_lon = numpy.asarray(lam_power_along_lon)
        glb_power_along_lon = numpy.asarray(glb_power_along_lon)

        lam_power_along_lon_mean = numpy.mean(
            lam_power_along_lon, 0)[0: max_wavenumber]
        glb_power_along_lon_mean = numpy.mean(
            glb_power_along_lon, 0)[0: max_wavenumber]

        power_spectrum_welch['lam'][fcst_length].append(
            numpy.mean(
                (lam_power_along_lat_mean + lam_power_along_lon_mean)/2.0))
        power_spectrum_welch['glb'][fcst_length].append(
            numpy.mean(
                (glb_power_along_lat_mean + glb_power_along_lon_mean)/2.0))

    power_spectrum_welch['lam'][fcst_length] = numpy.asarray(
        power_spectrum_welch['lam'][fcst_length])
    power_spectrum_welch['glb'][fcst_length] = numpy.asarray(
        power_spectrum_welch['glb'][fcst_length])

    return power_spectrum_welch


def run_fft(griddata):
    def _obtain_fft_power(fft_coeff_all):
        power_spectrum_all = (numpy.abs(fft_coeff_all)**2)**(0.5)
        power_spectrum_useful = power_spectrum_all
        fft_coef_useful = fft_coeff_all

        freq_rows_useful = numpy.fft.fftfreq(fft_coeff_all.shape[0])
        freq_cols_useful = numpy.fft.fftfreq(fft_coeff_all.shape[1])

        return (fft_coef_useful, power_spectrum_useful,
                freq_rows_useful, freq_cols_useful)

    fft_coeff_all = numpy.fft.fft2(griddata)
    fft_coeff_all = fftpack.fftshift(fft_coeff_all)

    (fft_coef_useful, power_spectrum_useful,
     freq_rows_useful, freq_cols_useful) = _obtain_fft_power(fft_coeff_all)

    # https://www.researchgate.net/post/How_to_find_
    # Power_Spectral_Density_PSD_of_any_image
    power_spectrum_useful = abs(
        numpy.fft.fftshift(numpy.fft.fft2(griddata)))**2

    return (fft_coef_useful, power_spectrum_useful,
            freq_rows_useful, freq_cols_useful)


def obtain_fftdata(cur_gridata,
                   modeltypes, fcst_length,
                   start_level, end_level):
    def _init_out_dict(modeltypes, fcst_lengt):
        fft_coef_useful = {}
        power_spectrum_useful = {}
        freq_rows_useful = {}
        freq_cols_useful = {}
        for modeltype in modeltypes:
            fft_coef_useful[modeltype] = {}
            power_spectrum_useful[modeltype] = {}
            freq_rows_useful[modeltype] = {}
            freq_cols_useful[modeltype] = {}

            fft_coef_useful[modeltype][fcst_length] = []
            power_spectrum_useful[modeltype][fcst_length] = []
            freq_rows_useful[modeltype][fcst_length] = []
            freq_cols_useful[modeltype][fcst_length] = []

        fft_coef_useful['ana'] = {}
        power_spectrum_useful['ana'] = {}
        freq_rows_useful['ana'] = {}
        freq_cols_useful['ana'] = {}

        fft_coef_useful['ana'][0] = []
        power_spectrum_useful['ana'][0] = []
        freq_rows_useful['ana'][0] = []
        freq_cols_useful['ana'][0] = []

        return (fft_coef_useful, power_spectrum_useful,
                freq_rows_useful, freq_cols_useful)

    (cur_fft_coef_useful, cur_power_spectrum_useful,
     cur_freq_rows_useful, cur_freq_cols_useful) = _init_out_dict(
        modeltypes, fcst_length)

    for modeltype in modeltypes:
        cur_griddata = cur_gridata[modeltype][fcst_length]
        for i in range(start_level, end_level):
            cur_griddata_lvl = cur_griddata[i, :, :]
            (cur_fft_coef_useful_lvl, cur_power_spectrum_useful_lvl,
             cur_freq_rows_useful_lvl, cur_freq_cols_useful_lvl) = run_fft(
                cur_griddata_lvl)

            cur_fft_coef_useful[modeltype][fcst_length].append(
                cur_fft_coef_useful_lvl)
            cur_power_spectrum_useful[modeltype][fcst_length].append(
                cur_power_spectrum_useful_lvl)
            cur_freq_rows_useful[modeltype][fcst_length].append(
                cur_freq_rows_useful_lvl)
            cur_freq_cols_useful[modeltype][fcst_length].append(
                cur_freq_cols_useful_lvl)

    cur_griddata = cur_gridata['ana'][0]
    for i in range(start_level, end_level):
        cur_griddata_lvl = cur_griddata[i, :, :]
        (cur_fft_coef_useful_lvl, cur_power_spectrum_useful_lvl,
         cur_freq_rows_useful_lvl, cur_freq_cols_useful_lvl) = run_fft(
            cur_griddata_lvl)

        cur_fft_coef_useful['ana'][0].append(
            cur_fft_coef_useful_lvl)
        cur_power_spectrum_useful['ana'][0].append(
            cur_power_spectrum_useful_lvl)
        cur_freq_rows_useful['ana'][0].append(
            cur_freq_rows_useful_lvl)
        cur_freq_cols_useful['ana'][0].append(
            cur_freq_cols_useful_lvl)

    return (cur_fft_coef_useful, cur_power_spectrum_useful,
            cur_freq_rows_useful, cur_freq_cols_useful)


def obtain_inv_fft_rms(fft_coef_useful, fcst_length,
                       model_variable,
                       total_model_levels,
                       max_wavenumber_threshold=None,
                       fft_coef_wind_spd=None,
                       smooth_err_over_freq=True):
    def _run_ifft(fft_coef, coef_threshold):
        data_shape = fft_coef.shape
        fft_wave0 = numpy.zeros(data_shape, dtype='complex')
        fft_wave1 = numpy.zeros(data_shape, dtype='complex')

        fft_wave0[:, data_shape[1]/2-coef_threshold:data_shape[1]/2+coef_threshold, :] = \
            fft_coef[:, data_shape[1]/2-coef_threshold:data_shape[1]/2+coef_threshold, :]
        fft_wave1[:, :, data_shape[2]/2-coef_threshold:data_shape[2]/2+coef_threshold] = \
            fft_coef[:, :, data_shape[2]/2-coef_threshold:data_shape[2]/2+coef_threshold]

        fft_wave = fft_wave0 + fft_wave1
        xx = numpy.fft.ifftshift(fft_wave, axes=1)
        xx = numpy.fft.ifftshift(xx, axes=2)

        rev_fft = numpy.fft.ifftn(xx, axes=(1, 2))

        return rev_fft

    def _create_err_array(err_dist1, cur_coef_thres, pre_coef_thres, cur_rms):
        data_shape = err_dist1[0, :, :].shape
        # err_dist1 = numpy.zeros(err_dist.shape)
        err_dist2 = deepcopy(err_dist1)
        for i in range(0, err_dist1.shape[0]):
            err_dist1[i, data_shape[0]/2-cur_coef_thres: data_shape[0]/2-pre_coef_thres, :] = cur_rms[i]
            err_dist1[i, data_shape[0]/2+pre_coef_thres: data_shape[0]/2+cur_coef_thres, :] = cur_rms[i]
            err_dist1[i, :, data_shape[1]/2-cur_coef_thres: data_shape[1]/2-pre_coef_thres] = cur_rms[i]
            err_dist1[i, :, data_shape[1]/2+pre_coef_thres: data_shape[1]/2+cur_coef_thres] = cur_rms[i]
        x = err_dist2 > 0.0
        err_dist1[x] = err_dist2[x]

        return err_dist1

    def _create_mean_profile(cur_lam, cur_glb, cur_ana):
        cur_lam_mean_profile = numpy.mean(numpy.mean(cur_lam.real, 1), 1)
        cur_glb_mean_profile = numpy.mean(numpy.mean(cur_glb.real, 1), 1)
        cur_ana_mean_profile = numpy.mean(numpy.mean(cur_ana.real, 1), 1)

        return (cur_lam_mean_profile, cur_glb_mean_profile,
                cur_ana_mean_profile)

    def _obtain_user_external_err_ratio(model_variable):
        # user can set a ratio to time the calculated the errors
        # default: for both global and LAM are both 1.0
        try:
            if model_variable in ['U', 'V']:
                bld_var_err_lam_ratio = BLD_VAR_ERR_RATIO['UV']['lam']
                bld_var_err_glb_ratio = BLD_VAR_ERR_RATIO['UV']['glb']
            else:
                bld_var_err_lam_ratio = BLD_VAR_ERR_RATIO[
                    model_variable]['lam']
                bld_var_err_glb_ratio = BLD_VAR_ERR_RATIO[
                    model_variable]['glb']
        except KeyError:
            bld_var_err_lam_ratio = 1.0
            bld_var_err_glb_ratio = 1.0

        return bld_var_err_lam_ratio, bld_var_err_glb_ratio

    def _obtain_model_errors(cur_lam, cur_glb, cur_ana,
                             bld_var_err_lam_ratio, bld_var_err_glb_ratio,
                             cur_lam_mean_profile, cur_glb_mean_profile,
                             cur_ana_mean_profile,
                             total_model_levels):
        """obtain the model errors for LAM and global models"""

        # step 1: the errors directly calculated
        #         between LAM/global and analysis
        lam_err = numpy.sqrt((cur_lam.real - cur_ana.real)**2)
        glb_err = numpy.sqrt((cur_glb.real - cur_ana.real)**2)

        lam_err_profile_mean = numpy.mean(numpy.mean(
            lam_err, 1), 1)
        glb_err_profile_mean = numpy.mean(numpy.mean(
            glb_err, 1), 1)

        # step 2: update error if the analysis is not in the
        #         middle of LAM and global models
        for cur_lvl in range(0, total_model_levels):
            cur_lam_mean = cur_lam_mean_profile[cur_lvl]
            cur_glb_mean = cur_glb_mean_profile[cur_lvl]
            cur_ana_mean = cur_ana_mean_profile[cur_lvl]

            # step 2.1: if analysis is not in the middle
            #           of LAM and global model
            if not ((cur_lam_mean < cur_ana_mean < cur_glb_mean) or
                    (cur_glb_mean < cur_ana_mean < cur_ana_mean)):
                # step 2.2: if LAM has smaller error than the global model
                if (abs(cur_lam_mean - cur_ana_mean) <
                        abs(cur_glb_mean - cur_ana_mean)):
                    lam_err_profile_mean[cur_lvl] = 0.0

                # step 2.2: if LAM has bigger error than the global model
                if (abs(cur_lam_mean - cur_ana_mean) >
                        abs(cur_glb_mean - cur_ana_mean)):
                    glb_err_profile_mean[cur_lvl] = 0.0

        return lam_err_profile_mean, glb_err_profile_mean

    def _smooth_model_error_over_frequencies(
            lam_err, glb_err, total_model_levels,
            smooth_number=3):
        """smoothing error matrix over frequencies"""
        for lvl in range(0, total_model_levels):
            lam_err[lvl, :, :] = smooth2d(
                deepcopy(lam_err[lvl, :, :]),
                smooth_number)
            glb_err[lvl, :, :] = smooth2d(
                deepcopy(glb_err[lvl, :, :]),
                smooth_number)

        return lam_err, glb_err

    """-----------------------------
    # main function starts
    -----------------------------"""
    inv_fft_err = {'lam': {}, 'glb': {}}
    lam_err = numpy.zeros(numpy.asarray(
        fft_coef_useful['lam'][fcst_length]).shape)
    glb_err = numpy.zeros(numpy.asarray(
        fft_coef_useful['lam'][fcst_length]).shape)
    pre_glb_err_profile = None
    pre_lam_err_profile = None
    pre_thres = 0
    if not max_wavenumber_threshold:
        max_wavenumber_threshold = min(
            numpy.asarray(fft_coef_useful['lam'][fcst_length]).shape[1],
            numpy.asarray(fft_coef_useful['lam'][fcst_length]).shape[2])

    for istep, cur_thres in enumerate(range(1, max_wavenumber_threshold)):
        if cur_thres > max_wavenumber_threshold:
            break

        if not fft_coef_wind_spd:
            cur_lam = _run_ifft(numpy.asarray(
                fft_coef_useful['lam'][fcst_length]), cur_thres)
            cur_glb = _run_ifft(numpy.asarray(
                fft_coef_useful['glb'][fcst_length]), cur_thres)
            cur_ana = _run_ifft(numpy.asarray(
                fft_coef_useful['ana'][0]), cur_thres)
        else:
            cur_lam = _run_ifft(numpy.asarray(
                fft_coef_wind_spd['lam'][fcst_length]), cur_thres)
            cur_glb = _run_ifft(numpy.asarray(
                fft_coef_wind_spd['glb'][fcst_length]), cur_thres)
            cur_ana = _run_ifft(numpy.asarray(
                fft_coef_wind_spd['ana'][0]), cur_thres)

        # user can set a ratio to time the calculated the errors
        # default: for both global and LAM are both 1.0
        (bld_var_err_lam_ratio, bld_var_err_glb_ratio) = \
            _obtain_user_external_err_ratio(model_variable)

        (cur_lam_mean_profile, cur_glb_mean_profile,
         cur_ana_mean_profile) = _create_mean_profile(
            cur_lam, cur_glb, cur_ana)

        (cur_lam_err_profile, cur_glb_err_profile) = \
            _obtain_model_errors(
                     cur_lam, cur_glb, cur_ana,
                     bld_var_err_lam_ratio, bld_var_err_glb_ratio,
                     cur_lam_mean_profile, cur_glb_mean_profile,
                     cur_ana_mean_profile,
                     total_model_levels)

        if istep >= 1:
            glb_err_diff_profile = cur_glb_err_profile - pre_glb_err_profile
            lam_err_diff_profile = cur_lam_err_profile - pre_lam_err_profile
            err_diff_profile = numpy.abs(glb_err_diff_profile -
                                         lam_err_diff_profile)
            # set a very big error for glb
            err_diff_profile[err_diff_profile < 0.001] = 999.0
            if numpy.all(err_diff_profile == 999.0):
                print ' -- all levels has absolute bias smaller than 0.001'
                break

        if istep > 1:
            lam_err = _create_err_array(
                lam_err, cur_thres, pre_thres,
                (cur_lam_err_profile-pre_lam_err_profile))
            glb_err = _create_err_array(
                glb_err, cur_thres, pre_thres,
                (cur_glb_err_profile-pre_glb_err_profile))
        else:
            lam_err = _create_err_array(
                lam_err, cur_thres,
                pre_thres, cur_lam_err_profile)
            glb_err = _create_err_array(
                glb_err, cur_thres,
                pre_thres, cur_glb_err_profile)

        pre_thres = cur_thres
        pre_glb_err_profile = cur_glb_err_profile
        pre_lam_err_profile = cur_lam_err_profile

    lam_err[lam_err == 0.0] = 0.0
    glb_err[glb_err == 0.0] = 1.0

    if smooth_err_over_freq:
        lam_err, glb_err = _smooth_model_error_over_frequencies(
            lam_err, glb_err, total_model_levels,
            smooth_number=3)

    inv_fft_err['lam'][fcst_length] = lam_err
    inv_fft_err['glb'][fcst_length] = glb_err

    return inv_fft_err
