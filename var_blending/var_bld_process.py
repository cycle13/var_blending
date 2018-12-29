from shutil import copyfile
import numpy
import os
from netCDF4 import Dataset
import gzip
import cPickle
import fft_process
from var_blending import MODELTYPES

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def obtain_cur_power_diff(
        lam_data, glb_data, max_wavenumber, total_model_levels):
    cur_gdata_dict = {}
    cur_gdata_dict['lam'] = {}
    cur_gdata_dict['glb'] = {}
    cur_gdata_dict['lam'][12] = lam_data
    cur_gdata_dict['glb'][12] = glb_data
    power_spectrum_welch_lam = fft_process.obtain_power_spectrum_using_welch(
        cur_gdata_dict, 12, MODELTYPES,
        total_model_levels, max_wavenumber)
    cur_power_diff = abs(power_spectrum_welch_lam['lam'][12] -
                         power_spectrum_welch_lam['glb'][12])
    return cur_power_diff


def obtain_hist_power_diff(work_dir, cur_model_var):
    """obtain the historical:
        * wet_ratio: rainfall area precentage over lands
                     [obtained later from constants]
        * lam_large_scale_power_ratio_list: rainfall power ratio
          for the wavenumber less than max_wavenumber
          [obtained later from constants]
        * power_diff: minimum power difference between
                      global and LAM over levels
    These 3 variables are used to calculate model errors
    """
    # get power_diff
    cur_hist_power_path = os.path.join(
            work_dir, 'lam_powers_{}.tar.gz'.format(cur_model_var))
    lam_power_list = []
    glb_power_list = []
    cur_hist_power = load_power(cur_hist_power_path)
    for cur_datetime, _ in cur_hist_power.iteritems():
        lam_power_list.append(cur_hist_power[cur_datetime]['lam'][12])
        glb_power_list.append(cur_hist_power[cur_datetime]['glb'][12])

    lam_power = numpy.asarray(lam_power_list)
    glb_power = numpy.asarray(glb_power_list)

    power_diff = abs(lam_power - glb_power)
    power_diff = numpy.min(power_diff, 0)

    return power_diff


def copy_netcdf_file(nc_in, nc_out):
    copyfile(nc_in, nc_out)


def load_error(error_data_path):
    fp = gzip.open(error_data_path, 'rb')
    model_err = cPickle.load(fp)
    fp.close()
    return model_err


def load_power(power_data_path):
    fp = gzip.open(power_data_path, 'rb')
    model_power = cPickle.load(fp)
    fp.close()
    return model_power


def load_wet_ratio(wet_ratio_path):
    fp = gzip.open(wet_ratio_path, 'rb')
    lam_large_scale_power_ratio = cPickle.load(fp)
    lam_wet_ratio = cPickle.load(fp)
    fp.close()
    return lam_large_scale_power_ratio, lam_wet_ratio


def freq2grid(fft_coef_vector, grid_shape):
    cur_analysis_fft_coef_reshaped = fft_coef_vector.reshape(grid_shape)
    cur_analysis_reshiftted = numpy.fft.ifftshift(
        cur_analysis_fft_coef_reshaped)
    cur_analysis = numpy.fft.ifft2(cur_analysis_reshiftted)

    return cur_analysis.real


def plot_analyzed_data(work_dir, model_var, lam_fcst,
                       glb_fcst, ana_fcst, lvl,
                       only_analysis=False):
    max_value = max(numpy.max(lam_fcst),
                    numpy.max(glb_fcst), numpy.max(ana_fcst))
    min_value = max(numpy.min(lam_fcst),
                    numpy.min(glb_fcst), numpy.min(ana_fcst))

    work_plot_dir = os.path.join(work_dir, 'figure', model_var)
    if not os.path.exists(work_plot_dir):
        os.makedirs(work_plot_dir)

    if only_analysis:
        plt.figure(figsize=(8.5, 10))
        plt.pcolormesh(ana_fcst)
        plt.colorbar()
        plt.clim(min_value, max_value)
        plt.title('analyzed fcst, lvl: {}'.format(lvl))
        plt.grid()
        plt.savefig(
            "analysis_results_L{}.png".format(lvl), bbox_inches='tight')
        plt.close()
    else:
        plt.suptitle('{}, level: {}'.format(model_var, lvl))
        plt.figure(figsize=(15, 15))
        plt.subplot(231)
        plt.pcolormesh(lam_fcst)
        plt.colorbar()
        plt.clim(min_value, max_value)
        plt.title('LAM fcst')
        plt.grid()

        plt.subplot(232)
        plt.pcolormesh(glb_fcst)
        plt.colorbar()
        plt.clim(min_value, max_value)
        plt.title('GLB fcst')
        plt.grid()

        plt.subplot(233)
        plt.pcolormesh(ana_fcst)
        plt.colorbar()
        plt.clim(min_value, max_value)
        plt.title('analyzed fcst')
        plt.grid()

        diff0 = ana_fcst - lam_fcst
        diff1 = ana_fcst - glb_fcst
        diff_max0 = max(abs(numpy.max(diff0)), abs(numpy.min(diff0)))
        diff_max1 = max(abs(numpy.max(diff1)), abs(numpy.min(diff1)))
        diff_max = max(diff_max0, diff_max1)

        plt.subplot(234)
        plt.pcolormesh(diff0, cmap=plt.get_cmap('seismic'))
        plt.clim(-diff_max, diff_max)
        plt.colorbar()
        plt.title('analysis increment (ana - lam)')
        plt.grid()

        plt.subplot(235)
        plt.pcolormesh(diff1, cmap=plt.get_cmap('seismic'))
        plt.clim(-diff_max, diff_max)
        plt.colorbar()
        plt.title('analysis increment (ana - glb)')
        plt.grid()

        plt.savefig(os.path.join(work_plot_dir,
                    "all_results_L{}.png".format(lvl)), bbox_inches='tight')
        plt.close()


def plot_amb(analysis_minus_lam, analysis_minus_glb,
             model_var, work_dir):
    plt.title(model_var)
    plt.plot(analysis_minus_lam, range(0, len(analysis_minus_lam)),
             'b', label='RMS: LAM ~ ANA')
    plt.plot(analysis_minus_glb, range(0, len(analysis_minus_lam)),
             'r', label='RMS: GLB ~ ANA')
    plt.legend()
    plt.savefig(os.path.join(work_dir, "{}_results.png".format(model_var)),
                bbox_inches='tight')
    plt.close()


def write_netcdf_out(tmp_nc, bld_fcst, model_var_list, analysis_data):
    dset = Dataset(tmp_nc, 'r+')
    for cur_model_var in model_var_list:
        try:
            dset.variables[cur_model_var][:] = numpy.asarray(
                analysis_data[cur_model_var]).reshape(
                (dset[cur_model_var][:].shape))
        except KeyError:
            pass
    dset.close()
    copyfile(tmp_nc, bld_fcst)
    os.remove(tmp_nc)
