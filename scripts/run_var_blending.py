#!/usr/bin/env python

import argparse
from datetime import datetime
import numpy

from var_blending import var_bld_err_processing, var_bld_process
from var_blending import fft_process
from var_blending import (PLOT_LEVS, MODEL_WEIGHTS_SMOOTHING_SIGMA)
import os
import sys

"""
return parser.parse_args([
    '-l', '/home/jzanetti/Desktop/ncar/testdata/new/hist_data/2018112112/lam_l12_d02.nc',
    '-g', '/home/jzanetti/Desktop/ncar/testdata/new/hist_data/2018112112/glb_l12_d02.nc',
    #'-l', '/home/jzanetti/Downloads/wrf_hourly_nz4kmN-NCEP-var-bld-fix_d02_2018-11-24_21_00_00',
    #'-l', '/home/jzanetti/Desktop/ncar/testdata/wrf_fcst/wrf_hourly_nz4kmN-NCEP_d02_2018-11-06_12_00_00',
    #'-g', '/home/jzanetti/Desktop/ncar/testdata/global_fcst/global_a2018110900_v2018110912.nc',
    #'-l', '/home/jzanetti/workspace/test/testdata/hist_data/2018101600/lam_l24.nc',
    #'-g', '/home/jzanetti/workspace/test/testdata/hist_data/2018101600/glb_l24.nc',
    '-b', '/home/jzanetti/workspace/test/testdata/zsj.nc',
    '-v', 'T', #'U', 'V', 'QVAPOR', 'P',
    '-t', '50',
    '--smooth_vertical_model_weights',
    '--use_wet_power_constrain_from_hist',
    '--generate_plot',
    '-mw', '4',
    '-w', '/home/jzanetti/workspace/test/testdata/new',
    #'--use_fix_wavenumber', '-fw', '4', '-ge0', '0.3', '-ge1', '0.7',
    ])
"""


def valid_datetime(timestamp):
    '''turn a timestamp into a datetime object'''
    try:
        return datetime.strptime(timestamp, "%Y%m%d%H%M")
    except ValueError:
        msg = "Not a valid date: '{}'.".format(timestamp)

    raise argparse.ArgumentTypeError(msg)


def setup_parser():
    parser = argparse.ArgumentParser(
        description='create errors for VAR blending')
    parser.add_argument('-l', '--lam_fcst', type=str,
                        required=True, help="lam_fcst (on s3 or local)")
    parser.add_argument('-g', '--glb_fcst', type=str,
                        required=True, help="glb_fcst (on s3 or local)")
    parser.add_argument('-b', '--bld_fcst', type=str,
                        required=True, help="bld_fcst (output)")
    parser.add_argument('-t', '--total_model_levels', type=str,
                        required=True, help="total_model_levels")
    parser.add_argument('-v', '--model_var_list', nargs='+',
                        required=True,
                        help="model_var_list, e.g., T, U, V, ...")
    parser.add_argument('-w', '--work_dir', type=str,
                        required=True, help="work directory")
    parser.add_argument('-mw', '--max_wavenumber',
                        type=str, required=True,
                        help="max_wavenumber from global model")

    parser.add_argument('--use_wet_power_constrain_from_hist',
                        dest='use_wet_power_constrain_from_hist',
                        default=False,
                        help='use wet power constrain from the '
                             'data used for calculating BE',
                        action='store_true')

    # ----------------------
    # optional: generate plots
    # ----------------------
    parser.add_argument('--generate_plot', dest='generate_plot',
                        default=False, help='generate_plot',
                        action='store_true')

    # ----------------------
    # optional: forecast length in model errors
    # ----------------------
    parser.add_argument('-f', '--forecast_length', type=str,
                        required=False, default=None,
                        help=('if forecast_length is not defined, '
                              'the first datatime '
                              'in model_error matrix will be used'))
    # ----------------------
    # optional: whether to smooth the vertical model weights
    # ----------------------
    parser.add_argument('--smooth_vertical_model_weights',
                        dest='smooth_vertical_model_weights',
                        default=False, help='smooth vertical model weights',
                        action='store_true')

    # ----------------------
    # optional: forecast length in model errors
    # ----------------------
    parser.add_argument('--use_fix_wavenumber', dest='use_fix_wavenumber',
                        default=False,
                        help='use_fix_wavenumber to compute errors',
                        action='store_true')
    parser.add_argument('-fw', '--fix_wavenumber', type=str, required=False,
                        default=None, help="fix_wavenumber")

    return parser.parse_args()


if __name__ == '__main__':
    args = setup_parser()

    # wet criteria check:
    # if the rainfall area over lands is too big, or
    #    the rainfall power ratio over large scale is too big:
    # we assume that the large scale forcing in LAM is big enough, so
    # we are not doing the blending
    if args.use_wet_power_constrain_from_hist:
        print '    step 0: wet criteria checks'
        passed_wet_criteria, large_scale_power_ratio, wet_ratio = \
            var_bld_err_processing.wet_criteria(
                args.lam_fcst, int(args.max_wavenumber))
        if not passed_wet_criteria:
            print('    wet criteria check failed: '
                  'large_scale_power_ratio: {}, wet_ratio: {}'.format(
                    large_scale_power_ratio, wet_ratio))
            var_bld_process.copy_netcdf_file(args.lam_fcst, args.bld_fcst)
            sys.exit()

    var_bld_process.copy_netcdf_file(args.lam_fcst, 'tmp.nc')

    analysis_data = {}
    for cur_model_var in args.model_var_list:
        print 'current model variable: {}'.format(cur_model_var)

        print '    step 1: load model fcsts to be analyzed'
        lam_data, _, _ = var_bld_err_processing.obtain_griddata(
            args.lam_fcst, cur_model_var)
        glb_data, _, _ = var_bld_err_processing.obtain_griddata(
            args.glb_fcst, cur_model_var)

        # * load historical power difference and
        #   calculate the current power difference
        # * if the current power difference is smaller than the
        #   smallest historical difference (usually at 12h fcst), that
        #   particular level is not included in the analysis
        if args.use_wet_power_constrain_from_hist:
            print '    step 1.1: load historical model power'
            hist_power_diff = var_bld_process.obtain_hist_power_diff(
                args.work_dir, cur_model_var)
            cur_power_diff = var_bld_process.obtain_cur_power_diff(
                lam_data, glb_data, int(args.max_wavenumber),
                int(args.total_model_levels))

        print '    step 2: load model errors'
        if not args.use_fix_wavenumber:
            cur_err_path = os.path.join(
                args.work_dir, 'model_errs_{}.tar.gz'.format(cur_model_var))
            cur_err = var_bld_process.load_error(cur_err_path)
        else:
            cur_err = None

        print '    step 3: obtain model error ratio'
        glb_model_err_ratio = var_bld_err_processing.obtain_bld_var_err_ratio(
                cur_model_var, int(args.total_model_levels), 'glb')
        lam_model_err_ratio = var_bld_err_processing.obtain_bld_var_err_ratio(
                cur_model_var, int(args.total_model_levels), 'lam')

        print '    step 4: convert model error from matrix to ' + \
              'vector, and calculated model weights'
        model_weights = []
        for lvl in range(0, int(args.total_model_levels)):
            # 3.1 convert model error matrix to vectors
            (lvl_glb_fcst_err, lvl_lam_fcst_err) = \
             var_bld_err_processing.obtain_model_errors(
                cur_err, lvl, args.forecast_length,
                (lam_data.shape[1], lam_data.shape[2]),
                cur_model_var,
                glb_model_err_ratio, lam_model_err_ratio,
                args.use_fix_wavenumber,
                args.fix_wavenumber,
                float(args.global_model_err_ratio0),
                float(args.global_model_err_ratio1))

            # 3.2 calculating the model weights
            model_weights = var_bld_err_processing.setup_blue(
                model_weights, lvl_glb_fcst_err, lvl_lam_fcst_err)

        # smooth vertical errors between LAM and global data
        if args.smooth_vertical_model_weights:
            print '    step 5: smoothing the calculated model weights'
            model_weights = var_bld_err_processing.smooth_model_weights(
                model_weights, int(args.total_model_levels),
                gaussian_sigma=MODEL_WEIGHTS_SMOOTHING_SIGMA)

        print '    step 6: running analysis'
        analysis_minus_lam = []
        analysis_minus_glb = []
        analysis_data_list = []
        for lvl in range(0, int(args.total_model_levels)):
            cur_lvl_lam_data = lam_data[lvl, :, :]
            cur_lvl_glb_data = glb_data[lvl, :, :]
            (cur_lam_fft_coef, power_spectrum_lam, _, _) = \
                fft_process.run_fft(cur_lvl_lam_data)
            (cur_glb_fft_coef, power_spectrum_glb, _, _) = \
                fft_process.run_fft(cur_lvl_glb_data)
            cur_lam_fft_coef_vector = numpy.asarray(
                var_bld_err_processing.matrix2vector(cur_lam_fft_coef))
            cur_glb_fft_coef_vector = numpy.asarray(
                var_bld_err_processing.matrix2vector(cur_glb_fft_coef))
            cur_model_weight = model_weights[lvl]

            if args.use_wet_power_constrain_from_hist and \
                    (cur_power_diff[lvl] < hist_power_diff[lvl]):
                print('    * skipped level {}: cur_power_diff ({}) '
                      '< hist_power_diff ({})'.format(
                        lvl, cur_power_diff[lvl], hist_power_diff[lvl]))
                cur_analysis = cur_lvl_lam_data
            else:
                cur_analysis_fft_coef = cur_model_weight*(
                    cur_glb_fft_coef_vector - cur_lam_fft_coef_vector) + \
                    cur_lam_fft_coef_vector

                cur_analysis = var_bld_process.freq2grid(
                    cur_analysis_fft_coef, cur_lvl_lam_data.shape)
            if args.generate_plot:
                if lvl in PLOT_LEVS:
                    var_bld_process.plot_analyzed_data(
                            args.work_dir, cur_model_var,
                            cur_lvl_lam_data, cur_lvl_glb_data,
                            cur_analysis, lvl, only_analysis=False)

            analysis_minus_lam.append(numpy.mean(
                numpy.sqrt((cur_analysis - cur_lvl_lam_data)**2)))
            analysis_minus_glb.append(
                numpy.mean(numpy.sqrt((cur_analysis - cur_lvl_glb_data)**2)))
            analysis_data_list.append(cur_analysis)

        analysis_data[cur_model_var] = analysis_data_list

        if args.generate_plot:
            var_bld_process.plot_amb(
                analysis_minus_lam, analysis_minus_glb,
                cur_model_var, args.work_dir)

    print 'saving model analysis'
    var_bld_process.write_netcdf_out(
        'tmp.nc', args.bld_fcst, args.model_var_list, analysis_data)

    print 'done'
