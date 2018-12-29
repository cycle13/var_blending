#!/usr/bin/env python

import argparse
from datetime import datetime, timedelta
from var_blending import var_bld_err_processing
from var_blending import MODELTYPES
import os
import gc

"""
return parser.parse_args([
    '-s', '2018100300', '-e', '2018100300',
    '-i', '12', '-f', '24', '-t', '50', '-w', '/home/jzanetti/workspace/test/testdata',
    '-v', 'QVAPOR', #'U', 'V', 'QVAPOR', 'P',
    '-g', 's3://metservice-research-us-west-2/research/experiments/sijin/ncar_trip/blending/L12',
    '-a', 's3://metservice-research-us-west-2/research/experiments/sijin/ncar_trip/blending/L0',
    '-l', 's3://metservice-research-us-east-2/prod/archive-data/wrf_archive/wrfout',
    '-n', 'nz4kmN-NCEP', '--max_wavenumber', '15',
    '--use_inv_fft_as_err'])
    # '--use_fix_wavenumber', '-fw', '4', '-ge0', '0.0', '-ge1', '1.0'])

return parser.parse_args([
    '--download_model_from_s3',
    '-s', '201811211200', '-e', '201811211200',
    '-i', '12', '-f', '12', '-t', '50', '-w', '/home/jzanetti/Desktop/ncar/testdata/new',
    '-v', 'T', #'U', 'V', 'QVAPOR', 'P',
    '-g', 's3://metservice-research-us-west-2/research/experiments/sijin/ncar_trip/blending/L12',
    '-a', 's3://metservice-research-us-west-2/research/experiments/sijin/ncar_trip/blending/L0',
    '-l', 's3://metservice-research-us-east-2/prod/archive-data/wrf_archive/wrfout',
    '-d', '2', '--use_wet_criteria',
    '-n', 'nz4kmN-NCEP', '--max_wavenumber', '4',
    '--use_inv_fft_as_err'])
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
    parser.add_argument('-s', '--start_datetime', type=valid_datetime,
                        required=True, help="start analysis date/time")
    parser.add_argument('-e', '--end_datetime', type=valid_datetime,
                        required=True, help="end analysis date/time")
    parser.add_argument('-i', '--analysis_interval', type=str,
                        required=True, help="analysis interval")
    parser.add_argument('-f', '--forecast_length', type=str,
                        required=True, help="forecast lengths")
    parser.add_argument('-v', '--model_var_list', nargs='+',
                        required=True,
                        help="model_var_list, e.g., T, U, V, ...")
    parser.add_argument('-w', '--work_dir', type=str,
                        required=True, help="work directory")
    parser.add_argument('-t', '--total_model_levels', type=str, required=True,
                        help="total model levels")
    parser.add_argument('-d', '--domain_id', type=str, required=True,
                        help="domain_id")

    # -------------------------------
    # use rainfall as a criteria to determine whether to calculate errors
    # -------------------------------
    parser.add_argument('--use_wet_criteria', dest='use_wet_criteria',
                        default=False, help='use_wet_criteria',
                        action='store_true')

    # -------------------------------
    # global/lam data downloads
    # -------------------------------
    parser.add_argument('--download_model_from_s3',
                        dest='download_model_from_s3',
                        default=False, help='download_model_from_s3',
                        action='store_true')
    parser.add_argument('-g', '--global_fcst_on_s3', type=str,
                        required=False, help="global fcst on S3 to download")
    parser.add_argument('-a', '--global_ana_on_s3', type=str,
                        required=False,
                        help="global analysis on S3 to download")
    parser.add_argument('-l', '--lam_fcst_on_s3', type=str,
                        required=False, help="lam fcst on S3 to download")
    parser.add_argument('-n', '--lam_model_name', type=str,
                        required=False, help="lam_model_name")

    # -------------------------------
    # use dynamic wave number for global/lam
    # -------------------------------
    parser.add_argument('--use_inv_fft_as_err', dest='use_inv_fft_as_err',
                        default=False, help='use_inv_fft to compute errors',
                        action='store_true')
    parser.add_argument('-mw', '--max_wavenumber', type=str, required=False,
                        help="max_wavenumber from global model")

    return parser.parse_args()


def download_hist_data(args):
    start_val_datetime = args.start_datetime
    end_val_datetime = args.end_datetime
    cur_val_datetime = start_val_datetime

    while cur_val_datetime <= end_val_datetime:
        var_bld_err_processing.download_data_from_s3_to_local(
            args.work_dir, cur_val_datetime,
            args.global_fcst_on_s3, args.global_ana_on_s3,
            args.lam_fcst_on_s3, args.lam_model_name,
            args.domain_id,
            int(args.forecast_length))
        cur_val_datetime = cur_val_datetime + timedelta(
            seconds=3600*int(args.analysis_interval))


def check_hist_data(args):
    all_data_available = {}
    start_val_datetime = args.start_datetime
    end_val_datetime = args.end_datetime
    cur_val_datetime = start_val_datetime

    while cur_val_datetime <= end_val_datetime:
        all_data_available[cur_val_datetime.strftime('%Y%m%d%H')] = \
            var_bld_err_processing.init_all_data_available(
                int(args.forecast_length))

        # check forecasts
        for modeltype in MODELTYPES:
            model_path = os.path.join(
                    args.work_dir,  'hist_data',
                    cur_val_datetime.strftime('%Y%m%d%H'),
                    '{}_l{}_d0{}.nc'.format(
                        modeltype,
                        int(args.forecast_length),
                        args.domain_id))

            if not os.path.exists(model_path):
                all_data_available[cur_val_datetime.strftime('%Y%m%d%H')][
                    modeltype][int(args.forecast_length)] = False

        # check analysis
        model_path = os.path.join(
            args.work_dir, 'hist_data',
            cur_val_datetime.strftime('%Y%m%d%H'),
            'ana_l0_d0{}.nc'.format(args.domain_id))
        if not os.path.exists(model_path):
            all_data_available[
                cur_val_datetime.strftime('%Y%m%d%H')]['ana'][0] = False

        cur_val_datetime = cur_val_datetime + timedelta(
            seconds=3600*int(args.analysis_interval))

    return all_data_available


if __name__ == '__main__':
    args = setup_parser()
    # 1.0: download historical data for calculating model errors
    if args.download_model_from_s3:
        download_hist_data(args)

    # 2.0 check if all required data are available
    all_data_available = check_hist_data(args)

    # 3.0 get the model error
    lam_wet_ratio_saved = False
    for cur_model_var in args.model_var_list:
        # 3.1 init the variables for models
        (gridata, fft_coef_useful, power_spectrum_useful,
         freq_rows_useful, freq_cols_useful, data_number,
         inv_fft_err, griddata_shape,
         lam_wet_ratio, lam_large_scale_power_ratio) = \
                var_bld_err_processing.init_err_vars(args.use_inv_fft_as_err)

        # 3.2 get model data and error matrix over different frequencies
        (fft_coef_useful, power_spectrum_useful,
         freq_rows_useful, freq_cols_useful, inv_fft_err,
         griddata_shape, lam_large_scale_power_ratio,
         lam_wet_ratio) = var_bld_err_processing.extract_model_data_and_err(
            args.start_datetime, args.end_datetime,
            args.work_dir, cur_model_var,
            all_data_available,
            int(args.forecast_length), int(args.analysis_interval),
            gridata, data_number,
            int(args.max_wavenumber),
            fft_coef_useful, power_spectrum_useful,
            freq_rows_useful, freq_cols_useful, inv_fft_err,
            int(args.total_model_levels),
            griddata_shape, args.domain_id,
            lam_wet_ratio, lam_large_scale_power_ratio,
            adjust_lam_err_use_topo=False,
            use_wet_criteria=args.use_wet_criteria,
            use_inv_fft_as_err=args.use_inv_fft_as_err)

        # 3.3 save model error
        print 'saving model errors for {}'.format(cur_model_var)
        model_err_coef_path = os.path.join(
            args.work_dir, 'model_errs_{}.tar.gz'.format(cur_model_var))
        var_bld_err_processing.save_error(model_err_coef_path, inv_fft_err)

        # 3.4 save model powers
        print 'saving model powers for {}'.format(cur_model_var)
        model_powers_coef_path = os.path.join(
            args.work_dir, 'lam_powers_{}.tar.gz'.format(cur_model_var))
        var_bld_err_processing.save_error(
            model_powers_coef_path, power_spectrum_useful)

        # 3.4 save wet ratio
        if (not lam_wet_ratio_saved) and (args.use_wet_criteria):
            print 'saving LAM wet ratio'
            wet_ratio_path = os.path.join(
                args.work_dir, 'lam_wet_ratio.tar.gz')
            var_bld_err_processing.save_wet_ratio(
                wet_ratio_path, lam_large_scale_power_ratio, lam_wet_ratio)
            lam_wet_ratio_saved = True

        gc.collect()

    print 'done'
