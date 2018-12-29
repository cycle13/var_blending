import os
import errno
from datetime import timedelta
import subprocess
from AMPSAws import resources


def download_lam_fcst(args, lam_fcst):
    """download LAM forecasts from S3"""
    lam_fcst_on_s3 = args.lam_fcst_on_s3_for_download

    cur_analysis_datetime = args.start_analysis_date

    while cur_analysis_datetime <= args.end_analysis_date:
        cur_lam_fcst_on_s3_dir = os.path.join(
            lam_fcst_on_s3, args.model_name,
            cur_analysis_datetime.strftime('%y/%m/%d/%H'))

        cur_lam_fcst_dir = os.path.join(
            lam_fcst, cur_analysis_datetime.strftime('%Y%m%d%H'))

        if not os.path.join(cur_lam_fcst_dir):
            os.makedirs(cur_lam_fcst_dir)

        cur_fcst_filename = ('wrf_hourly_{model_name}_'
                             'd0{domain_id}_{model_time}').format(
            model_name=args.model_name, domain_id=args.domain_id,
            model_time=(cur_analysis_datetime +
                        timedelta(seconds=3600*int(
                            args.fcst_length))).strftime('%Y-%m-%d_%H:%M:%S'))

        cur_remote_fcst = os.path.join(
            cur_lam_fcst_on_s3_dir, cur_fcst_filename)
        cur_local_fcst = os.path.join(cur_lam_fcst_dir, cur_fcst_filename)

        if not os.path.exists(cur_local_fcst):
            try:
                resources.copy(cur_remote_fcst, cur_local_fcst)
            except:
                print 'download failed: check the remote at {}'.format(
                        cur_remote_fcst)
                cur_analysis_datetime = cur_analysis_datetime + timedelta(
                    seconds=int(args.analysis_interval)*3600)
                pass

        cur_analysis_datetime = cur_analysis_datetime + timedelta(
            seconds=int(args.analysis_interval)*3600)


def download_global_analysis(args, global_analysis):
    """download global analysis/forecasts"""
    global_analysis_on_s3 = args.global_analysis_on_s3_for_download

    cur_analysis_datetime = args.start_analysis_date

    while cur_analysis_datetime <= args.end_analysis_date:
        cur_global_analysis_filename = \
            'global_a{cur_analysis}_v{cur_valid}.nc'.format(
                cur_analysis=cur_analysis_datetime.strftime('%Y%m%d%H'),
                cur_valid=(cur_analysis_datetime +
                           timedelta(
                            seconds=3600*int(args.fcst_length))).strftime(
                                '%Y%m%d%H'))
        cur_remote_fcst = os.path.join(global_analysis_on_s3,
                                       cur_global_analysis_filename)
        cur_local_fcst = os.path.join(global_analysis,
                                      cur_global_analysis_filename)

        if not os.path.exists(cur_local_fcst):
            try:
                resources.copy(cur_remote_fcst, cur_local_fcst)
            except:
                cur_analysis_datetime = cur_analysis_datetime + timedelta(
                    seconds=int(args.analysis_interval)*3600)
                pass

        cur_analysis_datetime = cur_analysis_datetime + timedelta(
            seconds=int(args.analysis_interval)*3600)


def download_global_bias(args, gmbias_dir):
    """download bias files for the global models"""
    remote_gmbias = os.path.join(
        args.global_bias_on_s3, 'L{}'.format(args.fcst_length),
        'gm_bias.GFS.15kmMIDUS.asc')
    remote_warm_nl = os.path.join(
        args.global_bias_on_s3, 'L{}'.format(args.fcst_length),
        'warm.nml')
    local_gmbias = os.path.join(gmbias_dir, 'gm_bias.GFS.15kmMIDUS.asc')
    local_warm_nl = os.path.join(gmbias_dir, 'warm.nml')

    resources.copy(remote_gmbias, local_gmbias)
    resources.copy(remote_warm_nl, local_warm_nl)


def run_blending(args, lam_fcst, global_analysis, gmbias, blending):
    """run blending"""
    cur_analysis_datetime = args.start_analysis_date

    while cur_analysis_datetime <= args.end_analysis_date:
        cur_blending = os.path.join(
            blending, cur_analysis_datetime.strftime('%Y%m%d%H'))
        if not os.path.exists(cur_blending):
            os.makedirs(cur_blending)

        # link lam forecast
        src_lam_fcst = os.path.join(
            lam_fcst, cur_analysis_datetime.strftime('%Y%m%d%H'),
            'wrf_hourly_{model_name}_d0{domain_id}_{model_time}'.format(
                model_name=args.model_name, domain_id=args.domain_id,
                model_time=(
                    cur_analysis_datetime +
                    timedelta(seconds=3600*int(args.fcst_length))).strftime(
                                '%Y-%m-%d_%H:%M:%S')))
        dst_lam_fcst = os.path.join(
            cur_blending, 'lam.nc')
        symlink_force(src_lam_fcst, dst_lam_fcst)

        # link global fcst
        src_global_analysis = os.path.join(
            global_analysis, 'global_a{cur_analysis}_v{cur_valid}.nc'.format(
                cur_analysis=cur_analysis_datetime.strftime('%Y%m%d%H'),
                cur_valid=(cur_analysis_datetime +
                           timedelta(seconds=3600*int(
                            args.fcst_length))).strftime('%Y%m%d%H')))
        dst_global_analysis1 = os.path.join(cur_blending, 'gm.t0.nc')
        dst_global_analysis2 = os.path.join(cur_blending, 'gm.t1.nc')
        symlink_force(src_global_analysis, dst_global_analysis1)
        symlink_force(src_global_analysis, dst_global_analysis2)

        # link gmbias and warm.nl
        src_gmbias = os.path.join(gmbias, 'gm_bias.GFS.15kmMIDUS.asc')
        src_nl = os.path.join(gmbias, 'warm.nml')
        dest_gmbias = os.path.join(cur_blending, 'gm_bias.GFS.15kmMIDUS.asc')
        dest_nl = os.path.join(cur_blending, 'warm.nml')
        symlink_force(src_gmbias, dest_gmbias)
        symlink_force(src_nl, dest_nl)

        # link warm.exe
        src_warm_exe = os.path.join(args.warm_dir, 'warm.exe')
        dst_warm_exe = os.path.join(cur_blending, 'warm.exe')
        symlink_force(src_warm_exe, dst_warm_exe)

        startprogram('./warm.exe', cur_blending, shell=True)

        cur_analysis_datetime = cur_analysis_datetime + timedelta(
            seconds=int(args.analysis_interval)*3600)


def upload_processed_lam(args, blending):
    """upload blended data back to S3"""
    cur_analysis_datetime = args.start_analysis_date

    while cur_analysis_datetime <= args.end_analysis_date:
        cur_blending = os.path.join(
            blending, cur_analysis_datetime.strftime('%Y%m%d%H'))
        if not os.path.exists(cur_blending):
            os.makedirs(cur_blending)
        cur_processed_lam_fcst = os.path.join(cur_blending, 'test.nc')
        remote_processed_lam_dir = os.path.join(
            args.lam_fcst_on_s3_for_upload, args.model_name,
            cur_analysis_datetime.strftime('%y/%m/%d/%H'))
        remote_processed_lam_filename = \
            'wrf_hourly_{model_name}_d0{domain_id}_{model_time}'.format(
                model_name=args.model_name, domain_id=args.domain_id,
                model_time=(cur_analysis_datetime + timedelta(
                    seconds=3600*int(args.fcst_length))).strftime(
                        '%Y-%m-%d_%H:%M:%S'))
        remote_processed_lam = os.path.join(remote_processed_lam_dir,
                                            remote_processed_lam_filename)

        try:
            resources.copy(cur_processed_lam_fcst, remote_processed_lam)
            print 'from {} to {}'.format(
                cur_processed_lam_fcst, remote_processed_lam)
        except:
            cur_analysis_datetime = cur_analysis_datetime + timedelta(
                seconds=int(args.analysis_interval)*3600)
            pass

        cur_analysis_datetime = cur_analysis_datetime + timedelta(
            seconds=int(args.analysis_interval)*3600)


def symlink_force(src, dest):
    try:
        os.symlink(src, dest)
    except OSError, e:
        if e.errno == errno.EEXIST:
            os.remove(dest)
            os.symlink(src, dest)
        else:
            raise e


def startprogram(cmd, cmd_dir, shell=False):
    process = subprocess.Popen(cmd, cwd=cmd_dir, shell=shell)
    output, err = process.communicate()
