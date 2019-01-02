from scripts import run_var_blending

# TEST 1: 
#    (1) using the estimated model errors; 
#    (2) using the GFS and WRF valid at the same time
test_data = {
    'bld_fcst': 'testdata/bld_d02.nc',
    'work_dir': 'testdata/working',
    'model_err_dir': 'testdata',
    'total_model_levels': '50',
    'model_var_list': 'T',
    'max_wavenumber': '4',
    'use_wet_power_constrain_from_hist': True,
    'smooth_vertical_model_weights': True
    }


# TEST 2: 
#    (1) using the estimated model errors; 
#    (2) using the GFS and WRF valid at different time
'''
test_data = {
    'work_dir': '/home/jzanetti/workspace/test/testdata/new',
    'var_blending': {
        'wrf': '/home/jzanetti/Desktop/ncar/testdata/wrf_fcst/wrf_hourly_nz4kmN-NCEP_d02_2018-11-06_12_00_00',
        'gfs': '/home/jzanetti/Desktop/ncar/testdata/global_fcst/global_a2018110900_v2018110912.nc'}, 
    }
'''

test_para = [
    #'-l', test_data['lam_fcst'],
    #'-g', test_data['glb_fcst'],
    '-b', test_data['bld_fcst'],
    '-w', test_data['work_dir'],
    '-v', test_data['model_var_list'],
    '-t', test_data['total_model_levels'],
    '-mw', test_data['max_wavenumber'],
    ]

optional_para = ['use_wet_power_constrain_from_hist', 
                 'smooth_vertical_model_weights', 
                 'generate_plot']

for citem in optional_para:
    try:
        if test_data[citem]:
            test_para.append('--{}'.format(citem))
    except KeyError:
        pass

'''
run_var_blending.run_var_blending(
    {'para': test_para,
     'mdbz_path': 'testdata//cur_mdbz.data',
     'topo_path': 'testdata/cur_topo.data',
     'cur_lam': 'testdata/cur_lam_{var}.data',
     'cur_glb': 'testdata/cur_glb_{var}.data'})
'''
    
run_var_blending.run_var_blending(
    {'para': test_para,
     'mdbz_path': 'testdata//cur_mdbz.data',
     'topo_path': 'testdata/cur_topo.data',
     'cur_lam': '/home/jzanetti/Downloads/wrf_hourly_nz4kmN-NCEP-var-bld-fix_d02_2018-11-24_21_00_00',
     'cur_glb': '/home/jzanetti/Desktop/ncar/testdata/global_fcst/global_a2018110900_v2018110912.nc'})    




