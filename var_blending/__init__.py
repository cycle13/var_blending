"""var_blending package

See the README file in subversion .
"""

MODELTYPES = ['lam', 'glb']

BLD_VAR_ERR_RATIO = {
    # e.g., [1.0, 1.0] means the error magnifier
    #        ratios at [bottom_level, top_level]
    # if a fixed wavenumber is used, this section only
    #        affects the large scale spectrum
    'T': {'lam': [1.0, 1.0], 'glb': [1.0, 1.0]},
    # 'T': {'lam': [1.0, 1.0], 'glb': [3.0, 0.75]},
    'UV': {'lam': [1.0, 1.0], 'glb': [1.0, 1.0]},
    'QVAPOR': {'lam': [1.0, 1.0], 'glb': [1.0, 1.0]}}

PLOT_LEVS = [0]
# PLOT_LEVS = [0, 5, 10, 15, 25,
#             30, 35, 40, 45, 50]

MODEL_WEIGHTS_SMOOTHING_SIGMA = 3.0

# * if the power of LAM in large scale is larger
#     than 40% of the total scale, give up
# * if the rainfall area of LAM over lands is over 30%, give up
WET_CRITERIAS = {
    'max_allowed_large_scale_power': 0.4,
    'max_allowed_wet_ratio': 0.3
    }
