# -*- coding: utf-8 -*-
"""
QDIV-gui configuration parameters

Created on Sat Jul 11 07:59:24 2026

@author: marchalj
"""
version="11.07.2026"

import numpy as np
import detectors
detectors.show()
    
debug_mode=False
chatty=False
path = './data/'
filename = 'data'
det=detectors.CSPEC_module
nb_of_pixels_per_tube=1024   
GUI_refresh_period=0.25
hardware="R5560"
firmware="cspec_rmm 2026.03.11"
IP=["10.128.0.50","10.128.0.51"]
registers_def_file='RegisterFile_cspec_rmm_2026_03_11.json'

#register_filename='registers_file_TAU_1p1us_PZ_1p1us_5_4.txt'
register_filename='registers_file_ArCO2.txt'

nb_of_bins_in_spectrum=2**9


