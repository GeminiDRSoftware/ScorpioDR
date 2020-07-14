
# If the filters name and central wavelength are different from 
# the global definitions in gemini_instruments/gemini/lookup.py
# redefine them here in filter_wavelengths.
 
# 20200514 - Not sure if we need these defined yet
filter_wavelengths = {
    #'Ks' : 2.1570,
    'g'   : 0.469,
    'r'   : 0.622, 
    'i'   : 0.755, 
    'z'   : 0.889, 
    'Y'   : 1.040, 
    'J'   : 1.250, 
    'H'   : 1.630, 
    'Ks'  : 2.175,
}

# 20200514 - For testing
array_properties = {
    "gain"  :   1.7,    # electrons/ADU
    "read_noise"    :   4.0,      # electrons
}
