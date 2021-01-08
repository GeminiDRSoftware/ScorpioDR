
# If the filters name and central wavelength are different from 
# the global definitions in gemini_instruments/gemini/lookup.py
# redefine them here in filter_wavelengths.

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

# Added 20200827 - Landon Gelman
amplifier_count = {
    "G"   :   4,
    "R"   :   4,
    "I"   :   4,
    "Z"   :   4,
    "Y"   :   32,
    "J"   :   32,
    "H"   :   32,
    "K"   :   32,
}

amplifier_order = {
    "G" :   [1,2,4,3],
    "R" :   [1,2,4,3],
    "I" :   [1,2,4,3],
    "Z" :   [1,2,4,3],
    "Y" :   [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32],
    "J" :   [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32],
    "H" :   [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32],
    "K" :   [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32],
}

# 20200615 - For testing, put near-average values in
array_properties = {
    "gainG"  :   1.7,    # electrons/ADU
    "gainR"  :   1.7,
    "gainI"  :   1.7,
    "gainZ"  :   1.7,
    "gainY"  :   1.7,
    "gainJ"  :   1.7,
    "gainH"  :   1.7,
    "gainK"  :   1.7,
    "read_noiseG"    :   2.63,      # electrons
    "read_noiseR"    :   2.95,
    "read_noiseI"    :   2.50,
    "read_noiseZ"    :   2.62,
    "read_noiseY"    :   5.0,
    "read_noiseJ"    :   5.3,
    "read_noiseH"    :   4.6,
    "read_noiseK"    :   4.4,
}
