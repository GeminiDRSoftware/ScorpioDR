import argparse
import os

from astropy.io import fits

#from astrodata import open as adopen
#import scorpio_instruments

phudict = {
    'common' : {
        'OBSMODE'  : 'spect',
        'GEMPRGID' : 'G-2019B-0207-Q',
        'CHANNEL'  : '%chan',
        'FILTER'   : '%chan',
        'FILTERN'  : None,
        'INPORT'   : 5,
        'IMTYPE'   : 'full',  # delete when Rodolfo confirms change
        'NEXTEND'  : 1,
        'SLITSIZE' : 4.32,
        'UTR*'     : None,
    },
    'bias' : {
        'OBJECT'   : 'Bias',
        'OBSMODE'  : 'image',
        'OBSTYPE'  : 'BIAS',
        'OBSCLASS' : 'dayCal',
        'OBSID'    : 'G-2019B-0207-Q-0014',
        'DATALAB'  : 'G-2019B-0207-Q-0014-0002-0001',
        'GCAL*'    : None,
        'SHUTTER'  : 'CLOSED',
    },
    'flat' : {
        'OBJECT'   : 'GCALflat',
        'OBSTYPE'  : 'FLAT',
        'OBSCLASS' : 'dayCal',
        'OBSID'    : 'G-2019B-0207-Q-0015',
        'DATALAB'  : 'G-2019B-0207-Q-0015-0002-0001',
        'GCALLAMP' : 'QH',
        'SHUTTER'  : 'CLOSED',
    },
    'arc' : {

    },
    'std' : {

    },
    'sci' : {

    }
}

extdict = {
    'common' : {
        'ARRNAM1'  : 'e2v_CCD231-84_%chan-1',
        'ARRNAM2'  : 'e2v_CCD231-84_%chan-2',
        'ARRNAM3'  : 'e2v_CCD231-84_%chan-3',
        'ARRNAM4'  : 'e2v_CCD231-84_%chan-4',
        'ARRNAMN'  : None,
        'ARRSEC1'  : '[1:2048,1545:2056]',
        'ARRSEC2'  : '[2049:4096,1545:2056]',
        'ARRSEC3'  : '[2049:4096,2057:2568]',
        'ARRSEC4'  : '[1:2048,2057:2568]',
        'ARRSECN'  : None,
        'CCD*'     : None,
        'CCDSUM'   : '4 1',
        'DATSEC1'  : '[11:522,11:522]',
        'DATSEC2'  : '[523:1034,11:522]',
        'DATSEC3'  : '[523:1034,523:1034]',
        'DATSEC4'  : '[11:522,523:1034]',
        'DATSECN'  : None,
        'DETSEC'   : '[1:4096,1545:2568]',
        'GAIN1'    : 1.7,
        'GAIN2'    : 1.7,
        'GAIN3'    : 1.7,
        'GAIN4'    : 1.7,
        'GAINN'    : None,
        'INHERIT'  : False,
        'OVRSECS1' : '[1:10,11:522]',
        'OVRSECS2' : '[1035:1044,11:522]',
        'OVRSECS3' : '[1035:1044,523:1034]',
        'OVRSECS4' : '[1:10,523:1034]',
        'OVRSECSN' : None,
        'OVRSECP1' : '[11:522,1:10]',
        'OVRSECP2' : '[523:1034,1:10]',
        'OVRSECP3' : '[523:1034,1035:1044]',
        'OVRSECP4' : '[11:522,1035:1044]',
        'OVRSECPN' : None,
        'RDNOIS1'  : 2.7,
        'RDNOIS2'  : 2.7,
        'RDNOIS3'  : 2.7,
        'RDNOIS4'  : 2.7,
        'RDNOISN'  : None,
        'RDOUTSPD' : 50000,
        'WAT*'     : None,
    },
    'bias' : {
    },
    'flat' : {
    },
    'arc' : {
    },
    'std' : {
    },
    'sci' : {
    }
}

parser = argparse.ArgumentParser(description="Update headers")

parser.add_argument("filenames", nargs="+", type=str, help="Files to update")

args = parser.parse_args()

for filename in args.filenames:

    froot, fext = os.path.splitext(filename)
    fbits = froot.split('.')
    funder = fbits[0].split('_')
    if fbits[0].endswith('_hdr'):  # ignore files already processed
        continue
    ftype = funder[funder.index('raw')+1]  # deduce obstype; not in header
    chan = fbits[-1]
    funder.append('hdr')
    fbits[0] = '_'.join(funder)
    fout = '.'.join(fbits) + fext

    hdulist = fits.open(filename)

    phu = hdulist[0].header
    exthdrs = [ext.header for ext in hdulist[1:]]

    for hdrs, kwdict in zip(([phu], exthdrs), (phudict, extdict)):
        for hdr in hdrs:
            for sect in ('common', ftype):
                for kw, val in kwdict[sect].items():
                    if val is None:
                        try:
                            del hdr[kw]
                        except KeyError:
                            pass
                    else:
                        if isinstance(val, str):
                            val = val.replace('%chan', chan)
                        hdr[kw] = val

    hdulist.writeto(fout, overwrite=True)

