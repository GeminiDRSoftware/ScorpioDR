#!/bin/env python

# A hacky script to help James convert the original 4D simulated data into the
# expected new 2D format (kept in case they need regenerating for some reason).

import argparse
import os
from datetime import timedelta

from astropy.io import fits

import astrodata
import gemini_instruments
import scorpio_instruments
from astrodata import open as adopen

def copyhdr(kw, refhdr, newhdr):
    if kw in refhdr:
        newhdr[kw] = (refhdr[kw], refhdr.comments[kw])
    elif kw in newhdr:
        try:
            del newhdr[kw]
        except KeyError:
            pass

parser = argparse.ArgumentParser(description="Convert 4D sims to 2D")

parser.add_argument("inlist", type=str, help="File list to convert")

args = parser.parse_args()

refs, data, outfn, dts = [], [], [], []

with open(args.inlist) as f_inlist:
    for f in f_inlist:

        f = f.strip()
        path, fn = os.path.split(f)
        fbase, fext = os.path.splitext(fn)

        refs.append(fits.open(f))
        if 'CCD' in refs[-1][0].header['DETECTOR']:
            arm = 'vis'
            data.append(fits.open(f))
        else:
            arm = 'nir'
            data.append(fits.open(fbase + '_slopeDetermined.fits'))
        outfn.append(fbase)
        dts.append(adopen(f).ut_datetime())

for expnum, (ref, dat, fn, dt) in enumerate(zip(refs, data, outfn, dts),
                                            start=1):

    ref_phu, dat_phu = ref[0].header, dat[0].header

    for kw in ['ADUTOELE']:
        copyhdr(kw, dat_phu, ref_phu)
    if arm == 'nir':
        ref_phu['LINCORR'] = dat_phu['ADDVAR']  # skipped (no coeffs)

    ref_exth, dat_exth = ref['SCI'].header, dat['SCI'].header

    for kw in ['SATLEVEL', 'NONLINEA', 'BSCALE', 'BZERO', 'BUNIT']:
        copyhdr(kw, dat_exth, ref_exth)

    if 'ADUTOELE' in dat_phu:
        copyhdr('GAIN', dat_exth, ref_exth)
        for n in range(1, 17):
            ref_exth[f'GAINOR{n}'] = (dat_exth[f'GAIN{n}'], f'Original GAIN{n}')
            try:
                # del ref_exth[f'GAIN{n}']
                del ref_exth[f'GAIN']
            except KeyError:
                pass
            ref_exth[f'GAIN{n}'] = 1.0

    ref_phu['EXPTIME'] = ref_exth['INTTIME']  # split to total exp into 1 int
    ref_phu['EXPTREQ'] = ref_exth['INTTREQ']

    tdelt = timedelta(seconds=ref_exth.get('INTTOTT' or 0.))

    # The sims provided by GWU have both the wrong GEMPRGID format for GPP and
    # also identical dummy values for every exposure/integration, so add some
    # usable values:
    keys = ref_phu['OBJECT'], ref_phu.get('GCALSHUT', None)
    obsnum = (('Bias', None), ('Dark', None), ('GCALflat', 'CLOSED'),
              ('GCALflat', 'OPEN'), ('Object', None)).index(keys) + 1
    if '_CCD' not in ref_phu['DETECTOR']:
        obsnum += 10

    try:
        del ref_phu['GEMPGRID']  # orig sims have typo
    except KeyError:
        pass
    ref_phu['GEMPRGID'] = 'G-2019B-1234-Q'
    ref_phu['OBSID'] = ref_phu['GEMPRGID'] + f'-{obsnum:04d}'

    sci_arr = dat['SCI'].data
    if len(sci_arr.shape) > 3:  # raw ccd
        sci_arr = sci_arr[:, 0, :, :]  # bin unit read group axis
    elif len(sci_arr.shape) < 3:  # cals
        sci_arr = [sci_arr]

    for n, plane in enumerate(sci_arr, start=1):

        # Fix this to deal with dimensionality of cals:
        ref['SCI'].data = plane

        ref_phu['DATALAB'] = ref_phu['OBSID'] + f'-{expnum:04d}-{n:04d}'
        ref_phu['DATE-OBS'] = dt.date().isoformat()
        ref_phu['TIME-OBS'] = dt.time().isoformat()

        # ref.write(fn, overwrite=True)
        ref.writeto(fn + f'-{n}.fits', output_verify='warn', overwrite=True)

        dt += tdelt

