#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as pl
import lightkurve as lk
from astropy import time
from astropy.timeseries import LombScargle
from scipy.optimize import curve_fit
from argparse import ArgumentParser

parser = ArgumentParser("Given a star name that `lightkurve` can resolve, "
                        "download the lightcurve data from MAST and plot the "
                        "lightcurve (or derived data).")
parser.add_argument("plot_kind", type=str, help="type of plot (case insensitive), one of "
                    "* 'lightcurve' or 'lc';"
                    "* 'fold', 'folded' or 'phased';"
                    "* 'ampspectrum' or 'as';"
                    "* 'ps', 'powerspectrum', 'ls' or 'lombscargle'; or"
                    "* 'peek': a multipanel plot combining 'lc', 'fold', 'as' and 'ps'")
parser.add_argument("target", type=str, nargs='+',
                    help="target identifier (e.g. gam Dor, HD27290)")
parser.add_argument("-P", "--period", type=float, default=-1.,
                    help="period on which to phase-fold lightcuve. "
                    "if negative, use that time period of the highest "
                    "peak in LS periodogram. e.g. -2 uses twice the "
                    "period of the tallest peak (default=-1.0)")
parser.add_argument("-d", "--delta", type=float, default=0,
                    help="displace each segment of phase-folded lightcurve by this much "
                    "(default=0)")
parser.add_argument("--phase-min", type=float, default=0.5,
                    help="phase of minimum brightness for phase-folded plots "
                    "(default=0.5)")
parser.add_argument("--asos", type=float, default=10.0,
                    help="overample amplitude spectrum by this factor "
                    "(default=10.0)")
parser.add_argument("--psos", type=float, default=1.0,
                    help="overample power spectrum by this factor "
                    "(default=1.0)")
parser.add_argument("--t0", type=float, default=None,
                    help="set first timestamp to this value "
                    "(default=don't change first timestamp)")
parser.add_argument("--t-min", type=float, default=-np.inf,
                    help="only use data with raw t > this "
                    "(default=-inf, i.e. use all data)")
parser.add_argument("--t-max", type=float, default=np.inf,
                    help="only use data with raw t < this "
                    "(default=inf, i.e. use all data)")
parser.add_argument("--date", action='store_true',
                    help="use actual date for lightcurve x-axis")
parser.add_argument("-a", "--axis", type=float, nargs=4, default=None,
                    help="axis parameters, as passed to matplotlib.axis")
parser.add_argument("-t", "--title", type=str, nargs='+', default=None,
                    help="specify title for plot, otherwise use target name")
parser.add_argument("--no-title", action='store_true',
                    help="don't show a title")
parser.add_argument("--annotate", type=int, default=None,
                    help="mark the nth period: positive for maxima, "
                    "negative for minima")
parser.add_argument("--xlabel", type=str, nargs='+', default=None,
                    help="xlabel for plot")
parser.add_argument("--ylabel", type=str, nargs='+', default=None,
                    help="ylabel for plot")
parser.add_argument("--scale-x", type=float, default=1.0,
                    help="multiply x-axis by this factor (default=1.0)")
parser.add_argument("--scale-y", type=float, default=1.0,
                    help="multiply y-axis by this factor (default=1.0)")
parser.add_argument("-o", "--output", type=str, default=None,
                    help="save plot to this file instead of showing")
parser.add_argument("--cache-dir", type=str, default='.',
                    help="cache data in this directory "
                    "(default=current directory)")
parser.add_argument("--no-cache", action='store_true', help="don't use cache")
parser.add_argument("--mission", type=str, default='TESS',
                    choices=['TESS', 'K2', 'kepler'])
parser.add_argument("--cadence", type=str, default='short')
parser.add_argument("-r", "--radius", type=float, default=None,
                    help="radius for lightkurve's search, in arcseconds "
                    "(default=None, i.e. lightkurve default)")
parser.add_argument("--sap", action='store_true',
                    help="use SAP flux instead of PDCSAP flux")
parser.add_argument("--alpha", type=float, default=1.0,
                    help="opacity for plotted points (default=1.0)")
parser.add_argument("-q", "--quiet", action='store_true')
args = parser.parse_args()

def vprint(*print_args, **print_kwargs):
    if not args.quiet:
        print(*print_args, **print_kwargs)

def mjd_to_datetime(mjd):
    return time.Time(mjd+2457000, format='jd').to_datetime()

def PS():
    ppm = 1e6*(y/np.nanmedian(y)-1)
    f, p = LombScargle(0.0864*t, ppm).autopower(
        nyquist_factor=1, samples_per_peak=args.psos,
        normalization='psd')
    p = p*np.var(ppm)/np.trapz(p, x=f)
    return f, p

def AS():
    f, p = LombScargle(t, T).autopower(
        normalization='psd', samples_per_peak=args.asos, nyquist_factor=1)
    I = np.argmax(p) + np.arange(-2,3) # 5 points
    f0 = curve_fit(lambda z, z0, A, w: A*np.sinc((z-z0)/w)**2,
                   f[I], p[I], (f[np.argmax(p)], np.max(p), 10*(f[2]-f[1])))[0][0]
    a = np.sqrt(p*4/len(t))
    P = -args.period/f0 if args.period < 0 else args.period

    return f, a, P

cache_file = '%s/%s.npy' % (
    args.cache_dir, ' '.join(args.target).lower().replace(' ', '_'))

try:
    if args.no_cache:
        vprint("Not using cache: forcing failure... ", end='')
        raise FileNotFoundError("Not using cache: forcing failure...")
    else:
        vprint('Loading cache file %s... ' % cache_file, end='')
        data = np.load(cache_file)
except (FileNotFoundError, ValueError):
    vprint('Failed!\nDownloading data using lightkurve... ', end='')

    lcs = lk.search_lightcurvefile(
        ' '.join(args.target), mission=args.mission, cadence=args.cadence, radius=args.radius).download_all()

    # no data
    if lcs is None:
        raise SystemExit(1)

    vprint('Done.\nCorrecting medians... ', end='')

    m_sap = np.nanmedian(lcs[0].SAP_FLUX.flux)
    m_pdcsap = np.nanmedian(lcs[0].PDCSAP_FLUX.flux)

    data = np.hstack([np.vstack([
        lc.PDCSAP_FLUX.time,
        lc.SAP_FLUX.flux-np.nanmedian(lc.SAP_FLUX.flux) + m_sap,
        lc.SAP_FLUX.flux_err,
        lc.PDCSAP_FLUX.flux-np.nanmedian(lc.PDCSAP_FLUX.flux) + m_pdcsap,
        lc.PDCSAP_FLUX.flux_err,
        lc.PDCSAP_FLUX.quality])
                      for lc in lcs])
    data = data[:,np.all(np.isfinite(data), axis=0)]

    vprint('Done.\nCaching data to file %s... ' % cache_file, end='')

    np.save(cache_file, data)

vprint('Done.')

if args.sap:
    t, y, dy, q = data[[0,1,2,5]]
else:
    t, y, dy, q = data[[0,3,4,5]]

I = (q == 0) & np.isfinite(t*y) & (t > args.t_min) & (t < args.t_max)

# hardcoded cleanup
if args.mission.lower() == 'tess':
    I = I & (np.abs(t-1348.35) > 1.05) # sector 1

t, y, dy = t[I], y[I], dy[I]
if args.mission.lower() in ['tess']:
    T = -2.5*np.log10(y) + 20.54
elif args.mission.lower() in ['kepler', 'k2']:
    # https://github.com/KeplerGO/lightkurve/issues/51
    T = 12 - 2.5*np.log10(y/1.74e5)

if args.t0 is not None:
    t = t - np.nanmin(t) + args.t_min

if args.plot_kind.lower() in ['lc', 'lightcurve', 'ts', 'timeseries']:
    if args.date:
        pl.plot_date(mjd_to_datetime(t), T*args.scale_y, '.', alpha=args.alpha)
    else:
        pl.plot(t*args.scale_x, T*args.scale_y, '.', alpha=args.alpha)

    pl.axis(np.array(pl.axis())[[0,1,3,2]])
    pl.ylabel('TESSmag')
    pl.xlabel('days')

    if args.annotate is not None:
        if args.period < 0:
            f, a, P = AS()
        else:
            P = args.period

        I = np.floor((t-t[0])/P).astype(int)
        if args.annotate > 0:
            T_extrema = np.array([np.min(T[i==I]) for i in np.unique(I)])
            t_extrema = np.array([t[i==I][np.argmin(T[i==I])] for i in np.unique(I)])
            T0 = np.mean(T_extrema[[abs(args.annotate)-1, abs(args.annotate)]]) - 0.01
            dT = 0.0
            va = 'bottom'
        else:
            T_extrema = np.array([np.max(T[i==I]) for i in np.unique(I)])
            t_extrema = np.array([t[i==I][np.argmax(T[i==I])] for i in np.unique(I)])
            T0 = np.mean(T_extrema[[abs(args.annotate)-1, abs(args.annotate)]]) + 0.01
            dT = 0.003
            va = 'top'

        t0, t1 = t_extrema[[abs(args.annotate)-1, abs(args.annotate)]]

        if args.date:
            pl.plot_date(list(map(mjd_to_datetime, [t0, t1])), T0*np.ones(2), 'k-')
        else:
            pl.plot([t0, t1], T0*np.ones(2), 'k-')

        t0 = (t0+t1)/2.
        if args.date:
            t0 = mjd_to_datetime(t0)

        pl.text(t0, T0+dT, 'P = %.3fd' % P, va=va, ha='center')

elif args.plot_kind.lower() in ['ps', 'powerspectrum', 'ls', 'lombscargle']:
    f, p = PS()
    pl.plot(f*args.scale_x, p*args.scale_y)
    pl.xlabel('frequency (µHz)')
    pl.ylabel('power density (ppm²/µHz)')

elif args.plot_kind.lower() in ['as', 'ampspectrum']:
    from astropy.timeseries import LombScargle
    f, a, P = AS()
    pl.plot(f*args.scale_x, a*args.scale_y)

    pl.xlabel('frequency (1/d)')
    pl.ylabel('amplitude (mag)')

elif args.plot_kind.lower() in ['fold', 'folded', 'phased']:
    from astropy.timeseries import LombScargle

    if args.period < 0:
        f, a, P = AS()
    else:
        P = args.period

    vprint('Folding period = %.8f d' % P)

    pl.plot(((t-t[np.argmax(T)] + args.phase_min*P)%P)*args.scale_x,
            (T + args.delta*np.floor((t-t[np.argmax(T)] + args.phase_min*P)/P))*args.scale_y,
            '.', alpha=args.alpha)
    pl.xlim(0.0, P)
    pl.axis(np.array(pl.axis())[[0,1,3,2]])
    pl.xlabel('days mod %.2fd' % P)
    pl.ylabel('TESSmag')

elif args.plot_kind.lower() in ['peek']:
    from astropy.timeseries import LombScargle

    f, a, P = AS()
    vprint('Folding period = %.8f d' % P)

    pl.subplot(2,2,1)
    pl.plot(t, T, '.', alpha=args.alpha)
    pl.axis(np.array(pl.axis())[[0,1,3,2]])
    pl.xlabel('days')
    pl.ylabel('TESSmag')
    ax = pl.gca()
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')

    pl.subplot(2,2,2)
    pl.plot((t-t[np.argmin(y)] + args.phase_min*P)%P, T, '.', alpha=args.alpha)
    pl.axis(np.array(pl.axis())[[0,1,3,2]])
    pl.xlabel('days mod %.2fd' % P)
    ax = pl.gca()
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position('right')
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')

    pl.subplot(2,2,3)
    pl.plot(f, a)
    pl.xlabel('frequency (c/d)')
    pl.xlim([-1., 51.])

    f, p = PS()
    pl.subplot(2,2,4)
    pl.loglog(f/0.0864, p)
    pl.xlabel('frequency (uHz)')
    ax = pl.gca()
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position('right')
else:
    print('%s is not a valid choice of plot' % args.plot_kind)

if args.axis:
    pl.axis(args.axis)

if args.date:
    pl.gcf().autofmt_xdate()
    pl.xlabel('date')

if args.xlabel is not None:
    pl.xlabel(' '.join(args.xlabel))

if args.ylabel is not None:
    pl.ylabel(' '.join(args.ylabel))

if not args.no_title:
    if args.title is None:
        pl.suptitle(' '.join(args.target))
    else:
        pl.suptitle(' '.join(args.title))

if args.output:
    pl.savefig(args.output)
else:
    pl.show()
