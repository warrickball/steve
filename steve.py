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
                    "* 'query': just print what's available from MAST via lightkurve;"
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
                    help="oversample amplitude spectrum by this factor "
                    "(default=10.0)")
parser.add_argument("--psos", type=float, default=1.0,
                    help="oversample power spectrum by this factor "
                    "(default=1.0)")
parser.add_argument("--nyquist-factor", type=float, default=1.0,
                    help="Nyquist factor for Astropy's Lomb-Scargle method "
                    "(default=1)")
parser.add_argument("--smooth", type=int, default=0,
                    help="smooth power spectrum with top-hat this many bins wide")
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
parser.add_argument("--style", type=str, default=None,
                    help="load this Matplotlib style file")
parser.add_argument("--annotate", type=int, default=None,
                    help="mark the nth period: positive for maxima, "
                    "negative for minima")
parser.add_argument("--xlabel", type=str, nargs='+', default=None,
                    help="xlabel for lower x-axis")
parser.add_argument("--xlabel2", type=str, nargs='+', default=None,
                    help="xlabel for upper x-axis, if present")
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
parser.add_argument("-U", "--update", action='store_true', help="update cache")
parser.add_argument("--mission", type=str, default='TESS',
                    choices=['TESS', 'K2', 'kepler'])
parser.add_argument("--author", type=str, default='SPOC',
                    help="author argument for lightkurve's search "
                    "(default='SPOC')")
parser.add_argument("--exptime", type=int, default=120,
                    help="exposure time argument for lightkurve's search "
                    "(default=120)")
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

def jd_to_datetime(t):
    return time.Time(t, format='jd').to_datetime()

def PS():
    ppm = 1e6*(y/np.nanmedian(y)-1)
    f, p = LombScargle(0.0864*t, ppm).autopower(
        nyquist_factor=args.nyquist_factor,
        samples_per_peak=args.psos,
        normalization='psd')
    p = p*np.var(ppm)/np.trapz(p, x=f)

    if args.smooth > 0:
        p = np.convolve(p, np.ones(args.smooth)/args.smooth, mode='same')

    return f, p

def AS():
    f, p = LombScargle(t, T).autopower(
        normalization='psd',
        samples_per_peak=args.asos,
        nyquist_factor=args.nyquist_factor)
    I = np.argmax(p) + np.arange(-2,3) # 5 points
    f0 = curve_fit(lambda z, z0, A, w: A*np.sinc((z-z0)/w)**2,
                   f[I], p[I], (f[np.argmax(p)], np.max(p), 10*(f[2]-f[1])))[0][0]
    a = np.sqrt(p*4/len(t))
    P = -args.period/f0 if args.period < 0 else args.period

    return f, a, P

if not args.style is None:
    pl.style.use(args.style)

cache_file = '%s/%s-%s-%s-%i.npy' % (
    args.cache_dir, ' '.join(args.target).lower().replace(' ', '_'),
    args.mission.lower(), args.author.lower(), args.exptime)

if args.plot_kind == 'query':
    print(lk.search_lightcurve(' '.join(args.target), radius=args.radius))
    exit(0)

try:
    if args.no_cache or args.update:
        vprint("Not loading from cache: forcing failure... ", end='')
        raise FileNotFoundError("Not using cache: forcing failure...")
    else:
        vprint('Loading cache file %s... ' % cache_file, end='')
        data = np.load(cache_file)
except (FileNotFoundError, ValueError):
    vprint('Failed!\nDownloading data using lightkurve... ', end='')

    lc = lk.search_lightcurve(
        ' '.join(args.target),
        mission=args.mission,
        exptime=args.exptime,
        author=args.author,
        radius=args.radius).download_all().stitch()

    if len(lc) == 0:
        raise ValueError("no data found")

    data = np.zeros((6, len(lc)))

    data[0], data[-1] = lc.time.jd, lc.quality.value
    data[1] = lc.sap_flux.value
    data[2] = -2.5*np.log10(data[1]/np.nanmedian(data[1]))
    data[3] = lc.flux.value
    data[4] = -2.5*np.log10(data[3]/np.nanmedian(data[3]))

    if args.mission.lower() in ['tess']:
        data[[2,4]] += lc.meta['TESSMAG']
    elif args.mission.lower() in ['kepler', 'k2']:
        data[[2,4]] += lc.meta['KEPMAG']

    if not args.no_cache:
        vprint('Done.\nCaching data to file %s... ' % cache_file, end='')
        np.save(cache_file, data)

vprint('Done.')

if args.mission.lower() in ['tess']:
    mag_label = 'TESSmag'
elif args.mission.lower() in ['kepler', 'k2']:
    mag_label = 'Kp'
else:
    mag_label = 'mag'

t, q = data[[0,-1]]
if args.sap:
    y, T = data[[1,2]]
else:
    y, T = data[[3,4]]

I = (q == 0) & np.isfinite(t*y) & (t > args.t_min) & (t < args.t_max)
t, y, T = t[I], y[I], T[I]

if args.t0 is not None:
    t = t - np.nanmin(t) + args.t0

if args.plot_kind.lower() in ['lc', 'lightcurve', 'ts', 'timeseries']:
    if args.date:
        pl.plot_date(jd_to_datetime(t), T*args.scale_y, '.', alpha=args.alpha)
    else:
        pl.plot(t*args.scale_x, T*args.scale_y, '.', alpha=args.alpha)

    pl.gca().invert_yaxis()
    pl.ylabel(mag_label)
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
            pl.plot_date(list(map(jd_to_datetime, [t0, t1])), T0*np.ones(2), 'k-')
        else:
            pl.plot([t0, t1], T0*np.ones(2), 'k-')

        t0 = (t0+t1)/2.
        if args.date:
            t0 = jd_to_datetime(t0)

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
    tmod = t-t[np.argmax(T)] + args.phase_min*P

    pl.plot((tmod%P)*args.scale_x,
            (T + args.delta*np.floor(tmod/P))*args.scale_y,
            '.', alpha=args.alpha)
    pl.xlim(0.0, P)
    pl.gca().invert_yaxis()
    pl.xlabel('days mod %.2fd' % P)
    pl.ylabel(mag_label)

    ax2 = pl.gca().secondary_xaxis('top', functions=(lambda x: x/P, lambda x: x*P))
    ax2.set_xlabel('orbital phase')

elif args.plot_kind.lower() in ['peek']:
    from astropy.timeseries import LombScargle

    f, a, P = AS()
    vprint('Folding period = %.8f d' % P)

    pl.subplot(2,2,1)
    pl.plot(t, T, '.', alpha=args.alpha)
    pl.gca().invert_yaxis()
    pl.xlabel('days')
    pl.ylabel(mag_label)
    ax = pl.gca()
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')

    pl.subplot(2,2,2)
    pl.plot((t-t[np.argmin(y)] + args.phase_min*P)%P, T, '.', alpha=args.alpha)
    pl.gca().invert_yaxis()
    pl.xlabel('days mod %.2fd' % P)
    ax = pl.gca()
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position('right')
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')

    pl.subplot(2,2,3)
    pl.plot(f, a)
    pl.xlabel('frequency (c/d)')
    pl.xlim([-1, min(f.max(), 50.)+1])

    f, p = PS()
    pl.subplot(2,2,4)
    pl.loglog(f, p)
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

if args.xlabel2 is not None:
    ax2.set_xlabel(' '.join(args.xlabel2))

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
