# steve

**steve** is the **ste**llar **v**ariability **e**xplorer: a small
Python script that plots lightcurves downloaded from the [Mikulski
Archive for Space Telescopes (MAST)](http://archive.stsci.edu/) using
[`lightkurve`](http://docs.lightkurve.org/).  It's an inelegant,
incomplete and unstable tool for quickly looking at a target's
variability.

## Requirements

`steve` requires the Python packages
[NumPy](https://numpy.org/),
[SciPy](https://scipy.org/),
[Matplotlib](https://matplotlib.org/),
[Astropy](https://www.astropy.org/) and
[Lightkurve](http://docs.lightkurve.org/).
All can be (probably) be installed however you install Python packages.

## Usage

To use `steve`, clone this repo and run `python3 steve.py`.  You can
see the various command line options by running `python3 steve.py -h`.
To get a feel for how `steve` works, let's try the bright eclipsing
binary [ζ Phe](https://en.wikipedia.org/wiki/Zeta_Phoenicis).  `steve`
can make five types of plot:

* a basic lightcurve (`lc` or `lightcurve`),
* a phase-folded lightcurve (`fold`, `folded` or `phased`), 
* an amplitude spectrum (`ampspectrum` or `as`) 
* a power spectrum (`powerspectrum` or `ps`) or
* all of the above on one four-panel plot (`peek`).

`steve`'s basic syntax is `python3 steve.py <plot type> <target
name>`.  So try

    python3 steve.py lc zet phe

You should get the (beautiful!) TESS lightcurve for ζ Phe.  Also
notice that wherever you ran `steve`, you now should have a file
called `zet_phe.npy`. This a crude cache file so that `steve` doesn't
have to keep querying MAST.¹ You can control where these are put using
the `--cache` option.

There are loads of eclipses, so let's try phase-folding the
lightcurve. We can also add some transparency to the points with
`--alpha`. Try

    python3 steve.py fold zet phe --alpha 0.1

`steve` fits a sinc function to the tallest peak in the amplitude
spectrum to estimate the principal period² and has clearly found half
the true period (which is a common problem). You can control the
phase-folding period with the `-P` or `--period`.  We could specify
the known period using `-P 1.66977` or we can use twice `steve`'s
estimate with `-P -2`, because negative values mean "multiply
`steve`'s estimate by this". So try

    python3 steve.py fold zet phe --alpha 0.1 -P -2

To see why `steve` gets the period wrong, try plotting the amplitude
spectrum with

    python3 steve.py as zet phe

You can also plot the power spectrum but that's mostly for people
(like me) interested in solar-like oscillations.³

Finally, it's sometimes easier to look at everything at once, so
there's an extra plot option called `peek` that makes a four-panel
plot with the lightcurve, phase-folded lightcurve, amplitude spectrum
*and* power spectrum (on log-log axes).

## FAQ

#### Why can't I get a lightcurve for my target?

First make sure it was actually observed by TESS! You can use the [Web TESS Viewing tool](https://heasarc.gsfc.nasa.gov/cgi-bin/tess/webtess/wtv.py) or query the star on MAST.  

Sometimes you need to widen `lightkurve`'s search from it's default (currently 0.0001 arcseconds) using the `-r` option, which specifies the search radius (also in arcseconds). For Algol (β Per), I found it necessary to increase the the search radius to 2 arcseconds.

#### Why does the lightcurve have unexpected features?

There are all sorts of systematic effects that might affect a lightcurve. `steve` defaults to the reduced PDCSAP flux so try using the SAP flux instead with the `--sap` option. 

If you want to see what I mean, try plotting the lightcurve of the Cepheid variable RT Aur and then try again with `--sap`...

## Examples

There might be slight differences because of version.

#### CO Cam, phase-folded lightcurve

    python3 steve.py fold co cam --sap -P -2 --title CO Cam --alpha 0.1 --phase-min 0.25

![CO Cam](https://pbs.twimg.com/media/EfYHAKrXsAA5TgV.png:small)

#### β Pic, amplitude spectrum

    python3 steve.py as bet pic --title β Pic --axis 0 95 0 1.1 --radius 2 --scale-y 1e3 --ylabel "amplitude (mmag)"

![β Pic](https://pbs.twimg.com/media/EfYFn0_XYAIxZth.png:small)

#### SX Phe, 50-period segments

    python3 steve.py fold sx phe --sap --title SX Phe -P -50 --delta 0.5 --ylabel TESSmag + offset

![SX Phe](https://pbs.twimg.com/media/EfYCLYgWoAIsKzJ.png:small)

## Footnotes

¹ At the moment, even though `lightkurve` also caches the data, it doesn't cache the search, so it still has to contact MAST. I could do something smarter but don't have the time.

² The amplitude and power spectra are computed using [Astropy](https://www.astropy.org/)'s [Lomb–Scargle method](https://docs.astropy.org/en/stable/api/astropy.timeseries.LombScargle.html).  The amplitude spectrum is by default oversampled by a factor of 10 (i.e. `steve` uses `samples_per_peak=10`).

³ For a nice, bright solar-like oscillator, try the red giant γ Men.
