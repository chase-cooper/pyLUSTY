This is a simple Python wrapper I wrote to facilitate the use of TLUSTY for 
producing models of DA(Z) white dwarf atmospheres and synthetic spectra, and 
for fitting photometry/spectroscopy with models. There are 4 main programs:

    - tlusty.py:    produce photospheric models and synthetic spectra
    - phot.py:      retrieve an object's photometry and fit it using pure-
                        hydrogen spectra
    - spec.py:      fit a user-supplied FITS file containing a spectrum using
                        pure-hydrogen spectra
    - linefit.py:   fit spectral lines for elemental abundances.

To use:

1)  Download TLUSTY/SYNSPEC. To access the most recent of TLUSTY and SYNSPEC, 
    TLUSTY208 and SYNSPEC54, you can follow the instructions given in Hubeny, 
    Allende Prieto, Osorio, and Lanz (2021), their Section 5. If you'd instead 
    prefer to use older versions of either code, TLUSTY205 and SYNSPEC51, these
    can be downloaded from https://as.arizona.edu/~hubeny/.

    Once the install is finished, I'd recommend placing all these files in a 
    folder called 'tlusty' in the top-level directory of the present code.

2)  Update 'consts.py': there are a handful of variables in 'consts.py' that 
    denote the location of important locations, such as where the TLUSTY and
    SYNSPEC executables are stored, that need to be edited to contain those
    paths.

3)  Try it out. The main function can be called as 

    $   python3 tlusty.py [name] [t_eff] [log_g] (synspec=1)

    where:
        -   'name' is the scenario name, which will denote the folder in which
            all outputs for this run are stored.
        -   't_eff' is the effective temperature of the white dwarf in Kelvin.
        -   'log_g' is the surface gravity of the white dwarf in cgs units

Necessary packages: astropy, astroquery, corner, dustmaps, emcee, matplotlib,
numpy, requests, scipy

    - for running TLUSTY/SYNSPEC: only base Python libraries
    - linefit.py: matplotlib, numpy, scipy
    - phot.py:  astroquery, dustmaps, corner, emcee, matplotlib, numpy, scipy,
        requests
    - spec.py:  matplotlib, numpy, scipy
    - utils.py: astropy, matplotlib, numpy, scipy

################################################################################

Directory:

pyLUSTY
    ├─ README.txt (you are here)
    ├─ data
    │   ├─ dustmaps
    │   │   ├─ sfd
    │   │       ├─ SFD_dust_4096_ngp.fits
    │   │       ├─ SFD_dust_4096_sgp.fits
    │   │
    │   ├─ filters
    │   │   ├─ 2mass
    │   │   │   ├─ j.dat, h.dat, k.dat
    │   │   ├─ galex
    │   │   │   ├─ fuv.dat, nuv.dat
    │   │   ├─ panstarrs
    │   │   │   ├─ g.dat, r.dat, i.dat, z.dat, y.dat
    │   │   ├─ sdss
    │   │   │   ├─ u.dat, g.dat, r.dat, i.dat, z.dat
    │   │   ├─ skymapper
    │   │       ├─ u.dat, g.dat, r.dat, i.dat, z.dat
    │   │
    │   ├─ fits_files (empty)
    │   │
    │   ├─ galex.dat
    │   ├─ hydrogengrid
    │   │   ├─ t10000_g6.00.spec, ..., t25000_g9.00.spec
    │   │
    │   ├─ mass_radius
    │       ├─ seq_020_thin.txt, ..., seq_130_thick.txt
    │
    ├─ hydrogengrid.zip
    ├─ linefit.py
    ├─ outputs (empty)
    ├─ phot.py
    ├─ spec.py
    ├─ tlusty
    │   ├─ scenarios (empty)
    │   ├─ cwd.flag
    │
    ├─ tlusty.py
    ├─ utils.py

