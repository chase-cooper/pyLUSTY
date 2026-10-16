pyLUSTY is a simple Python wrapper I wrote to facilitate the use of TLUSTY for 
producing models of DA(Z) white dwarf atmospheres and synthetic spectra, as 
well as for fitting photometry/spectroscopy with models. There are 4 main 
programs that constitute pyLUSTY:

    - tlusty.py:    produce photospheric models and synthetic spectra
    - phot.py:      retrieve an object's photometry and fit it using pure-
                        hydrogen spectra
    - spec.py:      fit a user-supplied FITS file containing a spectrum using
                        pure-hydrogen spectra
    - linefit.py:   fit spectral lines for elemental abundances.

Instructions for setting up the code environment, as well as the functionality
of the 4 programs, is explained below

################################################################################

#############
### SETUP ###
#############

Due to file size limitations, not everything that is necessary to run this code
is provided in the GitHub repository. Critically, the actual program files for
TLUSTY and SYNSPEC, as well as the many data files needed for their operation,
must be retrieved from the source and moved into this directory. To do so:

1)  Download TLUSTY/SYNSPEC. To access the most recent of TLUSTY and SYNSPEC, 
    TLUSTY208 and SYNSPEC54, you can follow the instructions given in Hubeny, 
    Allende Prieto, Osorio, and Lanz (2021), their Section 5 (available at
    https://arxiv.org/pdf/2104.02829). Only the program files for TLUSTY and
    SYNSPEC, the line list, and the iron data files are necessary for pyLUSTY; 
    the use of opacity tables is not yet implemented, and so downloading the
    related files is not necessary. If you'd instead prefer to use older 
    versions of either code, TLUSTY205 and SYNSPEC51, these can be downloaded 
    from https://as.arizona.edu/~hubeny/. 

    Once the install is finished, I'd recommend placing all these files in a 
    folder called 'tlusty' within the top-level directory of the present code.

2)  Update 'utils.py': near the top of 'utils.py', there is a section where a
    handful of constants are defined. These constants store the locations of
    important files and directories, and these need to be set properly before
    running any of the programs in pyLUSTY. The following constants need to be
    set:

        - HOME: set to the path of the top-level directory for this code, e.g.
            "/Users/chasecooper/pyLUSTY/"
        - PATH: set to the path of the directory containing all of the files
            from the installation of TLUSTY and SYNSPEC.
        - TLEXE: set to the path of the TLUSTY executable.
        - SYNEXE: set to the path of the SYNSPEC executable.
    
    The remaining constants are set relative to PATH, and do not necessarily
    need to be changed:

        - SCEN: the path of the directory where the outputs of each model run
            will be saved.
        - FLAG: the path of the keyword file used when running TLUSTY and
            SYNSPEC.

3)  Install the following Python packages:
        - astropy
        - astroquery
        - corner
        - dustmaps
        - emcee
        - matplotlib
        - numpy
        - requests
        - scipy

With that, the environment should now be setup and you can use pyLUSTY.

################################################################################

###########################
### PROGRAM DESCRIPTION ###
###########################

There are four main functions in pyLUSTY:

1) tlusty.py

tlusty.py is a Python wrapper that is meant to facilitate producing models of
DA(Z) white dwarfs with TLUSTY and SYNSPEC. The function is called as follows:

$   python3 tlusty.py [name] [t_eff] [log_g]

where:
    -   'name' is the scenario name, which is used to denote the folder in 
        which all outputs for this model run are stored. 
    -   't_eff' is the effective temperature of the white dwarf in Kelvin.
    -   'log_g' is the surface gravity of the white dwarf in cgs units

This is all that is strictly necessary to run a TLUSTY model. The results of
the run will be stored in a folder at the path PATH+name (e.g. 
/Users/chaseacooper/pyLUSTY/tlusty/scenarios/model_test). The final thermal
structure of the model is found at the end of the file 'lte.6'. For models 
without metals, the results are stored in a subfolder called 'lte'; for those 
with metals, they are stored in 'metals'.

In many cases, you'll also want the synthetic spectrum for the model. This can
be done by including the command line keyword 'synspec=1':

$   python3 tlusty.py [name] [t_eff] [log_g] synspec=1

The results are stored in the scenario subfolder 'synspec'. The final spectrum
in units of ergs/cm^2/s/AA is found in the file 'fort.7'.

2) phot.py

The photometric fitting code 'phot.py' uses the name of an object to query
different surveys for that object's photometry. Then, it runs a Monte Carlo-
Markov chain (MCMC) fitting algorithm that finds the best-fitting physical 
parameters for a given white dwarf. This is accomplished using a grid of pure-
hydrogen white dwarf models, which can be found in 'data/hydrogengrid'. The 
grid spans temperatures between 10,000K and 25,000K in steps of 500K, and 
log(g)'s between 6 and 9 in steps of 0.5 dex. The function is called as:

#   python3 phot.py [object_id]

where 'object_id' is the name of the object. Any name will work as long as
the name can be resolved by SIMBAD. The results of the MCMC fit will be stored
under the 'outputs' directory, with the relevant files having names that end
with the last 4 characters of the given name of the object.

3) spec.py

The spectroscopic fitting code 'spec.py' compares a given optical spectrum to
those from the hydrogen grid previously mentioned, once again using MCMC to get
the best-fitting temperature and log(g). The function is called as:

$   python3 spec.py name=[name] fname=[file name]

where
    - 'name' is the name of the object (only used for naming output files)
    - 'fname' is the file containing the object's spectrum*

(WIP)

4) linefit.py

(WIP)


################################################################################

#################
### DIRECTORY ###
#################

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

