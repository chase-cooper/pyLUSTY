from astropy.io import fits
from astropy.table import Table
from astropy.convolution import convolve
from scipy.interpolate import interp1d,LinearNDInterpolator
import matplotlib.pyplot as plt
import numpy as np
import os,urllib.request,time

###########################################################################################

# Top-level directory
HOME    = "/Users/cooper/Documents/Research/pyLUSTY/"    

# Paths to TLUSTY and SYNSPEC executables
TLEXE   = HOME+'tlusty/tlusty302/tlusty302.exe'
SYNEXE  = HOME+'tlusty/synspec51/synspec51.exe'

# Directory containing ALL TLUSTY-related files
PATH    = HOME+"tlusty/"    

# Scenario directory, where all output files will go
SCEN    = PATH+'scenarios'

# Path to the flag file containing keywords for TLUSTY/SYNSPEC
FLAG    = PATH+'cwd.flag'

###########################################################################################

### For TLUSTY/SYNSPEC operation:

def ZtoElement(Z):
    elements = {
        1:  'H',
        2:  'He',
        6:  'C',
        7:  'N',
        8:  'O',
        13: 'Al',
        14: 'Si',
        15: 'P',
        16: 'S'
    }
    return elements[Z]

def newSpeciesDict():
    """
    Return a species dictionary with every species' mode set to 0 except hydrogen.
    """
    return {
        'h':    {'mode': 2,'abn': 0.},       # Only hydrogen is on by default
        'he':   {'mode': 0,'abn': 0.},
        'li':   {'mode': 0,'abn': 0.},
        'be':   {'mode': 0,'abn': 0.},
        'b':    {'mode': 0,'abn': 0.},
        'c':    {'mode': 0,'abn': 0.},
        'n':    {'mode': 0,'abn': 0.},
        'o':    {'mode': 0,'abn': 0.},
        'fl':   {'mode': 0,'abn': 0.},
        'ne':   {'mode': 0,'abn': 0.},
        'na':   {'mode': 0,'abn': 0.},
        'mg':   {'mode': 0,'abn': 0.},
        'al':   {'mode': 0,'abn': 0.},
        'si':   {'mode': 0,'abn': 0.},
        'p':    {'mode': 0,'abn': 0.},
        's':    {'mode': 0,'abn': 0.},
        'cl':   {'mode': 0,'abn': 0.},
        'ar':   {'mode': 0,'abn': 0.},
        'k':    {'mode': 0,'abn': 0.},
        'ca':   {'mode': 0,'abn': 0.},
        'sc':   {'mode': 0,'abn': 0.},
        'ti':   {'mode': 0,'abn': 0.},
        'v':    {'mode': 0,'abn': 0.},
        'cr':   {'mode': 0,'abn': 0.},
        'mn':   {'mode': 0,'abn': 0.},
        'fe':   {'mode': 0,'abn': 0.},
        'co':   {'mode': 0,'abn': 0.},
        'ni':   {'mode': 0,'abn': 0.}
    }

def writeUnit5(path:str,teff:float,log_g:float,species_dict:dict=newSpeciesDict(),is_lte:bool=True,is_ltgray:bool=True):
    """
    Helper function that writes the unit 5 files for TLUSTY model runs.

    Inputs:
    - path (str):           the location of the new unit 5 file
    - teff (float):         the temperature value to be written
    - log_g (float):        the value of log(g) to be written
    - species_dict (dict):  the dictionary containing all elements, their modes of treatment, and
                            their abundances.
    - is_lte (bool):        whether or not the model to be run is an LTE model
    - is_ltgray (bool):     whether or not the model to be run is an LTGray model. If false, a 
                            previous model output must be provided as fort.8.
    """

    # First block: temp, gravity, LTE/LTGRAY flags, flag file, frequency points
    l = 'T' if is_lte else 'F'
    g = 'T' if is_ltgray else 'F'

    f = open(path+'/fort.5','w+')
    f.write(str(int(teff))+'\t'+np.format_float_positional(log_g,precision=3)+'\t! T_eff, log(g)\n')
    f.write(f'{l}\t{g}       ! LTE, LTGRAY\n')
    f.write('\'cwd.flag\' ! Add\'l params\n')
    f.write('200        ! frequency points(?)\n')

    # Second block: species, modes of computation, and abundances
    num_exp = len(species_dict.keys())
    f.write(str(num_exp)+'\t\t! atoms\n')
    f.write('* mode abn modpf\n')
    for key in species_dict.keys():
        # if key=='h':continue
        f.write(f"  {species_dict[key]['mode']}   {species_dict[key]['abn']}   0   ! {key}\n")
    
    # Third block: Explicit ions
    f.write('*\n* explicit ions\n*\n')
    f.write('* iat iz nlevs ilast ilvlin nonstd typion file\n*\n')

    # Include files for relevant species. Do not set mode=2 to any species but these!!!
    if species_dict['h']['mode']==2:
        f.write("""  1   0   16   0   0     0   ' H 1' './data/h1s16.dat'
  1   1   1    1   0     0   ' H 2' ' '\n""")
#     if species_dict['he']['mode']==2:
#         f.write("""  2   0   24   0   0     0   'He 1' './data/he1.dat'
#   2   1   20   0   0     0   'He 2' './data/he2.dat'
#   2   2   1    1   0     0   'He 3' ' '\n""")
    if species_dict['c']['mode']==2:
        f.write("""  6   0   40   0   0     0   ' C 1' './data/c1.dat'
  6   1   22   0   0     0   ' C 2' './data/c2.dat'
  6   2   46   0   0     0   ' C 3' './data/c3_34+12lev.dat'
  6   3   1    1   0     0   ' C 4' ' '\n""")
#   6   3   25   0   0     0   ' C 4' './data/c4.dat'
#   6   4   1    1   0     0   ' C 5' ' '\n""")
    if species_dict['n']['mode']==2:
        f.write("""  7   0   34   0   0     0   ' N 1' './data/n1.dat'
  7   1   1    1   0     0   ' N 2' ' '\n""")
#   7   1   42   0   0     0   ' N 2' './data/n2_32+10lev.dat'
#   7   2   32   0   0     0   ' N 3' './data/n3.dat'
#   7   3   48   0   0     0   ' N 4' './data/n4_34+14lev.dat'
#   7   4   16   0   0     0   ' N 5' './data/n5.dat'
#   7   5   1    1   0     0   ' N 6' ' '\n""")
    if species_dict['o']['mode']==2:
        f.write("""  8   0   33   0   0     0   ' O 1' './data/o1_23+10lev.dat'
  8   1   1    1   0     0   ' O 1' ' '\n""")
#   8   1   48   0   0     0   ' O 2' './data/o2_36+12lev.dat'
#   8   2   41   0   0     0   ' O 3' './data/o3_28+13lev.dat'
#   8   3   39   0   0     0   ' O 4' './data/o4.dat'
#   8   4   6    0   0     0   ' O 5' './data/o5.dat'
#   8   5   1    1   0     0   ' O 6' ' '\n""")
#     if species_dict['na']['mode']==2:
#         f.write("""  11  0   32   0   0     0   'Na 1' './data/na1.t'
#   11  1   8    1   0     0   'Na 2' ' '\n""")
    if species_dict['mg']['mode']==2:
        f.write("""  12  1   25   0   0     0   'Mg 2' './data/mg2.dat'
  12  2   1    1   0     0   'Mg 3' ' '\n""")
    if species_dict['al']['mode']==2:
        f.write("""  13  1   29   0   0     0   'Al 2' './data/al2_20+9lev.dat'
  13  2   23   0   0     0   'Al 3' './data/al3_19+4lev.dat'
  13  3   1    1   0     0   'Al 4' ' '\n""")
    if species_dict['si']['mode']==2:
        f.write("""  14  1   40   0   0     0   'Si 2' './data/si2_36+4lev.dat'
  14  2   30   0   0     0   'Si 3' './data/si3.dat'
  14  3   23   0   0     0   'Si 4' './data/si4.dat'
  14  4   1    1   0     0   'Si 5' ' '\n""")
    if species_dict['s']['mode']==2:
        f.write("""  16  1   33   0   0     0   ' S 2' 'data/s2_23+10lev.dat'
  16  2   41   0   0     0   ' S 3' './data/s3_29+12lev.dat'
  16  3   1    1   0     0   ' S 4' ' '\n""")
#   16  3   38   0   0     0   ' S 4' './data/s4_33+5lev.dat'
#   16  4   25   0   0     0   ' S 5' './data/s5_20+5lev.dat'
#   16  5   1    1   0     0   ' S 6' ' '\n""")
    if species_dict['ca']['mode']==2:
        f.write("""  20  1   32   0   0     0   'Ca 2' './data/ca2.t'
  20  2   1    1   0     0   'Ca 3' ' '\n""")
    if species_dict['fe']['mode']==2:
        f.write("""  26  1   36   0   0     -1  'Fe 2' './data/fe2va.dat'
   0  0                             './data/gf2601.gam'
                                    './data/gf2601.lin'
                                    './data/fe2p_14+11lev.rap'
  26  2   50   0   0     -1  'Fe 3' './data/fe3va.dat'
   0  0                             './data/gf2602.gam'
                                    './data/gf2602.lin'
                                    './data/fe3p_22+7lev.rap'
  26  3   1    1   0     0   'Fe 4' ' '\n""")
#   26  3   43   0   0     -1  'Fe 4' './data/fe4va.dat'
#    0  0                             './data/gf2603.gam'
#                                     './data/gf2603.lin'
#                                     './data/fe4p_21+11lev.rap'
#   26  4   42   0   0     -1  'Fe 5' './data/fe5va.dat'
#    0  0                             './data/gf2604.gam'
#                                     './data/gf2604.lin'
#                                     './data/fe5p_19+11lev.rap'
#   26  5    1   1   0     0   'Fe 6' ' '\n""")

    f.write("  0   0   0   -1   0     0   '    ' ' '")
    f.close()

def shrinkScenarios():
    """
    Shrink scenario folder sizes by removing large files that are not used.
    """
    for scen in os.listdir(SCEN):
        if scen[0]=='.': continue
        if 'fort.81' in os.listdir(f'{SCEN}/{scen}/lte'):
            os.remove(f'{SCEN}/{scen}/lte/fort.81')
        if 'fort.82' in os.listdir(f'{SCEN}/{scen}/lte'):
            os.remove(f'{SCEN}/{scen}/lte/fort.82')
        if 'fort.83' in os.listdir(f'{SCEN}/{scen}/lte'):
            os.remove(f'{SCEN}/{scen}/lte/fort.83')

##########################################################################################

# ### Functions to convolve model spectra with the proper COS LSF
def lsfParams(fname:str):
    """
    Get relevant parameters to determine which LSF file you need for a given
    FITS file. Borrowed from HST Notebook on LSF convolution
    """
    # Select the primary header
    fuvHeader0 = fits.getheader(fname, ext=0)
    # print(f"For the file {fname}, the relevant parameters are: ")

    # Make a dictionary to store what you find here
    param_dict = {}

    # We want data for the FUV detector, G130M grating at LP3 cenwave 1291
    keywords = ["DETECTOR", "OPT_ELEM", "LIFE_ADJ", "CENWAVE", "DISPTAB"]

    # Print out the relevant values:
    for hdrKeyword in keywords:
        # For DISPTAB
        try:
            # Save the key/value pairs to the dictionary
            value = fuvHeader0[hdrKeyword].split("$")[1]
            # DISPTAB needs the split here
            param_dict[hdrKeyword] = value
        # For other params
        except (IndexError, AttributeError):
            # Save the key/value pairs to the dictionary
            value = fuvHeader0[hdrKeyword]
            param_dict[hdrKeyword] = value

        # Print the key/value pairs
        # print(f"{hdrKeyword} = {value}")

    return param_dict

def fetch_files(det, grating, lpPos, cenwave, disptab):
    """
    Given all the inputs, this will download both
    the LSF and Disptab files to use in the convolution and return their paths.

    Input: 
    det (str): The detector used
    grating (str): Type of grating used
    lpPos (str): Lifetime position used
    cenwave (str): Central wavelength used
    disptab (str): DISPTAB used (will get the path in the function)

    Returns:
    LSF_file_name (str): filename of the new downloaded LSF file
    disptab_path (str): path to the new downloaded disptab file

    Borrowed from HST Notebook on LSF convolution
    """
    datadir = 'hstdata'

    # Link to where all the files live
    COS_site_rootname = (
        "https://www.stsci.edu/files/live/sites/www/files/"
        "home/hst/instrumentation/cos/"
        "performance/spectral-resolution/_documents/"
    )  
    # print(det)

    # Only one file for NUV
    if det == "NUV":
        LSF_file_name = "nuv_model_lsf.dat"

    # FUV files follow a naming pattern
    elif det == "FUV":
        LSF_file_name = f"aa_LSFTable_{grating}_{cenwave}_LP{lpPos}_cn.dat"

    # Where to find file online
    LSF_file_webpath = COS_site_rootname + LSF_file_name
    urllib.request.urlretrieve(
        LSF_file_webpath, f"{datadir}/{LSF_file_name}"
    )
    
    # Where to save file to locally
    # print(f"Downloaded LSF file to {f"{datadir}/{LSF_file_name}"}")

    # And we'll need to get the DISPTAB file as well
    disptab_path = f"{datadir}/{disptab}"
    urllib.request.urlretrieve(
        f"https://hst-crds.stsci.edu/unchecked_get/references/hst/{disptab}",
        disptab_path
    )
    
    # print(f"Downloaded DISPTAB file to {disptab_path}")

    return LSF_file_name, disptab_path

def read_lsf(filename):
    # This is the table of all the LSFs: called "lsf"
    # The first column is a list of the wavelengths corresponding to the line profile, so we set our header accordingly
    # If its an NUV file, header starts 1 line later
    # Borrowed from HST Notebook on LSF convolution
    if "nuv_" in filename:
        ftype = "nuv"
        # print(f"Detector used: {ftype}")
        hs = 1

    # Otherwise, assume its an FUV file
    else:
        ftype = "fuv"
        # print(f"Detector used: {ftype}")

    hs = 0
    lsf = Table.read(filename,
                     format="ascii",
                     header_start=hs)

    # This is the range of each LSF in pixels (for FUV from -160 to +160, inclusive)
    # The middle pixel of the lsf is considered zero; the center is relative zero

    # Integer division to yield whole pixels
    pix = np.arange(len(lsf)) - len(lsf) // 2

    # The column names returned as integers.
    lsf_wvlns = np.array([int(float(k)) for k in lsf.keys()])

    return lsf, pix, lsf_wvlns

def get_disp_params(disptab, cenwave, segment, x=[]):
    """
    Helper function to redefine_lsf(). Reads through a DISPTAB file and gives relevant
    dispersion relationship/wavelength solution over input pixels.

    Parameters:
    disptab (str): Path to your DISPTAB file.
    cenwave (str): Cenwave for calculation of dispersion relationship.
    segment (str): FUVA or FUVB?
    x (list): Range in pixels over which to calculate wvln with dispersion relationship (optional).
    
    Returns:
    disp_coeff (list): Coefficients of the relevant polynomial dispersion relationship
    wavelength (list; if applicable): Wavelengths corresponding to input x pixels 

    Borrowed from HST Notebook on LSF convolution
    """
    with fits.open(disptab) as d:
        wh_disp = np.where(
            (d[1].data["cenwave"] == cenwave)
            & (d[1].data["segment"] == segment)
            & (d[1].data["aperture"] == "PSA")
        )[0]
        # 0 is needed as this returns nested list [[arr]]
        disp_coeff = d[1].data[wh_disp]["COEFF"][0]
    
    # If given a pixel range, build up a polynomial wvln solution pix -> λ
    if len(x):
        wavelength = np.polyval(p=disp_coeff[::-1], x=np.arange(16384))
        return disp_coeff, wavelength
    
    # If x is empty:
    else:
        return disp_coeff

def redefine_lsf(lsf_file, cenwave, disptab, detector="FUV"):
    """
    Helper function to convolve_lsf(). Converts the LSF kernels in the LSF file from a fn(pixel) -> fn(λ)
    which can then be used by convolve_lsf() and re-bins the kernels.

    Parameters:
    lsf_file (str): path to your LSF file
    cenwave (str): Cenwave for calculation of dispersion relationship
    disptab (str): path to your DISPTAB file
    detector (str): FUV or NUV?

    Returns:
    new_lsf (numpy.ndarray): Remapped LSF kernels.
    new_w (numpy.ndarray): New LSF kernel's LSF wavelengths.
    step (float): first order coefficient of the FUVA dispersion relationship; proxy for Δλ/Δpixel.

    Borrowed from HST Notebook on LSF convolution
    """

    if detector == "FUV":
        xfull = np.arange(16384)

        # Read in the dispersion relationship here for the segments
        # FUVA is simple
        disp_coeff_a, wavelength_a = get_disp_params(disptab,
                                                     cenwave,
                                                     "FUVA",
                                                     x=xfull)
        
        # FUVB isn't taken for cenwave 1105, nor 800:
        if (cenwave != 1105) & (cenwave != 800):
            disp_coeff_b, wavelength_b = get_disp_params(disptab,
                                                         cenwave,
                                                         "FUVB",
                                                         x=xfull)
        elif cenwave == 1105:
            # 1105 doesn't have an FUVB so set it to something arbitrary and clearly not real:
            wavelength_b = [-99.0, 0.0]

        # Get the step size info from the FUVA 1st order dispersion coefficient
        step = disp_coeff_a[1]

        # Read in the lsf file
        lsf, pix, w = read_lsf(lsf_file)

        # Take median spacing between original LSF kernels
        deltaw = np.median(np.diff(w))

        # Resamples if the spacing of the original LSF wvlns is too narrow
        if (deltaw < len(pix) * step * 2):  
            # This is all a set up of the bins we want to use
            # The wvln difference between kernels of the new LSF should be about twice their width
            new_deltaw = round(len(pix) * step * 2.0)  

            # nw = Number of LSF wavelengths
            new_nw = (int(round((max(w) - min(w)) / new_deltaw)) + 1)  

            # New version of lsf_wvlns
            new_w = min(w) + np.arange(new_nw) * new_deltaw  

            # Populating the lsf with the proper bins:
            # Empty 2-D array to populate
            new_lsf = np.zeros((len(pix), new_nw)) 

            for i, current_w in enumerate(new_w):
                # Find closest original LSF wavelength to new LSF wavelength
                dist = abs(current_w - w)  
                lsf_index = np.argmin(dist)
                # Column name corresponding to closest orig LSF wvln
                orig_lsf_wvln_key = lsf.keys()[lsf_index]  
                # Assign new LSF wvln the kernel of the closest original lsf wvln
                new_lsf[:, i] = np.array(lsf[orig_lsf_wvln_key])  

        else:
            new_lsf = lsf
            new_w = w

        return new_lsf, new_w, step

    elif detector == "NUV":
        xfull = np.arange(1024)

        # Read in the dispersion relationship here for the segments
        disp_coeff_a, wavelength_a = get_disp_params(disptab,
                                                     cenwave,
                                                     "NUVA",
                                                     x=xfull)
        
        disp_coeff_b, wavelength_b = get_disp_params(disptab,
                                                     cenwave,
                                                     "NUVB",
                                                     x=xfull)
        
        disp_coeff_c, wavelength_c = get_disp_params(disptab,
                                                     cenwave,
                                                     "NUVC",
                                                     x=xfull)

        # Get the step size info from the NUVB 1st order dispersion coefficient
        step = disp_coeff_b[1]

        # Read in the lsf file
        lsf, pix, w = read_lsf(lsf_file)

        # Take median spacing between original LSF kernels
        deltaw = np.median(np.diff(w))

        # This section is a set up of the new bins we want to use:
        # The wvln difference between kernels of the new LSF should be about twice their width
        new_deltaw = round(len(pix) * step * 2.0) 

        # nw = Number of LSF wavelengths
        new_nw = (int(round((max(w) - min(w)) / new_deltaw)) + 1)  

        # New version of lsf_wvlns
        new_w = min(w) + np.arange(new_nw) * new_deltaw  

        # Populating the lsf with the proper bins:
        # Empty 2-D array to populate
        new_lsf = np.zeros((len(pix), new_nw))  

        for i, current_w in enumerate(new_w):
            # Find closest original LSF wavelength to new LSF wavelength
            dist = abs(current_w - w)  
            lsf_index = np.argmin(dist)

            # Column name corresponding to closest orig LSF wvln
            orig_lsf_wvln_key = lsf.keys()[lsf_index]  
            
            # Assign new LSF wvln the kernel of the closest original lsf wvln
            new_lsf[:, i] = np.array(lsf[orig_lsf_wvln_key])  
            
        return new_lsf, new_w, step

def convolve_lsf(wavelength, spec, cenwave, lsf_file, disptab, detector="FUV"):
    """
    Main function; Convolves an input spectrum - i.e. template or STIS spectrum - with the COS LSF.

    Parameters:
    wavelength (list or array): Wavelengths of the spectrum to convolve.
    spec (list or array): Fluxes or intensities of the spectrum to convolve.
    cenwave (str): Cenwave for calculation of dispersion relationship
    lsf_file (str): Path to your LSF file
    disptab (str): Path to your DISPTAB file
    detector (str) : Assumes an FUV detector, but you may specify 'NUV'.

    Returns:
    wave_cos (numpy.ndarray): Wavelengths of convolved spectrum.!Different length from input wvln
    final_spec (numpy.ndarray): New LSF kernel's LSF wavelengths.!Different length from input spec

    Borrowed from HST Notebook on LSF convolution.
    """
    # First calls redefine to get right format of LSF kernels
    new_lsf, new_w, step = redefine_lsf(lsf_file,
                                        cenwave,
                                        disptab,
                                        detector=detector)

    # Sets up new wavelength scale used in the convolution
    nstep = round((max(wavelength) - min(wavelength)) / step) - 1
    wave_cos = min(wavelength) + np.arange(nstep) * step

    # Resampling onto the input spectrum's wavelength scale:
    # Builds up interpolated function from input spectrum
    interp_func = interp1d(wavelength, spec)  

    # Builds interpolated initial spectrum at COS' wavelength scale for convolution
    spec_cos = interp_func(wave_cos)  

    # Initializes final spectrum to the interpolated input spectrum
    final_spec = interp_func(wave_cos)  

    # Loop through the redefined LSF kernels:
    for i, w in enumerate(new_w):  
        # First need to find the boundaries of each kernel's "jurisdiction": where it applies
        # The first and last elements need to be treated separately

        # First kernel:
        if i == 0:  
            diff_wave_left = 500
            diff_wave_right = (new_w[i + 1] - w) / 2.0

        # Last kernel
        elif i == len(new_w) - 1:  
            diff_wave_right = 500
            diff_wave_left = (w - new_w[i - 1]) / 2.0

        # All other kernels
        else:  
            diff_wave_left = (w - new_w[i - 1]) / 2.0
            diff_wave_right = (new_w[i + 1] - w) / 2.0

        # Splitting up the spectrum into slices around the redefined LSF kernel wvlns
        # Will apply the kernel corresponding to that chunk to that region of the spectrum - its "jurisdiction"
        chunk = np.where(
            (wave_cos < w + diff_wave_right) & (wave_cos >= w - diff_wave_left)
        )[0]
        if len(chunk) == 0:
            # Off the edge, go to the next chunk
            continue

        # Selects the current kernel
        current_lsf = new_lsf[:, i]  

        if len(chunk) >= len(
            current_lsf
        ):  # Makes sure that the kernel is smaller than the chunk
            final_spec[chunk] = convolve(
                spec_cos[chunk],
                # Applies the actual convolution
                current_lsf,  
                boundary="extend",
                normalize_kernel=True,
            )

    # Remember, not the same length as input spectrum data!
    return wave_cos, final_spec  

#####

def loadConvolutionFiles(fname):
    """
    A function to extract instrument parameters and necessary LSF files to enable
    convolution. Separate from the main function so it can just be called once.
    """
    params  = lsfParams(fname)
    lsf_file,disp_file  = fetch_files(*params.values())
    lsf_file = 'hstdata/'+lsf_file
    return params,lsf_file,disp_file

def convolveModelWithCOS(wvln,flux,fname,files=None):
    """
    A function that simply puts all the previous steps together. Returns the
    convolved model spectrum
    """

    if files is None:
        files = loadConvolutionFiles(fname)
    params,lsf_file,disp_file = files

    wvln_m,spec_m = wvln,flux
    wvln_m,spec_m = convolve_lsf(wvln_m,spec_m,cenwave=params['CENWAVE'],lsf_file=lsf_file,disptab=disp_file,detector=params['DETECTOR'])
    return wvln_m,spec_m

##########################################################################################

### Functions to open observed and model spectra

def openModel(fnamem:str):
    """
    Read model spectrum
    """
    dat = np.loadtxt(fnamem)
    return dat[:,0],dat[:,1]

def openXshooter(fname:str):
    """
    Open a fits file containing data from the XSHOOOTER echelle spectrograph
    """

    file = fits.open(fname)
    if file[0].header['FLUXCAL'] != "ABSOLUTE":
        print("Spectrum is not flux calibrated. Exiting...")
        exit()
    res  = file[0].header['SPEC_RES']
    data = file[1].data
    file.close()

    # The data
    wvln0 = data["WAVE"][0] * 10    # nm to AA
    flux0 = data["FLUX"][0]
    sigma0 = data["ERR"][0]
    return wvln0,flux0,sigma0,res

def openUVES(fname:str):
    """
    
    """

    file    = fits.open(fname)
    # print(repr(file[0].header))
    is_cal  = file[0].header['FLUXCAL']
    if is_cal != 'ABSOLUTE':
        print("WARNING: spectrum at \'"+fname+"\' is not flux calibrated.")
        exit(0)
    res     = file[0].header['SPEC_RES']
    data    = file[1].data
    file.close()

    wvln = data["WAVE"][0]
    flux = data["FLUX"][0]
    sigma = data["ERR"][0]

    return wvln,flux,sigma,res

def openCOSfile(fname:str):
    """
    Extract the spectrum from a single .fits file.

    I'm not sure that all files I'll be working with share the same format, so
    this function may be edited on the fly.
    """
    with fits.open(fname) as hdu:

        ### Extract data
        header      = fits.getheader(fname)
        detector    = header['DETECTOR']
        data        = Table.read(fname)

        # NUV
        if detector=='NUV':

            wvln_nuva   = data['WAVELENGTH'][0]
            wvln_nuvb   = data['WAVELENGTH'][1]
            wvln        = np.concat([wvln_nuva,wvln_nuvb])

            flux_nuva   = data['FLUX'][0]
            flux_nuvb   = data['FLUX'][1]
            flux        = np.concat([flux_nuva,flux_nuvb])

            err_nuva    = data['ERROR'][0]
            err_nuvb    = data['ERROR'][1]
            err         = np.concat([err_nuva,err_nuvb])

        elif detector=='FUV':
            wvln_fuva   = data['WAVELENGTH'][0]
            wvln_fuvb   = data['WAVELENGTH'][1]
            wvln        = np.concat([wvln_fuvb,wvln_fuva])

            flux_fuva   = data['FLUX'][0]
            flux_fuvb   = data['FLUX'][1]
            flux        = np.concat([flux_fuvb,flux_fuva])

            err_nuva    = data['ERROR'][0]
            err_nuvb    = data['ERROR'][1]
            err         = np.concat([err_nuvb,err_nuva])

        else:
            print("Invalid COS file: detector is not 'NUV' or 'FUV'. Exiting")
            exit()

        # Remove zeros
        nonzero     = np.where(flux != 0)
        wvln        = wvln[nonzero]
        flux        = flux[nonzero]
        err         = err[nonzero]

        hdu.close()

    return wvln,flux,err

################################################################################################

### Other

def chi2(y,ymodel,yerr):
    """
    Calculate the reduced chi-squared score of a model fit. Points with no corresponding 
    uncertainty are not considered.
    """

    # Remove points where there's no uncertainty value
    y       = y[yerr > 0]
    ymodel  = ymodel[yerr > 0]
    yerr    = yerr[yerr > 0]

    # Return chi^2 value
    res = np.sum(np.nan_to_num((y-ymodel)**2/yerr**2,nan=0))
    return res

def ABmag2flux(mag,wvln=None):
    """
    Calculate an object's flux given its AB magnitude. Flux can either be in units of
    Janskys (default) or W/m^2/nm. In the case of the latter, an effective wavelength
    is needed.
    """

    c = 3e8
    fluxJanskys = 10**(3.56 - 0.4*mag)                  # Jy = 10^-26 W/m^2/Hz
    # flux        = 1e-26*fluxJanskys*c/(wvln*wvln)   # W/m^2/nm
    return fluxJanskys
ABmag2flux = np.vectorize(ABmag2flux)

def convolveModelWithGaussian(wvln,flux,width=3.):
    """
    Convolve a provided model spectrum (wvln [AA] and flux) with a Gaussian. The width of the
    Gaussian is determined automatically from the spacing of the observed spectrum.
    """

    mid_wvln    = wvln[len(wvln)//2]

    # Build Gaussian over a small wavelength range
    range_gauss = np.where(np.abs(wvln-mid_wvln) <= 5*width)
    wvln_gauss  = wvln[range_gauss]
    gaussian    = np.exp(-0.5*( (wvln_gauss-mid_wvln) / width )**2)

    # Convolve
    conv_flux   = np.convolve(flux,gaussian,mode='same') / np.sum(gaussian)

    # fig,ax = plt.subplots()
    # ax.plot(wvln,flux)
    # ax.plot(wvln,conv_flux)
    # plt.show()
    # plt.close()

    return wvln,conv_flux

def loadSpectra(shrink_factor:int=1):
    # Spectrum dictionary
    start = time.time()
    print("Loading spectra...",end='')

    spec_dict = {}
    for t in np.arange(10000,25001,500):
        tstr = str(t)
        for g in np.arange(6,9.1,0.5):
            gstr = np.format_float_positional(g,precision=2,min_digits=2)

            scen = f"t{tstr}_g{gstr}"
            spec_dict[scen] = {}

            dat = np.genfromtxt(f'data/hydrogengrid/{scen}.spec')
            fwvln = dat[:,0][::shrink_factor]
            fflux = dat[:,1][::shrink_factor]
            spec_dict[scen] = {
                'wvln':fwvln,
                'flux':fflux
            }

    end = time.time()
    print(f"Loaded spectra in {round(end-start,2)}s")
    return spec_dict

def interp_spectrum_dict(t,logg,spec_dict,plotit:bool=False):
    """
    Given a temperature (anything) and log(g), linearly interpolate between
    adjacent spectra to produce a new spectrum. Instead of reading in data from files, the data
    are accessed through a dictionary where all the data has been preloaded.
    """
    t_lo        = t - (t % 500)
    t_lo_str    = str(int(t_lo))
    t_hi        = t_lo + 500
    t_hi_str    = str(int(t_hi))
    diff_t      = t - t_lo
    # print(t_lo_str,t_hi_str,diff_t)

    logg_lo     = logg - (logg % 0.5)
    logg_lo_str = np.format_float_positional(logg_lo,min_digits=2)
    logg_hi     = logg_lo + 0.5
    logg_hi_str = np.format_float_positional(logg_hi,min_digits=2)
    diff_logg   = logg - logg_lo
    # print(logg_lo_str,logg_hi_str,diff_logg)

    # Bilinear interpolation over model grid
    spec00  = spec_dict[f"t{t_lo_str}_g{logg_lo_str}"]
    spec01  = spec_dict[f"t{t_lo_str}_g{logg_hi_str}"]
    spec10  = spec_dict[f"t{t_hi_str}_g{logg_lo_str}"]
    spec11  = spec_dict[f"t{t_hi_str}_g{logg_hi_str}"]

    # Interpolate fluxes to lower model wavelength spectrum
    wvln    = spec00['wvln']

    spec0x  = spec00['flux']*(1-diff_logg/0.5) + np.interp(wvln,spec01['wvln'],spec01['flux']) * diff_logg/0.5
    spec1x  = np.interp(wvln,spec10['wvln'],spec10['flux'])*(1-diff_logg/0.5) + np.interp(wvln,spec11['wvln'],spec11['flux'])*diff_logg/0.5

    spec    = spec0x*(1-diff_t/500) + spec1x*(diff_t/500)

    if plotit:
        fig,ax = plt.subplots(ncols=3)

        ax[0].plot(spec00['wvln'],spec00['flux'],label='spec00') 
        ax[0].plot(spec01['wvln'],spec01['flux'],label='spec01')
        ax[0].plot(wvln,spec0x,label='spec0x')                     

        ax[1].plot(spec10['wvln'],spec10['flux'],label='spec10')
        ax[1].plot(spec11['wvln'],spec11['flux'],label='spec11')
        ax[1].plot(wvln,spec1x,label='spec1x')

        ax[2].plot(wvln,spec0x,label='spec0x')  
        ax[2].plot(wvln,spec1x,label='spec1x')  
        ax[2].plot(wvln,spec,label='final spec')  

        for a in ax:
            a.legend()
        plt.show()
        plt.close()

    return wvln,spec

################################################################################################

### interpolate over the Bédard+2020 models to get mass from a white dwarf's radius and T_eff.
###     For use in photometric fitting
start = time.time()

atmo_size   = 'thin'    # Mass of hydrogen layer: 'thin' (q=1e-10) or 'thick' (q=1e-4)
points      = None      # arrays of [t_eff,radius]
values_m    = []        # values of mass
values_logg = []        # values of log(g)

for mass in np.arange(20,131,5):
    mass_str = '0'*bool(mass<100) + str(int(mass))

    with open(f'data/mass_radius/seq_{mass_str}_{atmo_size}.txt') as file:
        lines = file.readlines()[5:-1]
        for modelnum in range(0,len(lines),3):
            entry = lines[modelnum].split()
            teff    = float(entry[1])
            logg    = float(entry[2])
            rad     = float(entry[3])

            if type(points) != np.ndarray:
                points = np.array([teff,rad])
            else:
                points  = np.vstack([points,np.array([teff,rad])])

            values_m.append(mass/100)
            values_logg.append(logg)
            # print(points[-1])
        file.close()

interpBedardMass = LinearNDInterpolator(points,values_m)        # derived value of log(g) is the same to <1%
interpBedardLogg = LinearNDInterpolator(points,values_logg)

print(f"Loaded Bédard model interpolators in {round(time.time()-start,4)}s")

################################################################################################


