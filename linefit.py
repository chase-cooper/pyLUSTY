from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import numpy as np
import subprocess,os,sys
from utils import *

##########################################################################################

### Produce a spectrum over a range [linecenter-width-2,linecenter+width+2] using the input files
###     of scenario 'scen_name' and varying the species [sp] to have abundances [abn]

# Good
def newModelSpectra(scen_name,linecenter,width,sp,abn,upper=0,ddex=0.5):
    """
    Produce a spectrum over a wavelength range of [linecenter-width-2,linecenter+width+2],
    where the abundance of element 'sp' is set to 'abn'.
    """

    # Move to scenario folder
    os.chdir(f'tlusty/scenarios/{scen_name}/synspec/')

    # Save original spectrum if not already done
    if 'base.7' not in os.listdir():
        subprocess.run(['cp','fort.7','base.7'])

    ### Write input files
    # Write new fort.55 -- spectral range 
    unit55 = f"""       0      50       0
            1       0       0       1
            0       1       1       1       0
            1       1       0       0       0
            2       0       0
        {linecenter-width-5}    {linecenter+width+5}      10       0  0.00001     0.01
            0       0"""
    f = open('fort.55','w')
    f.write(unit55)
    f.close()

    # Iteratively produce new spectra
    for abn in np.arange(abn,upper+1e-4,ddex):
        if upper==0: abn=abn

        # Write new fort.56 -- elemental abundance change
        f = open('fort.56','w+')
        f.write(f'1\n')
        f.write(f"{sp} {10**abn}\n")  
        f.close()

        # run SYNSPEC
        try:
            os.system(f'{SYNEXE} < fort.5 > fort.6')
        except:
            print("\n Model run has encountered an error. Exiting...")
            sys.exit(0)

        # Save new spectrum
        abnstr = np.format_float_positional(abn,precision=2,min_digits=2)
        spec_name = f"{sp}={abnstr}"
        subprocess.run(['cp','fort.7',spec_name])

        if upper==0: break

    # Return to pwd directory
    os.chdir(HOME)

# Good
def preload_spectra(scen_name,sp,lower,upper,ddex=0.1):
    """
    Load a range of spectra, with each spectrum specified the abundance of a given element 'sp'.
    This function assumes that models with abundances between 'lo' and 'hi', and with differences
    between them of 'ddex', all exist.

    By default, ddex=0.1. The rest of the code here is designed with that assumption. USING A
    VALUE OF 'ddex' THAT IS NOT 0.1 WILL NOT WORK.
    """

    spec_dict = {
        'ddex':ddex     # specify the step size between abundances, to be retrieved later
    }

    upper += 1e-4       # Avoid rounding errors
    for abn in np.arange(lower,upper,ddex):
        abnstr              = np.format_float_positional(abn,precision=2,min_digits=2)
        modelname           = f"{SCEN}/{scen_name}/synspec/{sp}={abnstr}"
        model_data          = np.genfromtxt(modelname)
        wvln,flux           = model_data[:,0],model_data[:,1]
        spec_dict[abnstr]   = [np.float64(wvln),np.float64(flux)]
    return spec_dict

# Good (?)
def forward_model(params,linecenter,width,spec_dict,conv_stuff=1.,ism_center=0.):
    """
    Given a set of parameters and a set of spectra, produce a normalizedsynthetic spectrum. The 
    base spectrum is first built by interpolating between spectra with similar abundances, then 
    the spectrum is red/blueshifted according to a radial velocity parameter. The spectrum is
    convolved, either with an instrumental LSF (COS) or with a Gaussian, then normalized, and if
    applicable an ISM line is inserted.
    """

    ### Unpack parameters
    if len(params)==2:
        do_ism = False
        abn,vrad = params
    elif len(params)==4:
        do_ism = True
        abn,vrad,ism_width,ism_amp = params
    else:
        do_ism = False
        exit()
    
    ### Interpolate spectrum at the given abundance

    # Get upper and lower bounding abundances (accurate to one decimal places)
    ddex        = spec_dict['ddex']
    abn_lo      = max(abn - (ddex/2) + 1e-4,-10.)
    abn_hi      = min(abn + (ddex/2) + 1e-4,-4.)
    abnstr_lo   = np.format_float_positional(abn_lo,precision=1,min_digits=1)+'0'
    abnstr_hi   = np.format_float_positional(abn_hi,precision=1,min_digits=1)+'0'

    # Read in spectra; interpolate to same wavelengths
    wvln_lo,flux_lo = spec_dict[abnstr_lo]
    wvln_hi,flux_hi = spec_dict[abnstr_hi]
    flux_hi = np.interp(wvln_lo,wvln_hi,flux_hi)

    # Linearly interpolate
    wvln_base = np.copy(wvln_lo)
    flux_base = (abn - float(abnstr_lo))/ddex * flux_hi + (float(abnstr_hi) - abn)/ddex * flux_lo

    ### Apply wavelength shift
    beta        = vrad / 3e5
    wvln_shift  = wvln_base * np.sqrt( (1 + beta) / (1 - beta) )

    ### Convolve with COS LSF or a Gaussian
    if type(conv_stuff)==tuple:
        # You are convolving this spectrum with the COS instrument LSF.
        wvln_conv,flux_conv = convolveModelWithCOS(wvln_shift,flux_base,files=conv_stuff)
    else:
        wvln_conv,flux_conv = convolveModelWithGaussian(wvln_shift,flux_base,width=conv_stuff)
    

    ### Normalize
    # Normalize the spectrum
    blue_ave_wvln   = linecenter-width
    red_ave_wvln    = linecenter+width
    blue_ave_range  = np.where(np.abs(wvln_conv-blue_ave_wvln) <= 0.5)
    red_ave_range   = np.where(np.abs(wvln_conv-red_ave_wvln) <= 0.5)
    blue_ave        = np.average(flux_conv[blue_ave_range])
    red_ave         = np.average(flux_conv[red_ave_range])
    norm            = np.interp(wvln_conv,[blue_ave_wvln,red_ave_wvln],[blue_ave,red_ave])
    flux_norm       = flux_conv / norm

    # Restrict to wavelengths in the region of interest
    norm_range  = np.where(np.abs(wvln_conv-linecenter) <= width)
    wvln_norm   = wvln_conv[norm_range]
    flux_norm   = flux_norm[norm_range]

    ### Apply ISM line
    if do_ism:
        ism_factor  = 1 - ism_amp*np.exp(-0.5*(wvln_norm-ism_center)**2/ism_width**2)
        flux_final   = flux_norm * ism_factor
    else:
        flux_final = flux_norm

    return np.float64(wvln_norm),np.float64(flux_final)


def lineFit(fname_spec,scen,linecenter,width,sp,lower,upper,conv_stuff=0.,
            do_ism=False,ism_center=0.,known_lines=[],):
    """
    
    """

    ### Load spectral dictionary
    spec_dict = preload_spectra(scen,sp,lower,upper,ddex=0.1)

    ### Load in and normalize observed spectrum
    # Load in
    if conv_stuff!=0:
        wvln_obs,flux_obs,sigma = openCOSfile(fname_spec)
    elif 'xsh' in fname_spec:
        wvln_obs,flux_obs,sigma,gaussian = openXshooter(fname_spec)
        conv_stuff = 5000/gaussian
    elif 'uves' in fname_spec:
        wvln_obs,flux_obs,sigma,gaussian = openUVES(fname_spec)
        conv_stuff = 5000/gaussian
    else:
        exit()

    # Normalize
    flux_obs = flux_obs[np.abs(wvln_obs-linecenter) <= width]
    sigma_obs= sigma[np.abs(wvln_obs-linecenter) <= width]
    wvln_obs = wvln_obs[np.abs(wvln_obs-linecenter) <= width]

    # Normalize the observed spectrum and uncertainties
    blue_ave = np.average(flux_obs[np.abs(wvln_obs-linecenter+width) < 0.5])
    red_ave  = np.average(flux_obs[np.abs(wvln_obs-linecenter-width) < 0.5])
    norm_obs = np.interp(wvln_obs,[wvln_obs[0],wvln_obs[-1]],[blue_ave,red_ave])
    sigma_obs/= flux_obs
    flux_obs /= norm_obs

    # fig,ax = plt.subplots()
    # ax.plot(wvln_obs,flux_obs)

    # ### Mask known ISM lines
    # for line in known_lines:
    #     ism_center,ism_width,ism_amp = line
    #     line_mask   = 1 - ism_amp*np.exp(-0.5*(wvln_obs-ism_center)**2/(ism_width)**2)
    #     wvln_obs /= line_mask

    # ax.plot(wvln_obs,flux_obs)
    # plt.show()

    ### Run least squares

    # Define parameter bounds and initial guess
    if do_ism:
        bounds = [
            (lower,-100,0.01,0.01),
            (upper,100,0.1,1)
        ]
        p0 = [-7.,0,0.03,0.3]
    else:
        bounds = [
            (lower,-100),
            (upper,100)
        ]
        p0 = [-7.,0]

    # Define residuals function
    def residuals(params):
        print(params)
        # Load in and normalize synthetic spectrum
        wvln_model,flux_model = forward_model(params,linecenter,width,spec_dict,conv_stuff,ism_center)
        flux_model = np.interp(wvln_obs,wvln_model,flux_model)

        # calculate residuals
        return flux_model - flux_obs

    # Run it
    result  = least_squares(residuals,p0,method='trf',xtol=1e-8,bounds=bounds)
    popt    = result.x

    ### Plot results
    fig,ax = plt.subplots()
    wvln_model,flux_model = forward_model(popt,linecenter,width,spec_dict,conv_stuff)

    # If applicable, add best-fitting ISM line
    if do_ism:
        ism_width,ism_amp = popt[2:]
        line_mask   = 1 - ism_amp*np.exp(-0.5*(wvln_model-ism_center)**2/(ism_width)**2)
        flux_model  *= line_mask
        ax.plot(wvln_model,line_mask,c='red',linestyle='dashed',zorder=100)

    ax.plot(wvln_model,flux_model,c='red',zorder=100)
    ax.plot(wvln_obs,flux_obs,c='black')
    plt.show()

if __name__=='__main__':
    pass

# newModelSpectra('test_linelist_kurucz',4481,5,12,-10,-4,0.1)
# newModelSpectra('test_linelist_kurucz',1335,3,6,-10,-4,0.1)

res = preload_spectra('test_linelist_kurucz',6,-10,-4,0.1)
# p0  = [-5.61,35,4479,0.02,0.4]
# forward_model(p0,4481,4,res,0.1)
fname_spec = 'fits files/gd56_cos_fuv.fits'
conv_stuff = loadConvolutionFiles(fname_spec)
# fname_spec = 'fits files/sample_balmer/Gaia3296_xsh_nuv.fits'
# conv_stuff = 0
lineFit(
    fname_spec=fname_spec,
    linecenter=1335,
    width=5,
    scen='test_linelist_kurucz',
    sp=6,
    lower=-10,
    upper=-4,
    conv_stuff=conv_stuff,
    do_ism=True,
    ism_center=1334.54,
    known_lines=[]
)



