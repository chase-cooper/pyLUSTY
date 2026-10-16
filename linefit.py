from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import numpy as np
import subprocess,os,sys
from utils import *

##########################################################################################

### Produce a spectrum over a range [linecenter-width-2,linecenter+width+2] using the input files
###     of scenario 'scen_name' and varying the species [sp] to have abundances [abn]

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
        # Load in and normalize synthetic spectrum
        wvln_model,flux_model = forward_model(params,linecenter,width,spec_dict,conv_stuff,ism_center)
        flux_model = np.interp(wvln_obs,wvln_model,flux_model)

        # calculate residuals
        return flux_model - flux_obs

    # Run it
    result  = least_squares(residuals,p0,method='trf',xtol=1e-10,bounds=bounds)
    popt    = result.x

    ### Probability distribution function
    other_params    = popt[1:]
    abn_range       = np.arange(lower,upper,0.01)
    probs           = np.zeros_like(abn_range)
    for i,abn in enumerate(abn_range):
        params = [abn]
        for p in other_params: params.append(p)
        wvln_model,flux_model = forward_model(params,linecenter,width,spec_dict,conv_stuff)
        score = chi2(flux_obs,flux_model,sigma_obs)
        probs[i] = np.exp(-0.5*score)

    pdf_norm    = probs/sum(probs)
    cdf         = np.zeros_like(pdf_norm)
    for i in range(len(cdf)-1):
        cdf[i] = sum(pdf_norm[:i+1])/sum(pdf_norm)
    cdf[-1] = 1

    # Limits from CDF
    most_prob       = abn_range[np.argmax(pdf_norm)]
    lower_3sigma    = np.interp(0.01,cdf,abn_range)
    upper_3sigma    = np.interp(0.99,cdf,abn_range)

    ### Plot results
    fig,ax = plt.subplots(ncols=2)

    ### Best-fit model
    params = [most_prob]
    for p in other_params: params.append(p)
    wvln_model,flux_model = forward_model(params,linecenter,width,spec_dict,conv_stuff)
    score = chi2(flux_obs,flux_model,sigma_obs)

    # If applicable, add best-fitting ISM line
    if do_ism:
        ism_width,ism_amp = popt[2:]
        line_mask   = 1 - ism_amp*np.exp(-0.5*(wvln_model-ism_center)**2/(ism_width)**2)
        flux_model  *= line_mask
        ax.plot(wvln_model,line_mask,c='red',linestyle='dashed',zorder=100)

    ax[0].plot(wvln_model,flux_model,c='red',zorder=100,label=r"Best fit, $\chi^2=$"+str(np.around(score,2)))
    ax[0].plot(wvln_obs,flux_obs,c='black',label="Observed Spectrum")
    ax[0].set_xlabel(r"Wavelength [$\AA$]",size='large')
    ax[0].set_ylabel("Normalized flux",size='large')
    ax[0].set(xlim=(np.min(wvln_model),np.max(wvln_model)),
              ylim=(0,1.5))
    ax[0].legend()

    ax[1].plot(abn_range,pdf_norm,c='black')
    ax[1].set_xlabel("Abundance",size='large')
    ax[1].set_ylabel("Relative probability",size='large')
    ax[1].vlines(most_prob,ymin=0,ymax=max(pdf_norm)*1.1,colors='grey',linestyle='solid',zorder=-100,
                 label="Best fit: "+np.format_float_positional(most_prob,precision=2))
    ax[1].vlines(lower_3sigma,ymin=0,ymax=max(pdf_norm)*1.1,colors='grey',linestyle='dashed',zorder=-100,
                 label=r"$3\sigma$ lower limit: "+np.format_float_positional(lower_3sigma,precision=2))
    ax[1].vlines(upper_3sigma,ymin=0,ymax=max(pdf_norm)*1.1,colors='grey',linestyle='dashed',zorder=-100,
                 label=r"$3\sigma$ upper limit: "+np.format_float_positional(upper_3sigma,precision=2))
    ax[1].set_xlim(
        most_prob + (upper_3sigma - most_prob)*2.5,
        most_prob + (lower_3sigma - most_prob)*2.5
    )
    ax[1].invert_xaxis()
    ax[1].set_ylim(0,max(pdf_norm)*1.1)
    ax[1].legend()

    fig.set_figwidth(11)
    plt.tight_layout()
    plt.savefig(f'outputs/linefit_elem{sp}_{linecenter}AA.png',dpi=250)
    # plt.legend()
    plt.show()

if __name__=='__main__':
    # Because this function needs a fair but of info to run, the user is prompted for each
    # parameter.

    fname_spec  = input("Please provide the path to the FITS spectrum to be fit:\t")

    scen        = input("Please provide the name of the scenario to use as a base:\t")

    linecenter  = input("Please provide the central wavelength of the line/region to be fit:\t")
    linecenter  = float(linecenter)

    width       = input("Please provide the width of the region to be fit, in Angstroms (default: 5AA):\t")
    width       = float(width)

    sp_str      = input("Please provide the symbol or number of the element whose abundance you are\nfitting for:\t")
    try:
        sp = int(sp_str)
    except:
        sp = ElementtoZ(sp_str)
        if sp==-1:
            print("Element not supported: "+sp_str)
            exit()

    lower       = input("Please provide the lower limit on the abundance of that species (log(Z/H);\ndefault: -10):\t")
    if lower=='': lower = -10
    else:
        lower       = np.around(float(lower),2)

    upper       = input("Please provide the upper limit on the abundance of that species (log(Z/H);\ndefault: -4):\t")
    if upper=='': upper = -4
    else:
        upper       = np.around(float(upper),2)

    do_ism      = input("Is there an ISM line in the region you're fitting? [y/n]\t").lower()
    if do_ism=='y':
        do_ism      = True
        ism_center  = input("Provide the central wavelength of the ISM feature:\t")
        ism_center  = float(ism_center)
    else:
        do_ism      = False
        ism_center  = 0.

    print("Generating range of spectra...")
    newModelSpectra(scen,linecenter,width,sp,lower,upper,0.1)
    spectra = preload_spectra(scen,sp,lower,upper,0.1)

    print("Checking spectrum file...")
    instr   = fits.getheader(fname_spec)['INSTRUME']
    if instr=='COS':
        conv_stuff = loadConvolutionFiles(fname_spec)
    else:
        conv_stuff = 1

    print("Fitting...")
    lineFit(
        fname_spec=fname_spec,
        linecenter=linecenter,
        width=width,
        scen=scen,
        sp=sp,
        lower=lower,
        upper=upper,
        conv_stuff=conv_stuff,
        do_ism=do_ism,
        ism_center=ism_center,
        known_lines=[]
    )



