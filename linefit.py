from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import numpy as np
import subprocess,os,sys
from utils import *

##########################################################################################

### Produce a spectrum over a range [linecenter-width-2,linecenter+width+2] using the input files
###     of scenario 'scen_name' and varying the species [sp] to have abundances [abn]

def newModelSpectra(scen_name,linecenter,width,sp,abn,upper=0.,ddex=0.1):
    """
    Produce a spectrum over a wavelength range of [linecenter-width-2,linecenter+width+2],
    where the abundance of element 'sp' is set to 'abn'.
    """

    # Move to scenario folder
    os.chdir(f'tlusty/scenarios/{scen_name}/synspec/')

    # # Save original spectrum if not already done
    # if 'base.7' not in os.listdir():
    #     subprocess.run(['cp','fort.7','base.7'])

    ### Write input files
    # Write new fort.55 -- spectral range 
    unit55 = f"""       0      50       0
            1       0       0       1
            0       1       1       1       0
            1       1       0       0       0
            2       0       0
        {linecenter-width-5}    {linecenter+width+5}      10       0  1.e-12     0.01
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
    if len(params)==3:
        do_ism = False
        abn,vrad,logscale = params
    elif len(params)==5:
        do_ism = True
        abn,vrad,ism_width,ism_amp,logscale = params
    else:
        do_ism = False
        print("Invalid number of parameters. Exiting.")
        exit()
    
    ### Interpolate spectrum at the given abundance

    # Get upper and lower bounding abundances (accurate to one decimal places)
    ddex        = spec_dict['ddex']
    abn_lo      = max(abn - (ddex/2) + 1e-4,-20.)
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

    ### Apply scale factor
    flux_base *= 10**logscale

    ### Apply wavelength shift
    beta        = vrad / 3e5
    wvln_shift  = wvln_base * np.sqrt( (1 + beta) / (1 - beta) )

    ### Convolve with COS LSF or a Gaussian
    if type(conv_stuff)==tuple:
        # You are convolving this spectrum with the COS instrument LSF.
        wvln_conv,flux_conv = convolveModelWithCOS(wvln_shift,flux_base,files=conv_stuff)
    else:
        wvln_conv,flux_conv = convolveModelWithGaussian(wvln_shift,flux_base,width=conv_stuff)
    
    ### Restrict to wavelengths in the region of interest
    norm_range  = np.where(np.abs(wvln_conv-linecenter) <= width)
    wvln_norm   = wvln_conv[norm_range]
    flux_norm   = flux_conv[norm_range]

    ### Apply ISM line
    if do_ism:
        ism_factor  = 1 - ism_amp*np.exp(-0.5*(wvln_norm-ism_center)**2/ism_width**2)
        flux_final   = flux_norm * ism_factor
    else:
        flux_final = flux_norm

    ### Determine scale factor by equalizing median fluxes
    # med_obs = np.median(flux_obs)
    # med_mod = np.median(flux_final)
    # flux_final *= med_obs/med_mod

    return np.float64(wvln_norm),np.float64(flux_final)

def lineFit(fname_spec,scen,sp,linecenter,width,lower=-10.,upper=-4.,conv_stuff=0.,
            p0=[],ism_center=0.,known_lines=[],
            nfits:int=3,dwidth:float=4.,uplim_ratio:float=1000.,
            showfig:bool=True,savefig:bool=False):
    """
    
    """

    ### Generate spectra
    newModelSpectra(scen,linecenter,width+(dwidth*(nfits-1)),sp,lower,upper,0.1)

    ### Load spectral dictionary
    spec_dict = preload_spectra(scen,sp,lower,upper,ddex=0.1)

    ### Load in and normalize observed spectrum

    # Load in
    try:
        inst = fits.getheader(fname_spec)['INSTRUME']
        if conv_stuff!=0:
            wvln_obs_base,flux_obs_base,sigma_obs_base = openCOSfile(fname_spec)
        # elif 'xsh' in fname_spec:
        elif inst=='XSHOOTER':
            wvln_obs_base,flux_obs_base,sigma_obs_base ,gaussian = openXshooter(fname_spec)
            conv_stuff = 5000/gaussian/2.355
        # elif 'uves' in fname_spec:
        elif inst=='UVES':
            wvln_obs_base,flux_obs_base,sigma_obs_base ,gaussian = openUVES(fname_spec)
            conv_stuff = 5000/gaussian/2.355
        else:
            print("Unexpected data format. Exiting...")
            exit()
    except:
        dat = np.genfromtxt(fname_spec,delimiter=',')
        wvln_obs_base   = dat[:,0]*10
        flux_obs_base   = dat[:,1]
        sigma_obs_base  = dat[:,2]
        conv_stuff      = 1/2.355

    ### Do N fits at different wavelengths
    fit_abundances  = []
    uncertainties   = []
    vrad_values     = []
    sigmas          = []
    scales          = []
    fit_type        = "Detection"

    # Define parameter bounds and initial guess for least squares
    mean_abn    = np.mean([lower,upper])
    if ism_center != 0.:
        bounds = [
            (lower+1e-2,-150,0.01,0.01,-30),
            (upper-1e-2,150,0.1,1,0)
        ]
        if p0==[]:
            p0 = [mean_abn,40,0.03,0.3,-15]
    else:
        bounds = [
            (lower+1e-2,-150,-30),
            (upper-1e-2,150,0)
        ]
        if p0==[]:
            p0 = [mean_abn,40,-15]

    fig,ax = plt.subplots(ncols=nfits,nrows=2)

    for i in range(nfits):
    
        # Trim
        wrange      = np.where(np.abs(wvln_obs_base-linecenter) <= width)
        wvln_obs    = wvln_obs_base[wrange]
        flux_obs    = flux_obs_base[wrange]
        sigma_obs   = sigma_obs_base[wrange]

        ### Run least squares

        # Define residuals function
        def residuals(params):
            # Load in and normalize synthetic spectrum
            wvln_model,flux_model = forward_model(params,linecenter,width,spec_dict,conv_stuff,ism_center)
            flux_model = np.interp(wvln_obs,wvln_model,flux_model)

            # Scale up residuals
            return (flux_model - flux_obs)*1e30

        # Run it
        result  = least_squares(residuals,p0,method='trf',xtol=1e-10,bounds=bounds)
        popt    = result.x
        # print(popt)

        ### Probability distribution function
        other_params    = popt[1:]
        abn_range       = np.arange(lower,upper,0.01)
        scores          = np.zeros_like(abn_range)

        # Vary abundance
        for j,abn in enumerate(abn_range):
            params = [abn]
            for p in other_params: params.append(p)
            # input(params)
            wvln_model,flux_model = forward_model(params,linecenter,width,spec_dict,conv_stuff,ism_center)
            flux_model = np.interp(wvln_obs,wvln_model,flux_model)

            score = chi2(flux_obs,flux_model,sigma_obs)
            scores[j] = score

        scores      -= min(scores)
        probs       = np.exp(-0.5*scores)
        pdf_norm    = probs/sum(probs)
        cdf         = np.zeros_like(pdf_norm)
        for k in range(len(cdf)-1):
            cdf[k] = sum(pdf_norm[:k+1])/sum(pdf_norm)
        cdf[-1] = 1

        # Limits from CDF
        most_prob                   = abn_range[np.argmax(pdf_norm)]
        lower_1sigma,upper_1sigma   = np.interp([0.341,0.659],cdf,abn_range)
        lower_3sigma,upper_3sigma   = np.interp([0.01,0.99],cdf,abn_range)
        best_vrad   = other_params[0]
        best_scale  = other_params[-1]

        # Save results
        fit_abundances.append(most_prob)

        unc = np.mean(np.abs([lower_1sigma,upper_1sigma]-most_prob))
        uncertainties.append(unc)

        sigmas.append([lower_1sigma,upper_1sigma,lower_3sigma,upper_3sigma])
        scales.append(best_scale)
        vrad_values.append(best_vrad)

        if np.max(pdf_norm)/pdf_norm[0] < uplim_ratio:
            fit_type = "Upper limit"

        ### Plot results
        if showfig or savefig:

            ### Best-fit model
            params = [most_prob]
            for p in other_params: params.append(p)
            wvln_model,flux_model = forward_model(params,linecenter,width,spec_dict,conv_stuff,ism_center)
            flux_model = np.interp(wvln_obs,wvln_model,flux_model)
            score = chi2(flux_obs,flux_model,sigma_obs)

            # If applicable, add best-fitting ISM line (WIP)
            # if ism_center != 0.:
                # ism_width,ism_amp = popt[2:4]
                # line_mask   = 1 - ism_amp*np.exp(-0.5*(wvln_obs-ism_center)**2/(ism_width)**2)
                # flux_ism    = line_mask*flux_model
                # ax[0].plot(wvln_obs,flux_ism,c='red',linestyle='dashed',zorder=100)

            ax[0,i].plot(wvln_obs,flux_model,c='red',zorder=100,label=r"Best fit, $\chi^2=$"+str(np.around(score,2)))
            ax[0,i].plot(wvln_obs,flux_obs,c='black',label="Observed Spectrum")
            ax[0,i].set_xlabel(r"Wavelength [$\AA$]",size='large')
            ax[0,i].set_ylabel("Flux",size='large')
            # ax[0].set(xlim=(np.min(wvln_obs),np.max(wvln_obs)))
            # ax[0,i].legend()

            ax[1,i].plot(abn_range,pdf_norm,c='black')
            ax[1,i].set_xlabel("Abundance",size='large')
            ax[1,i].set_ylabel("Relative probability",size='large')
            ax[1,i].vlines(most_prob,ymin=0,ymax=max(pdf_norm)*1.1,colors='grey',linestyle='solid',zorder=-100,
                        label="Best fit: "+np.format_float_positional(most_prob,precision=2))
            ax[1,i].vlines([lower_1sigma,upper_1sigma],ymin=0,ymax=max(pdf_norm)*1.1,colors='grey',linestyle='dashed',zorder=-100,
                        label=r"$1\sigma$ limits: +"+np.format_float_positional(upper_1sigma-most_prob,precision=2)
                        +', '+np.format_float_positional(lower_1sigma-most_prob,precision=2))
            ax[1,i].vlines([lower_3sigma,upper_3sigma],ymin=0,ymax=max(pdf_norm)*1.1,colors='grey',linestyle='dotted',zorder=-100,
                        label=r"$3\sigma$ limits: +"+np.format_float_positional(upper_3sigma-most_prob,precision=2)
                        +', '+np.format_float_positional(lower_3sigma-most_prob,precision=2))
            # ax[1].vlines(upper_3sigma,ymin=0,ymax=max(pdf_norm)*1.1,colors='grey',linestyle='dashed',zorder=-100,
            #             label=r"$3\sigma$ upper limit: "+np.format_float_positional(upper_3sigma,precision=2))
            
            ax[1,i].set_xlim(
                max(np.nan_to_num(most_prob + (upper_3sigma - most_prob)*2.5,nan=-10),lower),
                min(np.nan_to_num(most_prob + (lower_3sigma - most_prob)*2.5,nan=0),upper)
            )
            ax[1,i].invert_xaxis()
            ax[1,i].set_ylim(0,np.nanmax(pdf_norm)*1.1)
            # ax[1,i].legend()

        width += dwidth

    fig.set_figwidth(11)
    plt.tight_layout()
    if savefig:
        plt.savefig(f'outputs/linefit_{scen}_elem{sp}_lc{linecenter}_w{width}.png',dpi=200)
    # plt.legend()
    if showfig:
        plt.show()
    plt.close()

    # Averaged results
    ave_abn     = np.mean(fit_abundances)
    ave_sig_abn = np.mean(uncertainties)
    sig_width   = np.std(fit_abundances)
    sig         = np.sqrt(sig_width**2 + ave_sig_abn**2)

    print(ave_abn)
    print(ave_sig_abn)
    print(sig_width)
    print(sig)

    return ave_abn,sig,ave_sig_abn,sig_width

def testLineFit():

    # FUV observed spectrum to fit
    fname_spec = "data/fits_files/gd56_cos_fuv.fits"
    conv_stuff = loadConvolutionFiles(fname_spec)
    # fname_spec = "data/fits_files/Gaia3296_xsh_nuv.fits"
    # conv_stuff = 0.

    # Scenario folder. If not present, run it
    if 'test_gd56' not in os.listdir(SCEN):
        subprocess.run(['python3','tlusty.py','test_gd56','15000','8','chondrite=-10','synspec=1'])
    scen = 'test_gd56'

    # # Line parameters   -- Si II lines at 1260 and 1265 AA
    linecenter  = 1262
    width       = 7
    sp          = 14    
    ism         = 0.
    lower       = -10.
    upper       = -4.

    # Line parameters   -- C II lines at 1335 AA
    # linecenter  = 1335
    # width       = 8
    # sp          = 6    
    # ism         = 1334.53
    # lower       = -10.
    # upper       = -4.

    # Line parameters   -- Ca II line, 3933 AA
    # linecenter  = 3935
    # width       = 5
    # sp          = 20
    # ism         = 0.
    # lower       = -10.
    # upper       = -5.

    # Line parameters   -- Mg II lines, 4481 AA
    # linecenter  = 4481
    # width       = 5
    # sp          = 12
    # ism         = 0.
    # lower       = -7
    # upper       = -4

    res = lineFit(fname_spec=fname_spec,
                  scen=scen,
                  linecenter=linecenter,
                  width=width,
                  sp=sp,
                  lower=lower,
                  upper=upper,
                  conv_stuff=conv_stuff,
                  ism_center=ism,
                  showfig=True,
                  nfits=5,dwidth=2)

testLineFit()
input()

if __name__=='__main__':
    # Because this function needs a fair but of info to run, the user is prompted for each
    # parameter.

    ################################################################################

    ### Option to automate line fitting
    line_file = input("If you have a list of lines to fit, please type the path to the list. Otherwise,\npress enter:\t")
    if line_file != '':

        # Read line fit file
        f = open(line_file,'r')
        lines = f.readlines()
        f.close()

        for line in lines:

            line = line.split()

            if len(line) <= 1:
                continue

            if line[0][0] == '#': continue

            # Read spectrum name, scenario name
            if len(line) == 2:
                fname_spec,scen = line

                # Check fname_spec
                try:
                    f = fits.open(fname_spec)
                    f.close()
                except:
                    print("File not found. Exiting.")
                    exit()

                # print("Checking spectrum file...")
                instr   = fits.getheader(fname_spec)['INSTRUME']
                if instr=='COS':
                    conv_stuff = loadConvolutionFiles(fname_spec)
                else:
                    conv_stuff = 0.

                # Check scen
                if scen not in os.listdir('tlusty/scenarios/'):
                    print(f"Scenario \"{scen}\" not found. Exiting.")
                    exit()

            else:
                # Open results file
                out = open(f"outputs/linefit_{scen}.log",'a')

                sp_str,lower,upper,linecenter,width,ism_center = line[:6]
                if len(line) > 6:
                    pass

                # Check species
                try:
                    sp = int(sp_str)
                except:
                    sp = ElementtoZ(sp_str)
                    if sp==-1:
                        print("Element not supported: "+sp_str)
                        exit()

                # Format lower and upper
                lower = np.around(float(lower),2)
                upper = np.around(float(upper),2)

                # Check linecenter
                linecenter = float(linecenter)
                if (linecenter < 900) or (linecenter > 25000):
                    print("Wavelength is out of range. Exiting")
                    exit()

                # Format width
                try:
                    width = float(width)
                except ValueError:
                    width = 5.

                # Format ISM line
                ism_center  = float(ism_center)

                # print("Generating range of spectra...")
                newModelSpectra(scen,linecenter,width,sp,lower,upper,0.1)
                spectra = preload_spectra(scen,sp,lower,upper,0.1)

                res = lineFit(
                            fname_spec=fname_spec,
                            linecenter=linecenter,
                            width=width,
                            scen=scen,
                            sp=sp,
                            lower=lower,
                            upper=upper,
                            conv_stuff=conv_stuff,
                            ism_center=ism_center,
                            p0=[-7,40],
                            known_lines=[],
                            showfig=False,
                            savefig=True
                )

                out.write('\t'.join([str(sp),str(linecenter),str(width),res['type'],str(res['value'])])+'\n')
                out.close()


    ################################################################################

    else:
        ### The name of the FITS file containing the observed spectrum
        fname_spec  = input("Please provide the path to the FITS spectrum to be fit:\t")
        try:
            f = fits.open(fname_spec)
            f.close()
        except:
            print("File not found. Exiting.")
            exit()

        ### The scenario that will be used to generate spectra
        scen        = input("Please provide the name of the scenario to use as a base:\t")
        if scen not in os.listdir('tlusty/scenarios/'):
            print("Scenario not found. Exiting.")
            exit()

        linecenter  = input("Please provide the central wavelength of the line/region to be fit (AA):\t")
        linecenter  = float(linecenter)
        if (linecenter < 900) or (linecenter > 25000):
            print("Wavelength is out of range. Exiting")
            exit()

        width       = input("Please provide the width of the region to be fit, in Angstroms (default: 5AA):\t")
        try:
            width = float(width)
        except ValueError:
            width = 5.

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



