import matplotlib.pyplot as plt
import numpy as np
import sys

### MCMC
import corner,emcee

### Other
from utils import *

### This file will specifically be for spectroscopic fitting with flat priors, since that
###     will most likely be the fitting method we use

##########################################################################################

def bilinear_interp(teff,logg,spec_dict):

    # From values, get adjacent grid points
    t_lo        = teff - (teff % 500)
    t_lo_str    = str(int(t_lo))
    t_hi        = t_lo + 500
    t_hi_str    = str(int(t_hi))
    diff_t      = teff - t_lo

    logg_lo     = logg - (logg % 0.5)
    logg_lo_str = np.format_float_positional(logg_lo,precision=2,min_digits=2)
    logg_hi     = logg_lo + 0.5
    logg_hi_str = np.format_float_positional(logg_hi,precision=2,min_digits=2)
    diff_logg   = logg - logg_lo

    # Retrieve spectra
    spec00      = spec_dict[f"t{t_lo_str}_g{logg_lo_str}"]
    spec01      = spec_dict[f"t{t_lo_str}_g{logg_hi_str}"]
    spec10      = spec_dict[f"t{t_hi_str}_g{logg_lo_str}"]
    spec11      = spec_dict[f"t{t_hi_str}_g{logg_hi_str}"]
    wvln        = spec00['wvln']

    # Interpolate
    spec0x      = spec00['flux']*(1-diff_logg/0.5) + np.interp(wvln,spec01['wvln'],spec01['flux'])*diff_logg/0.5
    spec1x      = spec10['flux']*(1-diff_logg/0.5) + np.interp(wvln,spec11['wvln'],spec11['flux'])*diff_logg/0.5
    spec_final  = spec0x*(1-diff_t/500) + spec1x*diff_t/500

    return wvln,spec_final
    
def normalizeBalmerLines(wvln,flux,sigma=[],fit_halpha:bool=False):

    lines = []

    ### Fit first 5 Balmer lines
    line_centers    = [6562.8,4861.4,4340.5,4101.7,3970.1,3889.1][(not fit_halpha):]
    widths          = [150,110,80,50,30,20][(not fit_halpha):]
    labels          = [r"H$\alpha$",r"H$\beta$",r"H$\gamma$",r"H$\delta$",r"H$\epsilon$",r"H$\zeta$"][(not fit_halpha):]

    for i in range(len(line_centers)):
        # normalize spectra to average fluxes at a fixed distance from the line center
        blue_ave        = np.average(flux[np.abs(wvln-(line_centers[i]-widths[i])) <= 2])
        red_ave         = np.average(flux[np.abs(wvln-(line_centers[i]+widths[i])) <= 2]) #
        cut_wvln        = wvln[np.abs(wvln-line_centers[i]) <= widths[i]]
        cut_flux_obs    = flux[np.abs(wvln-line_centers[i]) <= widths[i]]
        
        norm_obs        = np.interp(cut_wvln,
                                    [cut_wvln[0],cut_wvln[-1]],
                                    [blue_ave,red_ave])
                                    # [cut_flux_obs[0],cut_flux_obs[-1]])

        cut_flux_norm = cut_flux_obs/norm_obs
        line = [cut_wvln-line_centers[i],cut_flux_norm]

        if np.any(sigma):
            cut_sigma   = sigma[np.abs(wvln-line_centers[i]) <= widths[i]]
            cut_sigma_norm       = cut_sigma/cut_flux_obs
            line.append(cut_sigma_norm)

        lines.append(line)

    return lines

def fitBalmerLinesToModel(wvln_obs,flux_obs,sigma_obs,wvln_model,flux_model,fit_halpha:bool=False):

    # Interpolate model spectrum to observed wavelengths
    flux_model = np.interp(wvln_obs,wvln_model,flux_model)

    # Get Balmer lines for each spectrum
    lines_obs   = normalizeBalmerLines(wvln_obs,flux_obs,sigma_obs,fit_halpha=fit_halpha)
    lines_mod   = normalizeBalmerLines(wvln_obs,flux_model,fit_halpha=fit_halpha)

    # Observed normalized uncertainties
    all_obs_line_wvlns  = np.concat([line[0] for line in lines_obs])

###############################################################################

### MCMC functions

def log_likelihood(theta,wvln_obs,flux_obs,sigma_obs,width,spec_dict):

    # start = time.time()

    # Extract model parameters
    teff,logg,vrad = theta

    # Get model spectrum
    wvln_model,flux_model = interp_spectrum_dict(teff,logg,spec_dict)
    range_model = np.where((wvln_model-5000) <= 2000)
    wvln_model = wvln_model[range_model]
    flux_model = flux_model[range_model]
    # print("Interp spectrum: "+str(time.time()-start))

    # Wavelength shift
    beta = vrad/3e5
    wvln_model *= np.sqrt( (1+beta) / (1-beta))

    # Vacuum wavelength to air wavelength (Morton 1991)
    wvln_model = wvln_model / (1 + 2.735e-4 + 131.4182/(wvln_model**2) + 2.7625e8/(wvln_model**4) ) 

    # Interpolate to observed wavelengths
    flux_model = np.interp(wvln_obs,wvln_model,flux_model)

    # Convolve with Gaussian to simulate point spread
    _,flux_model = convolveModelWithGaussian(wvln_obs,flux_model,width)
    # print("Convolve:        "+str(time.time()-start))

    # Normalize Balmer lines
    lines_model = normalizeBalmerLines(wvln_obs,flux_model)
    # print("Norm model:      "+str(time.time()-start))

    # Normalize observed Balmer lines
    lines_obs = normalizeBalmerLines(wvln_obs,flux_obs,sigma_obs)
    # print("Norm observed:   "+str(time.time()-start))

    # Calculate chi2
    total_score = 0
    for i in range(len(lines_model)):
        wvln_obs_i,flux_obs_i,sigma_obs_i = lines_obs[i]
        wvln_model_i,flux_model_i  = lines_model[i]
        subscore = chi2(flux_obs_i,flux_model_i,sigma_obs_i)
        total_score += subscore
        # plt.plot(wvln_obs_i,flux_obs_i+0.3*i,c='black')
        # plt.plot(wvln_obs_i,flux_model_i+0.3*i,c='red')
    # print("Chi2:            "+str(time.time()-start))
    # input()

    return -0.5*total_score

def log_prior(theta,t_prior_lo,t_prior_hi):

    teff,logg,vrad = theta

    if  (teff < t_prior_lo) or (teff > t_prior_hi)\
    or  (logg < 6) or (logg > 9)\
    or (vrad < -200) or (vrad > 200):
        return -np.inf

    else:
        return 0

def log_prob(theta,wvln_obs,flux_obs,sigma_obs,width,t_prior_lo,t_prior_hi,spec_dict):

    logprior = log_prior(theta,t_prior_lo,t_prior_hi)
    if not np.isfinite(logprior):
        return -np.inf
    else:
        return logprior + log_likelihood(theta,wvln_obs,flux_obs,sigma_obs,width,spec_dict)
    
def spectroscopic_fit(name,wvln_obs,flux_obs,sigma_obs,width,t_prior_lo,t_prior_hi,spec_dict,
                        nwalkers,nstep):

    ### Run MCMC

    # Initial parameter guess
    x0  = np.array([np.mean([t_prior_lo,t_prior_hi]),8,10])
    pos = [x0*(1+1e-2*np.random.randn((3))) for w in range(nwalkers)]

    # backup file
    fn = f'outputs/emcee_spec_{name[-4:]}.h5'
    backend = emcee.backends.HDFBackend( fn )
    backend.reset( nwalkers , ndim )

    # The big show
    sampler = emcee.EnsembleSampler(nwalkers=nwalkers,ndim=ndim,
                                    log_prob_fn=log_prob,
                                    args=(wvln_obs,flux_obs,sigma_obs,width,t_prior_lo,t_prior_hi,spec_dict),
                                    backend=backend)

    autocor_step = 100
    for i in np.arange(0,nstep,autocor_step):
        pos = sampler.run_mcmc(pos,nsteps=autocor_step,progress=True,skip_initial_state_check=True)
        autocor = sampler.get_autocorr_time(quiet=True)
        if i >= autocor_step*max(autocor):
            print(f"Finished after {i+autocor_step} steps")
            break

def get_spectroscopic_results(name,nburn,wvln_obs,flux_obs,width,spec_dict):

    # Get results
    reader  = emcee.backends.HDFBackend(f'outputs/emcee_spec_{name[-4:]}.h5')
    nthin   = 2
    samples = reader.get_chain(discard=nburn,thin=nthin,flat=True)

    # Get distributions of temperature, log(g), and vrad
    t_mean      = np.percentile(samples[:,0],50)
    t_sig       = np.percentile(samples[:,0],[15.9,84.1]) - t_mean
    logg_mean   = np.percentile(samples[:,1],50)
    logg_sig    = np.percentile(samples[:,1],[15.9,84.1]) - logg_mean
    vrad_mean   = np.percentile(samples[:,2],50)
    vrad_sig    = np.percentile(samples[:,2],[15.9,84.1]) - vrad_mean

    f = open('log.dat','a')
    f.write(name+'\tspec\t'+str(t_mean)+'\t'+str(t_sig[0])+'\t'+str(t_sig[1])+'\t')
    f.write(str(logg_mean)+'\t'+str(logg_sig[0])+'\t'+str(logg_sig[1])+'\t')
    f.write(str(vrad_mean)+'\t'+str(vrad_sig[0])+'\t'+str(vrad_sig[1])+'\n')
    f.close()
    
    # MCMC
    fig     = corner.corner(samples,quantiles=[0.16,0.5,0.84],show_titles=True,labels=["T",r"$\log(g)$",r"$v_{\rm rad}$"])#,range=[[t_prior_lo,t_prior_hi],[6,9],[-200,20]])
    plt.savefig(f"outputs/corner_spec_{name[-4:]}.png",dpi=250)
    plt.close()

    # Best fit Balmer lines
    if (len(wvln_obs) > 0) and (len(flux_obs) > 0):
        lines_obs = normalizeBalmerLines(wvln_obs,flux_obs)

        print(t_mean,logg_mean,vrad_mean)
        wvln_model,flux_model = interp_spectrum_dict(t_mean,logg_mean,spec_dict)

        # Wavelength shift
        beta = vrad_mean/3e5
        wvln_model *= np.sqrt( (1+beta) / (1-beta))
    
        # Vacuum wavelength to air wavelength (Morton 1991)
        wvln_model = wvln_model / (1 + 2.735e-4 + 131.4182/(wvln_model**2) + 2.7625e8/(wvln_model**4) ) 
    
        # Interpolate to observed wavelengths
        flux_model = np.interp(wvln_obs,wvln_model,flux_model)
    
        # Convolve with Gaussian to simulate point spread
        _,flux_model = convolveModelWithGaussian(wvln_obs,flux_model,width)
        # print("Convolve:        "+str(time.time()-start))
        lines_model = normalizeBalmerLines(wvln_obs,flux_model)

        line_centers    = [6562.8,4861.4,4340.5,4101.7,3970.1,3889.1][1:]
        widths          = [150,110,80,50,30,20][1:]
        labels          = [r"H$\alpha$",r"H$\beta$",r"H$\gamma$",r"H$\delta$",r"H$\epsilon$",r"H$\zeta$"][1:]

        fig,ax = plt.subplots()
        for i in range(len(lines_obs)):
            dlam,flux = lines_obs[i]
            if i==0:
                ax.plot(dlam,flux+0.3*i,c='black',label=f"T={int(t_mean)}K, log(g)={np.round(logg_mean,2)}")
            else:
                ax.plot(dlam,flux+0.3*i,c='black')

            dlam,flux = lines_model[i]
            ax.plot(dlam,flux+0.3*i,c='red',zorder=100)

            ax.text(-widths[i]-15,1+0.3*i,labels[i])

        ax.legend()
        ax.set_xlim(-130,130)
        ax.set_xlabel(r"$\Delta\lambda [\AA]$",size='x-large')
        ax.set_ylabel("Normalized flux + offset",size='x-large')
        plt.tight_layout()
        plt.savefig(f"outputs/balmer_{name[-4:]}.png",dpi=250)

if __name__ == '__main__':
    # Get function keywords
    name = ''
    mode = 'flatprior'
    fname = ''
    prior = ''
    t_prior_lo = 10000
    t_prior_hi = 25000

    ndim        = 3
    nwalkers    = 100
    nstep       = 200
    nburn       = 100

    if len(sys.argv) > 1:
        keywords = sys.argv[1:]
        for arg in keywords:
            arg = arg.split('=')

            if arg[0] == 'name':
                name = arg[1]
            if arg[0] in ['file','fname']:
                fname = arg[1]
            if arg[0]=='mode':
                mode = arg[1]

            if arg[0]=='nwalkers':
                nwalkers = int(arg[1])
            if arg[0]=='nstep':
                nstep = int(arg[1])
            if arg[0]=='nburn':
                nburn = int(arg[1])

            if arg[0]=='prior':
                vals = arg[1].split(',')
                t_prior_lo = float(vals[0])
                t_prior_hi = float(vals[1])

    # Read in observed spectrum
    if 'xsh' in fname:
        wvln_obs,flux_obs,sigma_obs,res = openXshooter(fname)

    elif 'uves' in fname:
        # need two files
        wvln_obs_s1,flux_obs_s1,sigma_obs_s1,res = openUVES(fname)
        fname2      = fname[:-8]+"vis.fits"
        wvln_obs_s2,flux_obs_s2,sigma_obs_s2,_   = openUVES(fname2)
        valid2      = np.where(wvln_obs_s2 > max(wvln_obs_s1))
        wvln_obs    = np.concat([wvln_obs_s1,wvln_obs_s2[valid2]])
        flux_obs    = np.concat([flux_obs_s1,flux_obs_s2[valid2]])
        sigma_obs   = np.concat([sigma_obs_s1,sigma_obs_s2[valid2]])

    else:
        exit()

    width = 5000/res

    # #######
    spec_dict = loadSpectra()

    spectroscopic_fit(name,wvln_obs,flux_obs,sigma_obs,width,t_prior_lo,t_prior_hi,spec_dict,
                      nwalkers,nstep)

    get_spectroscopic_results(name,nburn,wvln_obs,flux_obs,width,spec_dict)


