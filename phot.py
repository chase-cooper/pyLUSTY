import matplotlib.pyplot as plt
import numpy as np
import requests,sys
from utils import *

### Catalogs
from astroquery.gaia import Gaia
from astroquery.ipac.irsa import Irsa
from astroquery.mast import Catalogs
from astroquery.mast import Mast
from astroquery.sdss import SDSS
from astroquery.simbad import Simbad
mast = Mast()

### Dust maps
from dustmaps.sfd import SFDQuery
from dustmaps.config import config
config['data_dir'] = HOME+'data/dustmaps'
sfd     = SFDQuery()

### MCMC Fitting
import corner,emcee

# Constants
c           = 3e8           # [m/s]
pc2m        = 3.086e16
pc2cm       = 3.086e18
grav        = 6.67e-11      # [m^3 kg^-1 s^-2]
msol        = 1.989e30      # [kg]
R_earth     = 6.371e6       # [m]
R_earth_cm  = R_earth*100   # [cm]

# GALEX dta file path
GALEX   = 'data/galex.dat'

################################################################################

### Catalog query functions

def query_simbad(name:str):

    ### Get designations from SIMBAD
    """
    Query simbad to get the Gaia DR3 designation for a given object.
    """

    name_gaia    = None

    # query SIMBAD
    res = Simbad.query_objectids(name)
    for des in res[:]:
        des = des[0]
        if 'Gaia DR3' in des:
            name_gaia = des

    # if no Gaia counterpart is found, quit
    if name_gaia==None:
        print("No Gaia ID found. Exiting...")
        exit(0)

    name_gaia = name_gaia.split()[-1]
    return name_gaia

def query_gaia(name:str,verbose:bool=False):

    ### Query Gaia main table for parallax, position, and photometry
    job = Gaia.launch_job("SELECT TOP 100 source_id, parallax, parallax_error, phot_g_mean_mag, phot_g_mean_flux, phot_g_mean_flux_error, phot_bp_mean_mag, phot_bp_mean_flux, phot_bp_mean_flux_error, phot_rp_mean_mag, phot_rp_mean_flux, phot_rp_mean_flux_error, b "
                                f"from gaiadr3.gaia_source where source_id = {name}",
                                dump_to_file=False)
    res = job.get_results()
    

    # Parallax and position
    plx         = res['parallax'][0]
    plx_err     = res['parallax_error'][0]
    b           = res['b'][0]   # Galactic latitude

    # Distance [pc]
    dist        = 1000/plx
    dist_err    = [dist-1000/(plx+plx_err),dist-1000/(plx-plx_err)]

    # Height above galactic plane
    z = dist*np.sin(np.deg2rad(b))
    z = np.abs(z)

    # Extinction fraction
    if dist < 100:
        ext_frac = 0
    else:
        ext_frac    = min( (dist-100) / (250/np.abs(np.sin(np.deg2rad(b))) - 100), 1)

    ### Query SDSS crossmatch table for SDSS neighbour
    job_sdss    = Gaia.launch_job(f"SELECT TOP 10 source_id, original_ext_source_id, angular_distance FROM gaiadr3.sdssdr13_best_neighbour WHERE source_id = {name}",
                                dump_to_file=False)
    sdss    = job_sdss.get_results()
    try:
        id_sdss = sdss['original_ext_source_id'][0]
    except IndexError:
        print("No SDSS DR13 neighbour found")
        id_sdss = None

    ### Query SkyMapper crossmatch table for SkyMapper neighbour
    job_skym    = Gaia.launch_job(f"SELECT TOP 10 source_id, original_ext_source_id, angular_distance FROM gaiadr3.skymapperdr2_best_neighbour WHERE source_id = {name}",
                                dump_to_file=False)
    skym    = job_skym.get_results()
    try:
        id_skym = skym['original_ext_source_id'][0]
    except IndexError:
        print("No SkyMapper DR2 neighbour found")
        id_skym = None

    ### Query PanSTARRS crossmatch table for PanSTARRS neighbour
    job_pans    = Gaia.launch_job(f"SELECT TOP 10 source_id, original_ext_source_id, angular_distance FROM gaiadr3.panstarrs1_best_neighbour WHERE source_id = {name}",
                                dump_to_file=False)
    pans    = job_pans.get_results()
    try:
        id_pans = pans['original_ext_source_id'][0]
    except IndexError:
        print("No PanSTARRS DR1 neighbour found")
        id_pans = None

    ### Query 2MASS crossmatch table for 2MASS neighbour
    job_tmas    = Gaia.launch_job(f"SELECT TOP 10 source_id, original_ext_source_id, angular_distance FROM gaiadr3.tmass_psc_xsc_best_neighbour WHERE source_id = {name}",
                                  dump_to_file=False)
    tmas    = job_tmas.get_results()
    try:
        id_2mas = tmas['original_ext_source_id'][0]
    except IndexError:
        print("No 2MASS neighbour found")
        id_2mas = None

    if verbose:
        print("*"*100)
        print("Parallax:        "+str(plx)+" +\\-"+str(plx_err)+" arcsec")
        print("Distance:        "+str(dist)+" +"+str(dist_err[0])+"\\-"+str(np.abs(dist_err[1]))+" pc")
        print("Gal. Latitude:   "+str(b)+" degrees")
        print("|z|:             "+str(z)+" pc")
        print("Ext. fraction:   "+str(ext_frac))
        print("\nNeighbours:")
        print(" - SDSS:         "+str(id_sdss))
        print(" - PanSTARRS:    "+str(id_pans))
        print(" - SkyMapper:    "+str(id_skym))
        print(" - 2MASS:        "+str(id_2mas))

    return plx,plx_err,ext_frac,id_sdss,id_pans,id_skym,id_2mas

def query_sdss(name:int,pos):

    # Query
    fields  = ['objID','u','g','r','i','z','err_u','err_g','err_r','err_i','err_z','extinction_u','extinction_g','extinction_r','extinction_i','extinction_z']
    xid     = SDSS.query_region(pos, radius='5 arcsec', spectro=False, fields=fields)
    
    # Get data
    for i in xid:
        # Check each row to see if it contains the SDSS designation from Gaia cross-matching
        if i['objID']==name:
            res = np.array(list(i))
    
            # Apply SDSS-to-AB mag corrections
            res[1] -= 0.04
            res[4] += 0.015
            res[5] += 0.03
            
            return res[1:6],res[6:11],res[11:]  # AB mags, errors, extinctions
    
    return None

def query_panstarrs(name:int,pos):
    try:
        catalog_data = Catalogs.query_region(pos,
                                        radius=0.01,
                                        catalog="Panstarrs",
                                        data_release="dr1",
                                        table="mean")
        
        # Find line whose object ID is the same one returned by Gaia cross-matching
        for line in catalog_data:
            if line['objID']==name:
                pans_mags   = np.array([line['gMeanApMag'],line['rMeanApMag'],line['iMeanApMag'],line['zMeanApMag'],line['yMeanApMag']])
                pans_errs   = np.array([line['gMeanApMagErr'],line['rMeanApMagErr'],line['iMeanApMagErr'],line['zMeanApMagErr'],line['yMeanApMagErr']])
                
                # Extinctions
                # pans_rvals       = np.array([3.384,2.482,1.93,1.551,1.258])        # values of R, from Laura
                pans_rvals       = np.array([3.172,2.271,1.682,1.322,1.087])       # values of R assuming R_V=3.1, Schlafly & Finkbeiner 2011

                return pans_mags,pans_errs,pans_rvals
    
    except:
        print("\nFailed to query PanSTARRS for object "+str(name_pans))
        return np.array([]),np.array([])
    
def query_skymapper(name,pos):

    # Query parameters
    ra,dec = pos.to_string(style='decimal').split()
    params = {
        'RA':   ra,
        'DEC':  dec,
        'SR':   '0.1',
        'VERB': '3',
        'RESPONSEFORMAT':   'CSV'
    }

    # make query
    url = 'https://skymapper.anu.edu.au/sm-cone/public/query?'
    for key in params.keys():
        url  += key + '=' + params[key] + '&'
    url = url[:-1]
    req     = requests.get(url)

    # save data
    txt     = req.text
    txt     = txt.split('\n')
    for line in txt:
        line = line.split(',')

        # Find line with correct name
        if line[1] == str(name):
            
            # Important indices:
            #   1  --> object id
            #   43 --> u mag
            #   57 --> g mag
            #   64 --> r mag
            #   71 --> i mag
            #   78 --> z mag
            #   95 --> closest Gaia neighbor

            skym_mags   = np.array([line[43],line[57],line[64],line[71],line[78]],dtype=np.float32)
            skym_errs   = np.array([line[44],line[58],line[65],line[72],line[79]],dtype=np.float32)

            # R-values
            # skym_rvals  = np.array([5,3.4,2.5,1.8,1.4])     # from Laura
            skym_rvals  = np.array([5.075,3.407,2.685,2.03,1.618])      # from Huang et al. 2019
            
            return skym_mags,skym_errs,skym_rvals
    
    # If object isn't found
    return np.array([]),np.array([]),np.array([])

def query_2mas(name,pos):
    try:
        catalog_data = Irsa.query_region(coordinates=pos,
                                         catalog='fp_psc')
        for line in catalog_data:
            if line['designation']==name:
                tmas_mags = np.array([line['j_m'],line['h_m'],line['k_m']])
                tmas_errs = np.array([line['j_cmsig'],line['h_cmsig'],line['k_cmsig']])
                tmas_rvals = np.array([0.72,0.46,0.306])    # from Yuan, Liu, and Xiang (2013)

        return tmas_mags,tmas_errs,tmas_rvals
        # return np.array(tmas_mags[0]),np.array(tmas_errs[0]),np.array(tmas_rvals[0])
        # return np.array([]),np.array([]),np.array([])
    except:
        print("\nFailed to query 2MASS for object "+name)
        return np.array([]),np.array([]),np.array([])

def load_galex(name):   # in testing, worried about the data itself
    """
    Retrieve GALEX photometry for an object given its GAIA DR3 designation. Only possible for
    likely white dwarfs (pwd > 0.75 in Gentile-Fusillo et al. 2021). Data kindly provided by 
    Mukremin Kilic.

    Parameters:
    - name (string):    the object's GAIA DR3 designation, either with or without the leading
                            'GAIA DR3' string.
    """

    # Get just the GAIA DR3 designation number
    desig   = name.split()[-1]

    # Load GALEX table
    file = open(GALEX,'r')
    galex_data = file.readlines()
    file.close()

    # If object is present, return GALEX magnitudes and uncertainties
    for line in galex_data:
        line = line.split()
        if line[1] != desig: continue

        # Table data. If not present (999), return empty arrays
        fuv     = float(line[2])
        fuv_err = float(line[3])
        nuv     = float(line[4])
        nuv_err = float(line[5])
        if '999' in [fuv,fuv_err,nuv,nuv_err]:
            break

        mags    = np.array([fuv,nuv])
        errs    = np.array([fuv_err,nuv_err])

        # Extinction correction (Yuan, Liu, and Xiang 2013)
        # rvals   = np.array([4.37,7.06])

        # Extinction correction (Zhang and Yuan 2023)
        # rvals   = np.array([6.97,7.26])

        # Extinction correction (Wall et al. 2019)
        rvals   = np.array([8.01,6.72])


        return mags,errs,rvals

    # Nothing found
    return np.array([]),np.array([]),np.array([])

################################################################################

### Retrieve all extinction-corrected photometric fluxes

def photometry(name:str):
    """
    Given an object's Gaia DR3 designation, this function:
        1) queries MAST for the object's RA and Dec, then uses the dustmaps package to get the
            line-of-sight E(B-V) using maps from SFD98 and the 14% decrease correction following
            SF11.
        2) queries Gaia for the object's height above the galactic plane z and parallax, then
            uses that and its position to calculate the fraction ofextinction to apply, following
            Harris (2006). It also queries the Gaia cross-matching tables for the corresponding
            SDSS, PanSTARRS, SkyMapper, and 2MASS designations.
        3) retrieves photometry, if present, from the SDSS, PanSTARRS, SkyMapper, 2MASS, and GALEX*
            catalogues, and dereddens them. Reddening values are provided in the individual query
            functions and their sources cited therein.
    
    Parameters:
    - name (string): the Gaia DR3 designation of the object, ex. "Gaia DR3 XXXXXXXXXXXXXX"

    Outputs:
    - all_wvlns (np.array(float)):      the effective wavelengths of bandpasses with successfully
            retrieved photometry
    - all_fluxes (np.array(float)):     the retrieved photometric fluxes in [Jy]
    - all_fluxerrs (np.array(float)):   the uncertainties in the retrieved photometric fluxes,
            in [Jy]. Uncertainties are derived from magnitude uncertainties, which are symmetric
            and produce assymmetric uncertainties in flux, so the upper uncertainties in flux,
            which are greater, are returned.
    """

    ########################
    ###  Object queries  ###
    ########################

    ### Get object position (RA, DEC) and get line-of-sight E(B-V)
    pos     = mast.resolve_object(name,resolver=None)
    ebv     = sfd(pos) * 0.86
    # print("E(B-V):      "+str(ebv))


    ### Query SIMBAD for object's designation in Gaia DR3. If no such designation exists, exit.
    print("Querying SIMBAD...",end='')
    name_gaia   = query_simbad(name)
    print("Done!")
    print('\n' + '*'*80 + '\n')


    ### Query Gaia for parallax and galactic latitude to derive distances and height above/below
    ###     the galactic plane |z|; use these to determine fraction of extinction to apply 
    ###     following procedure in Harris (2006).
    ###
    ###     Additionally, query Gaia's crossmatching tables for object designations in SDSS DR12,
    ###     PanSTARRS DR1, and SkyMapper DR2 surveys.
    print("Querying Gaia DR3...")
    # name_gaia   =   name.split()[-1]
    plx,plx_err,ext_frac,name_sdss,name_pans,name_skymapper,name_2mass = query_gaia(name_gaia,verbose=True)
    print("\nFinished querying Gaia DR3")
    print("\n" + '*'*80 + '\n')


    ### Query SDSS DR19 if an SDSS object ID is returned from Gaia query. Returns object mags,
    ###     errors, and extinctions for each ugriz passband.
    if name_sdss:
        print("Querying SDSS...")
        sdss_res    = query_sdss(name_sdss,pos)
        if sdss_res!=None:
            sdss_mags,sdss_errs,sdss_exts   = sdss_res
            sdss_mags   = np.array(sdss_mags)
            sdss_errs   = np.array(sdss_errs)
            sdss_exts   = np.array(sdss_exts)

            # Correct SDSS mags for reddening
            sdss_mags   -= sdss_exts*ext_frac

            # Magnitude error lower bound from Bergeron+2019
            sdss_errs   = np.maximum(sdss_errs,0.03)

            print(sdss_mags)
            print(sdss_errs)
            print(sdss_exts)
        
        else:
            # If no SDSS object ID is found, create empty arrays for SDSS quantities
            sdss_mags = np.array([])
            sdss_errs = np.array([])
            sdss_exts = np.array([])

        print("Done!")
        print('\n' + '*'*80 + '\n')
    
    else:
        # If no SDSS object ID is found, create empty arrays for SDSS quantities
        sdss_mags = np.array([])
        sdss_errs = np.array([])
        sdss_exts = np.array([])


    ### Query PanSTARRS DR1 if a PanSTARRS object ID is returned from Gaia query. Returns object
    ###     magnitudes, errors, and R-values for each grizy passband. R-values are provided in the
    ###     function call.
    if name_pans:
        print("Querying PanSTARRS...")
        pans_mags,pans_errs,pans_rvals = query_panstarrs(name_pans,pos) # already numpy arrays

        # Extinction correction
        pans_exts   = pans_rvals * ebv
        pans_mags   -= pans_exts * ext_frac

        # Magnitude error lower bound from Bergeron+2019
        pans_errs   = np.maximum(pans_errs,0.03)

        print(pans_mags)
        print(pans_errs)
        print(pans_exts)

        print("Done!")
        print('\n' + '*'*80 + '\n')

    else:
        # If no PanSTARRS object ID is found, create empty arrays for PanSTARRS quantities
        pans_mags = np.array([])
        pans_errs = np.array([])
        pans_exts = np.array([])
    

    ### Query SkyMapper DR2 if a SkyMapper object ID is returned from the Gaia query. Returns
    ###     object magnitudes, errors, and R-values for each ugriz passband. R-values are provided
    ###     in the function call.
    if name_skymapper:
        print("Querying SkyMapper DR2...")

        skym_mags,skym_errs,skym_rvals  = query_skymapper(name_skymapper,pos)
        skym_exts   = skym_rvals * ebv
        skym_mags   -= skym_exts * ext_frac

        # Magnitude error lower bound from Bergeron+2019
        skym_errs   = np.maximum(skym_errs,0.03)

        print(skym_mags)
        print(skym_errs)
        print(skym_exts)

        print("Done!")
        print('\n' + '*'*80 + '\n')
    
    else:
        # If no SkyMapper object ID is found, create empty arrays for SkyMapper quantities
        skym_mags = np.array([])
        skym_errs = np.array([])
        skym_exts = np.array([])

    ### Query 2MASS if a 2MASS object ID is returned from the Gaia query. Returns object
    ###     magnitudes, errors, and R-values for each jhk passband. R-values are provided in the
    ###     function call.
    if name_2mass:
        print("Querying 2MASS...")

        tmas_mags,tmas_errs,tmas_rvals = query_2mas(name_2mass,pos)
        tmas_exts = tmas_rvals*ebv
        tmas_mags -= tmas_exts

        # Magnitude error lower bound from Bergeron+2019
        tmas_errs   = np.maximum(tmas_errs,0.03)
        tmas_errs   = np.nan_to_num(tmas_errs,nan=0.03)

        print(tmas_mags)
        print(tmas_errs)
        print(tmas_exts)
        print("Done!")
        print('\n' + '*'*80 + '\n')

    else:
        tmas_mags = np.array([])
        tmas_errs = np.array([])
        tmas_exts = np.array([])

    ### Load GALEX photometry table. Returns object magnitudes, errors, and R-values for both
    ###     passbands. R-values are provided in the function call
    print("Loading GALEX data...")

    galex_mags,galex_errs,galex_rvals = load_galex(name)
    if np.any(galex_mags==999) or np.any(galex_errs==999):
        galex_mags = np.array([])
        galex_errs = np.array([])
        galex_exts = np.array([])
    else:
        galex_exts = galex_rvals*ebv
        galex_mags -= galex_exts

        # Magnitude error lower bound from Bergeron+2019
        galex_errs   = np.maximum(galex_errs,0.03)

    print(galex_mags)
    print(galex_errs)
    print(galex_exts)

    # # Data is dubious, turning off for now
    # galex_mags = np.array([0,0])

    print("Done!")
    print('\n' + '*'*80 + '\n')

################################################################################

    ### Formatting of data and conversion of magnitudes to fluxes in Jy
    ###
    ### This kinda sucks but whatever it works

    # Filter effective wavelengths [nm] from SVO
    sdss_eff_wvlns  = np.array([360.8,467.18,614.11,745.79,892.28]) * 1e-9
    pans_eff_wvlns  = np.array([481.02,615.55,750.3,866.84,961.36]) * 1e-9
    skym_eff_wvlns  = np.array([350.02,501.6,607.69,773.28,912.03]) * 1e-9
    tmas_eff_wvlns  = np.array([1235,1662,2159]) * 1e-9
    galex_eff_wvlns = np.array([154.89,230.34]) * 1e-9

    # Mags to fluxes and their errors
    all_wvlns       = np.array([])
    all_fluxes      = np.array([])
    all_fluxerrs_hi = np.array([])
    all_fluxerrs_lo = np.array([])
    
    if np.any(sdss_mags):   
        sdss_fluxes     = ABmag2flux(sdss_mags,sdss_eff_wvlns)
        sdss_fluxes_hi  = ABmag2flux(sdss_mags - sdss_errs,sdss_eff_wvlns)
        sdss_fluxes_lo  = ABmag2flux(sdss_mags + sdss_errs,sdss_eff_wvlns)
        sdss_fluxerr_hi = sdss_fluxes_hi - sdss_fluxes
        sdss_fluxerr_lo = sdss_fluxes - sdss_fluxes_lo
    else:
        sdss_fluxes     = np.array([0,0,0,0,0])
        sdss_fluxerr_hi = np.array([1e10,1e10,1e10,1e10,1e10])
        sdss_fluxerr_lo = np.array([1e10,1e10,1e10,1e10,1e10])

    if np.any(pans_mags):
        pans_fluxes     = ABmag2flux(pans_mags,pans_eff_wvlns)
        pans_fluxes_hi  = ABmag2flux(pans_mags - pans_errs,pans_eff_wvlns)
        pans_fluxes_lo  = ABmag2flux(pans_mags + pans_errs,pans_eff_wvlns)
        pans_fluxerr_hi = pans_fluxes_hi - pans_fluxes
        pans_fluxerr_lo = pans_fluxes - pans_fluxes_lo
    else:
        pans_fluxes     = np.array([0,0,0,0,0])
        pans_fluxerr_hi = np.array([1e10,1e10,1e10,1e10,1e10])
        pans_fluxerr_lo = np.array([1e10,1e10,1e10,1e10,1e10])

    if np.any(skym_mags):
        skym_fluxes     = ABmag2flux(skym_mags,skym_eff_wvlns)
        skym_fluxes_hi  = ABmag2flux(skym_mags - skym_errs,skym_eff_wvlns)
        skym_fluxes_lo  = ABmag2flux(skym_mags + skym_errs,skym_eff_wvlns)
        skym_fluxerr_hi = skym_fluxes_hi - skym_fluxes
        skym_fluxerr_lo = skym_fluxes - skym_fluxes_lo
    else:
        skym_fluxes     = np.array([0,0,0,0,0])
        skym_fluxerr_hi = np.array([1e10,1e10,1e10,1e10,1e10])
        skym_fluxerr_lo = np.array([1e10,1e10,1e10,1e10,1e10])

    if np.any(tmas_mags):   
        # NOTE: 2MASS mags are not absolute (AB) magnitudes; fluxes are calculated using zero-
        #   point fluxes in Cohen+03
        tmas_fluxes     = np.array([1594,1024,666.7]) * 10**(-0.4*tmas_mags)
        tmas_fluxes_lo  = np.array([1594,1024,666.7]) * 10**(-0.4 * (tmas_mags+tmas_errs) )
        tmas_fluxes_hi  = np.array([1594,1024,666.7]) * 10**(-0.4 * (tmas_mags-tmas_errs) )
        tmas_fluxerr_hi = tmas_fluxes_hi - tmas_fluxes
        tmas_fluxerr_lo = tmas_fluxes - tmas_fluxes_lo
    else:
        tmas_fluxes     = np.array([0,0,0])
        tmas_fluxerr_hi = np.array([1e10,1e10,1e10])
        tmas_fluxerr_lo = np.array([1e10,1e10,1e10])

    if np.any(galex_mags):
        galex_fluxes    = ABmag2flux(galex_mags,galex_eff_wvlns)
        galex_fluxes_hi = ABmag2flux(galex_mags-galex_errs,galex_eff_wvlns)
        galex_fluxes_lo = ABmag2flux(galex_mags+galex_errs,galex_eff_wvlns)
        galex_fluxerr_hi = galex_fluxes_hi - galex_fluxes
        galex_fluxerr_lo = galex_fluxes - galex_fluxes_lo
    else:
        galex_fluxes    = np.array([0,0])
        galex_fluxerr_hi = np.array([1e10,1e10])
        galex_fluxerr_lo = np.array([1e10,1e10])

    all_wvlns       = np.concatenate([sdss_eff_wvlns,pans_eff_wvlns,skym_eff_wvlns,tmas_eff_wvlns,galex_eff_wvlns])
    all_fluxes      = np.concatenate([sdss_fluxes,pans_fluxes,skym_fluxes,tmas_fluxes,galex_fluxes])
    all_fluxerrs_hi = np.concatenate([sdss_fluxerr_hi,pans_fluxerr_hi,skym_fluxerr_hi,tmas_fluxerr_hi,galex_fluxerr_hi])
    all_fluxerrs_lo = np.concatenate([sdss_fluxerr_lo,pans_fluxerr_lo,skym_fluxerr_lo,tmas_fluxerr_lo,galex_fluxerr_lo])

    return all_wvlns,all_fluxes,all_fluxerrs_hi,plx,plx_err

################################################################################

### Convolve a model spectrum with a filter response function

def convolveModelWithBandpass(mwvln,mflux,band_fname):
    """
    Convolve a TLUSTY/SYNSPEC spectrum (at location 'model_fname') with a filter response function
    (at location 'band_fname').
    """

    # Speed of light in AA/s
    cAA = 3e18

    # Datatype change
    mwvln       = np.array(mwvln)
    mflux       = np.array(mflux)

    model_wvln0 = mwvln[1:]
    model_flux0 = mflux[1:]
    model_dlam0 = mwvln[1:] - mwvln[:-1]

    # TLUSTY/SYNSPEC spectra are Eddington fluxes, so require a factor of 4*pi to get fluxes
    model_flux0 *= 4*np.pi

    filt    = np.genfromtxt(band_fname)
    filt_wvlns  = filt[:,0]
    filt_resp   = filt[:,1]

    # Interpolate response to model wavelengths
    model_response  = np.interp(model_wvln0,filt_wvlns,filt_resp,left=0,right=0)

    # Resulting flux in [erg/s/cm^2/Hz]
    conv_flux   = (1/cAA) * np.sum(model_response*model_flux0*model_wvln0*model_dlam0) / np.sum(model_response * model_dlam0 / model_wvln0)
    conv_flux   *= 1e23     # Convert to Jy
    return conv_flux

### Convolve model spectrum with bandpasses to get synthetic photometry

def model_flux_to_phot(mwvln,mflux):
    """
    
    """
    sdss_model_flux = np.array([])
    pans_model_flux = np.array([])
    skym_model_flux = np.array([])
    tmas_model_flux = np.array([])
    galex_model_flux= np.array([])

    for f in ['u','g','r','i','z']:
        pass_fname  = f"data/filters/sdss/{f}.dat"
        model_flux  = convolveModelWithBandpass(mwvln,mflux,pass_fname)
        sdss_model_flux  = np.concatenate([sdss_model_flux,np.array([model_flux])])
        
    for f in ['g','r','i','z','y']:
        pass_fname  = f"data/filters/panstarrs/{f}.dat"
        model_flux  = convolveModelWithBandpass(mwvln,mflux,pass_fname)
        pans_model_flux = np.concatenate([pans_model_flux,np.array([model_flux])])
        
    for f in ['u','g','r','i','z']:
        pass_fname  = f"data/filters/skymapper/{f}.dat"
        model_flux  = convolveModelWithBandpass(mwvln,mflux,pass_fname)
        skym_model_flux = np.concatenate([skym_model_flux,np.array([model_flux])])

    for f in ['j','h','k']:
        pass_fname = f"data/filters/2mass/{f}.dat"
        model_flux = convolveModelWithBandpass(mwvln,mflux,pass_fname)
        tmas_model_flux = np.concatenate([tmas_model_flux,np.array([model_flux])])

    for f in ['fuv','nuv']:
        pass_fname = f"data/filters/galex/{f}.dat"
        model_flux = convolveModelWithBandpass(mwvln,mflux,pass_fname)
        galex_model_flux = np.concatenate([galex_model_flux,np.array([model_flux])])

    return np.concatenate([sdss_model_flux,pans_model_flux,skym_model_flux,tmas_model_flux,galex_model_flux])

##################################################################################################################

### MCMC

def log_likelihood(theta,x,y,y_err,spec_dict):
    """
    Log likelihood function
    """

    teff,logg,r,plx   = theta
    d = (1000 / plx) * pc2m
    rm = r*R_earth

    # Forward model
    wvln,flux       = interp_spectrum_dict(teff,logg,spec_dict)
    phot_model_emi  = model_flux_to_phot(wvln,flux)
    phot_model_obs  = phot_model_emi*rm*rm/d/d

    # Chi2
    score       = chi2(y,phot_model_obs,y_err)
    if np.nan_to_num(score,nan=-1)==-1:
        print(y)
        print(phot_model_obs)
        print(y_err)
        print("FUCKFUCKFUCKFUCKFUCKFUCK")

    return  -0.5*score

def log_prior(theta,plx0,plx0_err):           # contains constants dependent on the model grid bounds
    """
    Log prior function
    """

    teff,logg,r,plx   = theta

    # distance prior from gaia data

    # log prior -- containts hard-coded constants set by model grid
    if (teff < 10000)   or (teff > 25000)\
    or (logg < 6.0)    or (logg > 9.0)\
    or (r < 0.2)        or (r > 3):
        return -np.inf

    else:
        return -0.5 * ( (plx - plx0) / plx0_err )**2

def log_prob(theta,x,y,y_err,plx0,plx0_err,spec_dict):
    """
    Log probability
    """
    lp = log_prior(theta,plx0,plx0_err)
    if not np.isfinite(lp):
        return -np.inf
    else:
        return lp + log_likelihood(theta,x,y,y_err,spec_dict)

def photometric_fit(name,wvln_obs,flux_obs,err_obs,plx0,plx0_err,spec_dict,
                    ndim:int=4,nwalkers:int=100,nstep:int=200):

    # # Set parallax and error therein for prior probability calc
    # global P0;      P0 = plx
    # global sig_P0;  sig_P0 = plx_err

    # initial guess of parameters
    x0  = np.array([15000,8,1,10])
    pos = [x0*(1 + 1e-2*np.random.randn((4))) for w in range(nwalkers)]

    # results file
    fn = f'outputs/emcee_phot_{name[-4:]}.h5'
    backend = emcee.backends.HDFBackend( fn )
    backend.reset( nwalkers , ndim )

    # The big show
    sampler = emcee.EnsembleSampler(nwalkers,ndim,log_prob,args=(wvln_obs,flux_obs,err_obs,plx0,plx0_err,spec_dict),backend=backend)

    # Regularly calculate autocorrelation time; if N > 50*tau, you're done
    ac_step = 100
    for i in range(0,nstep,ac_step):
        pos = sampler.run_mcmc(pos,ac_step,progress=True,skip_initial_state_check=True)
        autocor = sampler.get_autocorr_time(quiet=True)
        if i >= ac_step*max(autocor):
            print(f"Finished after {i+ac_step} steps")
            break
    # sampler.run_mcmc(pos,nstep,progress=True,skip_initial_state_check=True)
    # autocorrelation_time = sampler.get_autocorr_time()
    # print(autocorrelation_time)

def get_photometric_results(name,nburn,wvln_obs=[],flux_obs=[],err_obs=[],
                            spec_dict:dict={},plotit:bool=False):

    # Load spectra if not already loaded
    if spec_dict=={}:
        spec_dict = loadSpectra()

    # Get results
    # reader  = emcee.backends.HDFBackend(f'outputs/emcee_phot_{name[-4:]}.h5')
    reader  = emcee.backends.HDFBackend(f'outputs/emcee_phot_{name[-4:]}.h5')
    nthin   = 2
    samples = reader.get_chain(discard=nburn,thin=nthin,flat=True)

    # Get distributions of temperature and R
    t_mean = np.percentile(samples[:,0],50)
    t_sig  = np.percentile(samples[:,0],[15.9,84.1]) - t_mean
    r_mean = np.percentile(samples[:,2],50) * R_earth_cm
    r_sig  = np.percentile(samples[:,2],[15.9,84.1]) * R_earth_cm - r_mean

    # Estimate distribution of log(g) from Bédard models
    logg_values = [
        interpBedardLogg(t_mean+t_sig[0],r_mean+r_sig[0]),
        interpBedardLogg(t_mean+t_sig[1],r_mean+r_sig[0]),
        interpBedardLogg(t_mean,         r_mean+r_sig[0]),
        interpBedardLogg(t_mean+t_sig[0],r_mean+r_sig[1]),
        interpBedardLogg(t_mean+t_sig[1],r_mean+r_sig[1]),
        interpBedardLogg(t_mean,         r_mean+r_sig[1]),
        interpBedardLogg(t_mean+t_sig[0],r_mean),
        interpBedardLogg(t_mean+t_sig[1],r_mean),
        interpBedardLogg(t_mean,         r_mean)
    ]
    imp_logg = logg_values[-1]
    imp_logg_sig = np.std(logg_values)

    f = open('log.dat','a')
    f.write(name+'\tphot\t'+str(t_mean)+'\t'+str(t_sig[0])+'\t'+str(t_sig[1])+'\t')
    f.write(str(imp_logg)+'\t'+str(imp_logg_sig)+'\t')
    f.write(str(r_mean)+'\t'+str(r_sig[0])+'\t'+str(r_sig[1])+'\n')
    f.close()

    ### Plot results

    # Only plot points with T within 5 sigma
    t_sig_mean = np.mean(np.abs(t_sig))
    valid   = np.where(np.abs(samples[:,0] - t_mean) <= 3*t_sig_mean)
    samples = samples[valid]
    
    # MCMC
    fig     = corner.corner(samples,quantiles=[0.16,0.5,0.84],show_titles=True,labels=["T",r"$\log(g)$","R",r"$\pi$"])#,range=[[13000,25000],[7.75,8.25],[3e6,1.5e7],[1,100]])

    plt.savefig(f"outputs/corner_phot_{name[-4:]}.png",dpi=250)
    if plotit:
        plt.show()
    plt.close()

    # Best-fit photometry
    if (len(wvln_obs) > 0) and (len(flux_obs) > 0) and (len(err_obs) > 0):

        fig2,ax2    = plt.subplots()

        d = (1000 / plx) * pc2m
        rm = r_mean / 100

        # photometric points
        wvln,flux       = interp_spectrum_dict(t_mean,imp_logg,spec_dict)
        phot_model_emi  = model_flux_to_phot(wvln,flux)
        phot_model_obs  = phot_model_emi*rm*rm/d/d

        # model spectrum
        wvln_nm     = wvln/10
        flux        *= 4*np.pi*rm*rm/d/d
        flux_conv   = flux * 3.3356e4*wvln*wvln

        model_label = r"Model: $T_{\rm eff}="+str(int(t_mean))\
            +"K,log(g)="+np.format_float_positional(imp_logg,precision=2,min_digits=2)+"$"
        ax2.plot(wvln_nm,flux_conv,c='grey',zorder=-100,label=model_label)
        ax2.scatter(wvln_obs[flux_obs!=0]*1e9,phot_model_obs[flux_obs!=0],c='black',s=20,label="Best-fit photometry")
        ax2.errorbar(wvln_obs[flux_obs!=0]*1e9,flux_obs[flux_obs!=0],err_obs[flux_obs!=0],
                        ecolor='red',capsize=5,linestyle='',label='Observed')
        ax2.set_xlabel("Wavelength [nm]",size='large')
        ax2.set_ylabel("Flux [Jy]",size='large')
        ax2.set_yscale('log')
        plt.legend()
        plt.xlim(0.9*min(wvln_obs[flux_obs!=0]*1e9),1.1*max(wvln_obs[flux_obs!=0]*1e9))
        plt.ylim(0.3*min(flux_obs[flux_obs!=0]),3.1*max(flux_obs[flux_obs!=0]))
        plt.tight_layout()
        plt.savefig(f'outputs/photometry_{name[-4:]}.png',dpi=250)
        if plotit:
            plt.show()
        plt.close()

    return t_mean,np.average(np.abs(t_sig)),imp_logg,imp_logg_sig

##################################################################################################################

if __name__ == '__main__':

    # Function keywords
    name    = sys.argv[1]
    name    = name.strip('"')

    ndim        = 4
    nwalkers    = 100
    nstep       = 200
    nburn       = 100

    if len(sys.argv) > 2:
        keywords = sys.argv[2:]
        for k in keywords:
            k = k.split('=')

            if k[0]=='nwalkers':
                nwalkers = int(k[1])
            if k[0]=='nstep':
                nstep = int(k[1])
            if k[0]=='nburn':
                nburn = int(k[1])

    # Get photometry and parallax
    wvln_obs,flux_obs,err_obs,plx,plx_err = photometry(name)
    print(wvln_obs)
    print(flux_obs)
    print(err_obs)
    print(plx)
    print(plx_err)

    # Spectrum dictionary
    spec_dict = loadSpectra(shrink_factor=5)

    photometric_fit(name,wvln_obs,flux_obs,err_obs,plx,plx_err,spec_dict,
                    ndim,nwalkers,nstep)

    get_photometric_results(name,nburn=nburn,
                            wvln_obs=wvln_obs,flux_obs=flux_obs,err_obs=err_obs,
                            spec_dict=spec_dict,plotit=True)



