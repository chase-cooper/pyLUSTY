import os,shutil,subprocess,sys,time
from utils import *

####################
###    Models    ###
####################

def pureHmodel(scen:str,teff:float,log_g:float):
    """
    Run a simple, LTE-LTGray model for a pure-H (DA) white dwarf.

    Inputs:
        - scen (string):        the name of the folder where results are stored.
        - teff (float):         the effective temperature of the white dwarf in Kelvin
        - log_g (float):        the surface gravity of the white dwarf, in cgs units

    Outputs:
        - run_finished (bool):  whether or not the model converged successfully (WIP)
    
    """

    ntries  = 10
    
    # Make directories and link files
    os.chdir(PATH)
    if scen in os.listdir('scenarios'):
        ans = input(f'Directory {PATH}/{scen}/ already exists. Would you like to replace it? [y/n] ')
        if ans.lower()=='n':
            print('Exiting...')
            sys.exit(0)
        elif ans.lower()=='y':
            print("Removing directory...")
            shutil.rmtree(f'scenarios/{scen}')
        else:
            print('Invalid response. Exiting...')
            sys.exit(0)
    os.chdir('scenarios')
    os.mkdir(scen)
    os.chdir(scen)
    os.mkdir('lte')
    os.chdir('lte')
    subprocess.run(['ln','-s','-f',PATH+'/data','data'])
    subprocess.run(['ln','-s','-f',FLAG,'cwd.flag'])

    writeUnit5(path=PATH+'scenarios/'+scen+'/lte',teff=teff,log_g=log_g,is_lte=True)

    model_broke = False
    for i in range(ntries):
        try:
            os.system(f"{TLEXE} < fort.5 > log.dat")
            # subprocess.Popen([TLEXE,'<','fort.5','>','lte.6'],shell=True,stdout=subprocess.PIPE,text=True).wait()

            out6    = open(PATH+f'scenarios/{scen}/lte/log.dat','r').read()
            out7    = open(PATH+f'scenarios/{scen}/lte/fort.7','r').read()
            if ('FINAL MODEL ATMOSPHERE' not in out6):
                print(f"Attempt {i+1}: Model did not finish.")
                model_broke = True
            elif ('NaN' in out7):
                print(f"Attempt {i+1}: NaNs in output structure.")
                model_broke = True
            else:
                model_broke = False
                break
        except:
            print(f"Attempt {i+1}:   Model run has encountered an error.")
            model_broke = True

    # Rename output files
    subprocess.run(['cp','fort.7','structure.dat'])
    subprocess.run(['rm','-f','fort.7'])

    if model_broke:
        print("Model did not successfully converge. Check 'log.dat' and 'structure.dat' for details.")
        run_finished = False
    else:
        print("LTE model completed!")
        run_finished    = True

    os.chdir(HOME)
    return run_finished
    
def metals(scen:str,teff:float,log_g:float,spdict:dict,is_lte:bool=True):
    """
    Produce an LTE WD model with metal lines. Requires an LTE model as a base.

    Inputs:
    - scen (str):       the name of the scenario where the pure-H model output is saved.
    - teff (float):     effective temperature
    - log_g (float):    surface gravity
    - spdict (dict):    species dictionary containing abundances for metals
    """

    # Make subdirectory and link files
    os.chdir(PATH)
    os.chdir('scenarios')
    try:
        os.chdir(scen)
    except:
        print("Scenario does not exist. Exiting...")
        sys.exit()
    os.mkdir('metals')
    os.chdir('metals')
    subprocess.run(['ln','-s','-f',PATH+'/data','data'])
    subprocess.run(['ln','-s','-f',FLAG,'cwd.flag'])

    # Write unit 5 file
    writeUnit5(path=os.getcwd(),
               teff=teff,
               log_g=log_g,
               species_dict=spdict,
               is_lte=is_lte,
               is_ltgray=False)         # LTGRAY=False means a base atmosphere file is needed

    # Copy output from pure-H model
    subprocess.run(['cp','../lte/structure.dat','fort.8'])

    ntries=10
    model_broke = False
    for i in range(ntries):
        try:
            os.system(f"{TLEXE} < fort.5 > log_metals.dat")
            out7    = open(PATH+f'scenarios/{scen}/metals/log_metals.dat','r').read()
            if ('FINAL MODEL ATMOSPHERE' not in out7) or ('NaN' in out7):
                print(f"Attempt {i+1}:  NaN's in output.")
                model_broke = True
            else:
                model_broke = False
                break
        except:
            print(f"Attempt {i+1}:   Model run has encountered an error.")
            model_broke = True

    # Rename output files
    subprocess.run(['cp','fort.7','structure.dat'])
    subprocess.run(['rm','-f','fort.7'])

    if model_broke:
        print("Model did not successfully converge. Check 'log.dat' and 'structure.dat' for details.")
        run_finished = False
    else:
        print("LTE metals model completed!")
        run_finished    = True

    os.chdir(HOME)
    return run_finished

def synspec(scen:str,species_dict,imode:int=0,lammin:int=900,lammax:int=10000):
    """
    
    """

    # Write input file
    unit5 = """       {1}      50       0
       1       0       0       0
       0       1       1       1       0
       1       1       0       0       0
       2       0       0
    {2}     {3}      10       0  0.0001     0.1
       0       0"""
    unit5 = unit5.replace('{1}',str(imode))
    unit5 = unit5.replace('{2}',str(lammin))
    unit5 = unit5.replace('{3}',str(lammax))

    os.chdir(PATH)
    os.mkdir('scenarios/'+scen+'/synspec')
    os.chdir('scenarios/'+scen+'/synspec')
    f = open('fort.55','w')
    f.write(unit5)
    f.close()

    subprocess.run(['ln','-s','-f',PATH+'/data','data'])
    subprocess.run(['ln','-s','-f',FLAG,'cwd.flag'])


    # imode determines whether or not we use a linelist
    if imode==0:
        subprocess.run(['cp',PATH+'/data/gfATO.dat','fort.19'])
        # subprocess.run(['cp',PATH+'/data/gfATO.dat',SPEC_LTE+'/fort.19'])
        
        # Remove lines for species not present in the model
        elem_ids={
            1: 'h',
            2: 'he',
            3: 'li',
            4: 'be',
            5: 'b',
            6: 'c',
            7: 'n',
            8: 'o',
            9: 'fl',
            10:'ne',
            11:'na',
            12:'mg',
            13:'al',
            14:'si',
            15:'p',
            16:'s',
            17:'cl',
            18:'ar',
            19:'k',
            20:'ca',
            21:'sc',
            22:'ti',
            23:'v',
            24:'cr',
            25:'mg',
            26:'fe',
            27:'co',
            28:'ni'
        }

        f = open('fort.19','r')
        lines = f.read().split('\n')
        f.close()

        f = open('fort.19','w')
        for line in lines[:-1]:
            try:
                elem = int(line[13:15])
                if species_dict[elem_ids[elem]]['mode']==2:
                    f.write(line+'\n')
            except KeyError:    # Line associated with an element not yet implemented; ignore line
                continue

    # fetch outputs from NLTE line run
    if 'metals' in os.listdir('..'):
        subprocess.run(['cp','../metals/structure.dat','fort.8'])
        subprocess.run(['cp','../metals/fort.5','fort.5'])
    else:
        subprocess.run(['cp','../lte/structure.dat','fort.8'])
        subprocess.run(['cp','../lte/fort.5','fort.5'])
        
    # run SYNSPEC
    try:
        os.system(f'{SYNEXE} < fort.5 > fort.6')
    except:
        print("\n Model run has encountered an error. Exiting...")
        sys.exit(0)
    print('SYNSPEC complete!')

    # Rename spectrum file
    subprocess.run(['cp','fort.7',scen+'_spectrum.dat'])
    subprocess.run(['cp','fort.7','../'+scen+'_spectrum.dat'])
    subprocess.run(['rm','-f','fort.7'])

    os.chdir(PATH)

###################################################################################################

def main():   # or main

    start = time.time()
    # process user inputs
    try:
        scen    = sys.argv[1]
        teff    = float(sys.argv[2])
        log_g   = float(sys.argv[3])
    except IndexError:
        print("Not enough arguments supplied. Call function as:\n    % python tlusty.py [name] [t_eff] [log_g]")
        sys.exit(0)
    
    # Dictionary that holds species treatment modes and abundances. No functionality for changing partition functions
    species_dict = newSpeciesDict()

    # SYNSPEC parameters
    do_synspec  = False
    imode       = 2         # means a pure-hydrogen run
    lammin      = 900
    lammax      = 25000

    # Check for additional user arguments
    if len(sys.argv) > 4:
        # check for keywords
        for arg in sys.argv[4:]:
            
            key = arg.split('=')[0]
            val = arg.split('=')[1]

            # metal species changes
            if key in species_dict.keys():
                # Implicit or explicit treatment
                if val[-1]!='!':
                    species_dict[key]['mode'] = 2
                else:
                    species_dict[key]['mode'] = 1
                
                # Abundances
                if val[-1]=='!':
                    val = val[:-1]
                if float(val)<0:
                    species_dict[key]['abn'] = 10**float(val)
                else:
                    species_dict[key]['abn'] = float(val)
                imode=0

            # setting chondrite-like abundances
            if key=='chondrite':
                imode = 0

                if float(val)<0:
                    scale = 10**float(val)
                else:
                    scale = float(val)
                    
                # update species abundances according to chondrite relative abundances
                sp_to_update    = ['c','n','o','na','mg','al','si','s','ca','cr','mn','fe','ni']
                rel_abn_si = {
                    'c':    0.7724,
                    'n':    0.0554,
                    'o':    7.552,
                    'na':   0.05747,
                    'mg':   1.04,
                    'al':   0.08308,
                    'si':   1,
                    's':    0.4449,
                    'ca':   0.05968,
                    'cr':   0.01313,
                    'mn':   0.00917,
                    'fe':   0.8632,
                    'ni':   0.0478
                }

                for s in sp_to_update:
                    # species_dict[s]['abn'] = np.format_float_scientific(scale*rel_abn_si[s],precision=5)
                    species_dict[s]['abn']  = str(scale*rel_abn_si[s])
                    if s in ['cr','mn','ni']:
                        species_dict[s]['mode'] = 1
                    else:
                        species_dict[s]['mode'] = 2

            # imode for SYNSPEC
            if key=='synspec':
                do_synspec = bool(val)

            # SYNSPEC lower bound
            if key=='lammin':
                lammin = int(val)

            # SYNSPEC upper bound
            if key=='lammax':
                lammax = int(val)

    # Initial LTE model
    # WARNING: do not do a pure He model at all, can't validate that it is correct
    success = pureHmodel(scen=scen,teff=teff,log_g=log_g)

    if success and (species_dict != newSpeciesDict()):
        metals(scen=scen,teff=teff,log_g=log_g,spdict=species_dict,is_lte=True)

    if success and do_synspec:
        synspec(scen=scen,species_dict=species_dict,imode=imode,lammin=lammin,lammax=lammax)

    end = time.time()
    print(f"Time elapsed:   {round(end-start,3)}s")

if __name__ == '__main__':
    main()
    shrinkScenarios()
