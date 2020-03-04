import momi        ## momi2 analysis
import logging        ## create log file
import numpy as np
import datetime
import matplotlib as plt

print("-----\n-----\n-----")

print("start logging\n")
logging.basicConfig(level=logging.INFO,
                    filename="momi_log_bs_asym.txt")

print("load sfs\n")
sfspath = "/home/kprovost/nas1/momi2/cardcard16_sfs_filtered_changelength_monomorphic.txt"
## this is a two-population sfs with monomorphic sites included in "length"
sfs = momi.Sfs.load(sfspath)
#print("Avg pairwise heterozygosity", sfs.avg_pairwise_hets[:5])
#print("populations", sfs.populations)
#print("percent missing data per population", sfs.p_missing)

print("\nPRIORS")
print("MUT RATE: 2.21e-9")
mutrate=2.21e-9

print("GEN TIME: 1")
gentime=1

print("ANCESTRAL NE: 300,000")
ancne=300000

print("DIV TIME RANGE: 500,000 to 2,500,000")
divtimelow=500000
divtimehigh=2500000

print("NE RANGE: 1,000 to 2,000,000")
nelow=1000
nehigh=2000000

print("MIGRATION RANGE: 0 to 0.5")
migratelow=0
migratehigh=0.5

print("MIGRATION DATE RANGE: 0 to 1,000,000\n")
migtimelow=0
migtimehigh=1000000

print("begin setting up models\n")
##### PURE ISOLATION MODEL #####
print("Pure Isolation model (base model)")
pure_isolation_model = momi.DemographicModel(N_e=ancne,muts_per_gen=mutrate,gen_time=gentime) ## why tho -- can you give it something?
pure_isolation_model.set_data(sfs)
## set up divergence times
pure_isolation_model.add_time_param("tdiv_sc",lower=divtimelow,upper=divtimehigh)
## set up effective population size
pure_isolation_model.add_size_param("ne_s",lower=nelow,upper=nehigh) ## this is from Brian's paper on cardinals
pure_isolation_model.add_size_param("ne_c",lower=nelow,upper=nehigh) 
## set up populations and phylogeny
pure_isolation_model.add_leaf("Son",N="ne_s")
pure_isolation_model.add_leaf("Chi",N="ne_c")
pure_isolation_model.move_lineages("Son", "Chi", t="tdiv_sc")
## randomize parameters and check them
#pure_isolation_model.set_params(randomize=True)
#print(pure_isolation_model.get_params())

## set up the rest of the models 

##### ASYMMETRIC MIGRATION #####
print("Asymmetric model (c2s as base)")
asym_model = pure_isolation_model.copy() ## copy isol
asym_model.add_pulse_param("mig_s2c",lower=migratelow,upper=migratehigh)
asym_model.add_pulse_param("mig_c2s",lower=migratelow,upper=migratehigh)
asym_model.add_time_param("tmig_asym",lower=migtimelow,upper=migtimehigh,upper_constraints=["tdiv_sc"])
asym_model.move_lineages("Chi","Son",t="tmig_asym",p="mig_s2c")
asym_model.move_lineages("Son","Chi",t="tmig_asym",p="mig_c2s")

## set parameters to be bootstrapped 
print("FINAL MODEL PARAMETERS:")
asym_model.set_params({'tdiv_sc': 2116497.645839919, 'ne_s': 187093.11311425938, 'ne_c': 411975.7314960601, 'mig_s2c': 0.2512902869760518, 'mig_c2s': 0.23942960914915856, 'tmig_asym': 437167.74218677887})
print(asym_model.get_params())

## Start bootstraps
## use stochastic then full likelihood 

n_bootstraps = 10
# make copies of the original models to avoid changing them
model = asym_model.copy()

bootstrap_results = []
for i in range(n_bootstraps):
    print("\n\n------")
    print("Fitting bootstrap "+str(i+1)+" out of "+str(n_bootstraps))

    print("Resampling data")
    # resample the data
    resampled_sfs = sfs.resample()
    # tell models to use the new dataset
    model.set_data(resampled_sfs)
    #add_pulse_copy.set_data(resampled_sfs)

    print("Randomizing submodel")
    # choose new random parameters for submodel, optimize
    model.set_params(randomize=True)
    
    now = datetime.datetime.now()
    print("Begin stochastic optimization"+str(now))
    model.stochastic_optimize(num_iters=10, n_minibatches=5, save_to_checkpoint="momi_checkpoint_priorset0_stoch_asym.txt", svrg_epoch=-1)
    
    now = datetime.datetime.now()
    print("Finished stochastic optimization"+str(now))
    print(model.get_params())

    now = datetime.datetime.now()
    print("Begin full optimization"+str(now))
    model.optimize(method="L-BFGS-B")
    
    now = datetime.datetime.now()
    print("Finished full optimization"+str(now))
    print(model.get_params())

    bootstrap_results.append(model.get_params())
    
    try:
        with open("cardinalis-bootstrap-results-asym.temp", 'a') as outfile:
            outfile.write("{}\n".format(bootstrap_results[i]))
    except:
        print("PROBLEM OPENING FILE FOR BOOTSTRAP"+str(i+1))

## Save bootstrappies
print("\n\n~~~~~\n")
print("Writing bootstraps to file -- CHECK FOR DUPLICATES!")
with open("cardinalis-bootstrap-results-asym.temp", 'a') as outfile:
    for bootstrap in bootstrap_results:
        outfile.write("{}\n".format(bootstrap))
