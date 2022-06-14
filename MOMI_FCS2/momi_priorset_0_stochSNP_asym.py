import momi        ## momi2 analysis
import logging        ## create log file
import numpy as np
import datetime

print("-----\n-----\n-----")
now = datetime.datetime.now()

print("start logging\n")
logging.basicConfig(level=logging.INFO,
                    filename="momi_log_priorset_0_stochSNP_asym.txt")

print("load sfs\n")
sfspath = "/home/kprovost/nas1/momi2/cardcard16_sfs_filtered_changelength_monomorphic.txt"
## this is a two-population sfs with monomorphic sites included in "length"
sfs = momi.Sfs.load(sfspath)
print("Avg pairwise heterozygosity", sfs.avg_pairwise_hets[:5])
print("populations", sfs.populations)
print("percent missing data per population", sfs.p_missing)

## set up two-population model with parameters
## use Pure Isolation model as base for all of the models 
## because only looking at this particular split, elected to change ranges to get 

print("\nPRIORS")
print("MUT RATE: 2.21e-9")
mutrate=2.21e-9

print("GEN TIME: 1")
gentime=1

print("ANCESTRAL NE: 300,000")
ancne=300000

print("DIV TIME RANGE: 800,000 to 1,000,000")
divtimelow=800000
divtimehigh=1000000

print("NE RANGE: 400,000 to 500,000")
nelow=400000
nehigh=500000

print("MIGRATION RANGE: 0 to 0.1")
migratelow=0
migratehigh=0.1

print("MIGRATION DATE RANGE: 100,000 to 200,000\n\n")
migtimelow=100000
migtimehigh=200000

print("begin setting up models\n")
##### PURE ISOLATION MODEL #####
print("\nPure Isolation model (base model)")
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
print("\nAsymmetric model (c2s as base)")
asym_model = pure_isolation_model.copy() ## copy isol
asym_model.add_pulse_param("mig_s2c",lower=migratelow,upper=migratehigh)
asym_model.add_pulse_param("mig_c2s",lower=migratelow,upper=migratehigh)
asym_model.add_time_param("tmig_asym",lower=migtimelow,upper=migtimehigh,upper_constraints=["tdiv_sc"])
asym_model.move_lineages("Chi","Son",t="tmig_asym",p="mig_s2c")
asym_model.move_lineages("Son","Chi",t="tmig_asym",p="mig_c2s")
## randomize and check parameters
asym_model.set_params(randomize=True)
print(asym_model.get_params())

## optimize each model once
print("#####")

models = [asym_model]
model_names = ["ASYM"]

AICs = []
count = 0
for model in models:
    now = datetime.datetime.now()
    name = str(model_names[count])
    print("Optimizing "+name+" model: "+str(now))
    model.stochastic_optimize(num_iters=10, snps_per_minibatch=1000, save_to_checkpoint="momi_checkpoint_priorset0_stochSNP_asym.txt", svrg_epoch=-1)
    now = datetime.datetime.now()
    print("Finished optimizing "+name+": "+str(now))
    print(model.get_params())
    print("Starting AIC likelihood for "+name)
    lik = model.log_likelihood()
    nparams = len(model.get_params())
    aic = 2*nparams - 2*lik
    print("AIC {}".format(aic))
    AICs.append(aic)
    count += 1
    print("-----")

minv = np.min(AICs)
delta_aic = np.array(AICs) - minv
print("Delta AIC per model: ", delta_aic)
print("AIC weight per model: ", np.exp(-0.5 * delta_aic))
    
## TODO: add searching multiple times
## TODO: add bootstrapping

