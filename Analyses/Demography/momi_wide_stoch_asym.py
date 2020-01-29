import momi        ## momi2 analysis
import logging        ## create log file
import numpy as np
import datetime

print("-----\n-----\n-----")
now = datetime.datetime.now()

print("start logging\n")
logging.basicConfig(level=logging.INFO,
                    filename="momi_wide_stoch_log_asym.txt")

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
print("ANCESTRAL NE: 300,000")
print("GEN TIME: 1")
print("DIV TIME RANGE: 500,000 to 2,500,000")
print("NE RANGE: 1,000 to 2,000,000")
print("MIGRATION RANGE: 0 to 0.5")
print("MIGRATION DATE RANGE: 0 to 1,000,000 [25,000 for SECC model]\n\n")


print("begin setting up models\n")
##### PURE ISOLATION MODEL #####
print("\nPure Isolation model (base model)")
pure_isolation_model = momi.DemographicModel(N_e=300000,muts_per_gen=2.21e-9,gen_time=1) ## why tho -- can you give it something?
pure_isolation_model.set_data(sfs)
## set up divergence times
pure_isolation_model.add_time_param("tdiv_sc",lower=500000,upper=2500000)
## set up effective population size
pure_isolation_model.add_size_param("ne_s",lower=1000,upper=2000000) ## this is from Brian's paper on cardinals
pure_isolation_model.add_size_param("ne_c",lower=1000,upper=2000000) 
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
asym_model.add_pulse_param("mig_s2c",lower=0,upper=0.5)
asym_model.add_pulse_param("mig_c2s",lower=0,upper=0.5)
asym_model.add_time_param("tmig_asym",lower=0,upper=1000000,upper_constraints=["tdiv_sc"])
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
    print("Stochastic optimizing "+name+" model: "+str(now))
    model.stochastic_optimize(num_iters=10, n_minibatches=5, save_to_checkpoint="momi_checkpoint_priorset3_stoch_asym.txt", svrg_epoch=-1)
    now = datetime.datetime.now()
    print("Finished stochastic optimizing "+name+": "+str(now))
    print(model.get_params())
    print("Starting AIC likelihood for stochastic "+name)
    lik = model.log_likelihood()
    nparams = len(model.get_params())
    aic = 2*nparams - 2*lik
    print("AIC {}".format(aic))
    AICs.append(aic)
    count += 1
    print("-----")

count = 0
for model in models:
	now = datetime.datetime.now()
	print("Fully optimizing "+name+" model: "+str(now))
	model.optimize(method="L-BFGS-B")
	now = datetime.datetime.now()
	print("Finished fully optimizing "+name+": "+str(now))
	print(model.get_params())
	print("Starting AIC likelihood for full "+name)
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

