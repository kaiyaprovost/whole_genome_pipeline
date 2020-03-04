import momi        ## momi2 analysis
import logging        ## create log file
import numpy as np
import datetime
import matplotlib as plt

print("-----\n-----\n-----")

print("start logging\n")
logging.basicConfig(level=logging.INFO,
                    filename="momi_log_pop3_secc.txt")

print("load sfs\n")
#sfspath = "/home/kprovost/nas1/momi2/cardcard16_sfs_filtered_changelength_monomorphic.txt"
sfspath = "/home/kprovost/nas1/momi2/cardcard16_sfs_monomorphic_fourpop.txt"

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

print("SON-CHI DIV TIME RANGE: 100,000 to 2,000,000")
divtimelow=100000
divtimehigh=2000000

print("SPLIT WITH CARNEUS (OUTGROUP): 2,000,000")
carnsplit=2000000

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
pure_isolation_model.add_size_param("ne_o",lower=nelow,upper=nehigh)
## set up populations and phylogeny
pure_isolation_model.add_leaf("Son",N="ne_s")
pure_isolation_model.add_leaf("Chi",N="ne_c")
pure_isolation_model.add_leaf("Out",N="ne_o")
pure_isolation_model.move_lineages("Son", "Chi", t="tdiv_sc")
pure_isolation_model.move_lineages("Chi", "Out", t=carnsplit)
## randomize parameters and check them
#pure_isolation_model.set_params(randomize=True)
#print(pure_isolation_model.get_params())

## set up the rest of the models 

##### seccETRIC MIGRATION #####
print("\nSecondary contact")
## constrain time of migration to be after LGM, ~25000 years or sooner
sec_contact_model = pure_isolation_model.copy() ## copy the Pure Isolation model
## add parameters
sec_contact_model.add_pulse_param("mig_secc",lower=0,upper=0.5)
sec_contact_model.add_time_param("tmig_secc",lower=0,upper=25000,upper_constraints=["tdiv_sc"])
sec_contact_model.move_lineages("Chi","Son",t="tmig_secc",p="mig_secc")
sec_contact_model.move_lineages("Son","Chi",t="tmig_secc",p="mig_secc")
## randomize and get parameters

## randomize parameters and check them
sec_contact_model.set_params(randomize=True)
print(sec_contact_model.get_params())

models = [sec_contact_model]
model_names = ["secc"]

AICs = []
count = 0
for model in models:
    now = datetime.datetime.now()
    name = str(model_names[count])
    print("Stochastic optimizing "+name+" model: "+str(now))
    model.stochastic_optimize(num_iters=10, n_minibatches=5, save_to_checkpoint="momi_3pop_checkpoint_secc.txt", svrg_epoch=-1)
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