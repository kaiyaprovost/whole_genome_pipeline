import momi		## momi2 analysis
import logging		## create log file
import numpy as np
import datetime

print("-----\n-----\n-----")

print("start logging\n")
logging.basicConfig(level=logging.INFO,
					filename="momi_log_RAND_isol.txt")

print("load sfs\n")
sfspath = "/home/kprovost/nas1/momi2/cardcard16_sfs_filtered_changelength_monomorphic.txt"
## this is a two-population sfs with monomorphic sites included in "length"
sfs = momi.Sfs.load(sfspath)
print("Avg pairwise heterozygosity", sfs.avg_pairwise_hets[:5])
print("populations", sfs.populations)
print("percent missing data per population", sfs.p_missing)

## from FSC: 
## Prior specifications were as follows: all effective population size estimates 
## were a log-uniform distribution with a range between 50,000 and 1,000,000 
## haploid individuals, all migration estimates were a log-uniform distribution 
## with a range between 0.001 and 20 individuals per generation, and all 
## divergence times were a uniform distribution with a range between 100,000 
## and 3,000,000 generations. In addition, we constrained the time of divergence 
## between the Sonoran and Chihuahuan deserts to be more recent than the time 
## of divergence between the desert populations of Northern Cardinal 
## and the C. c. carneus outgroup, and the time of secondary contact 
## to be more recent than either of those values. 

## set up two-population model with parameters
## use Pure Isolation model as base for all of the models 
## because only looking at this particular split, elected to change ranges to get 

print("\nPRIORS")
print("MUT RATE: 2.21e-9")
print("ANCESTRAL NE: 1,000,000")
print("GEN TIME: 1")
print("DIV TIME RANGE: 500,000 to 1,500,000")
print("NE RANGE: 1,000 to 2,000,000")
print("MIGRATION RANGE: 0 to 0.1")
print("MIGRATION DATE RANGE: 0 to 1,000,000 [25,000 for SECC model]\n\n")

print("begin setting up models\n")
##### PURE ISOLATION MODEL #####
print("\nPure Isolation model (base model)")
pure_isolation_model = momi.DemographicModel(N_e=1000000,muts_per_gen=2.21e-9,gen_time=1) ## why tho -- can you give it something?
pure_isolation_model.set_data(sfs)
## set up divergence times
pure_isolation_model.add_time_param("tdiv_sc",lower=500000,upper=1500000)
## set up effective population size
pure_isolation_model.add_size_param("ne_s",lower=1000,upper=2000000) ## this is from Brian's paper on cardinals
pure_isolation_model.add_size_param("ne_c",lower=1000,upper=2000000) 
## set up populations and phylogeny
pure_isolation_model.add_leaf("Son",N="ne_s")
pure_isolation_model.add_leaf("Chi",N="ne_c")
pure_isolation_model.move_lineages("Son", "Chi", t="tdiv_sc")
## randomize parameters and check them
pure_isolation_model.set_params(randomize=True)
print(pure_isolation_model.get_params())

## optimize each model once
print("#####")

models = [pure_isolation_model]
model_names = ["ISOL"]

AICs = []
n_runs = 10
count = 0
for model in models:
	copy = model.copy()
	results = []
	best_lik = []
	now = datetime.datetime.now()
	name = str(model_names[count])
	print("Testing "+name+" models for local optima: "+str(now))
	for i in range(n_runs):
		print(f"Starting run "+str(i+1)+" out of "+str(n_runs)+"...")
		now = datetime.datetime.now()
		copy.set_params(
		# parameters inherited from no_pulse_model are set to their previous values
		#no_migration_model.get_params(),
		# other parmaeters are set to random initial values
		randomize=True)
		copy.optimize(options={"maxiter":200})
		#results += copy.optimize(options={"maxiter":200})
		print("Run "+str(i+1)+" results:")
		print(copy.get_params())
		print("Get AIC")
		lik = copy.log_likelihood()
		nparams = len(copy.get_params())
		aic = 2*nparams - 2*lik
		print("AIC {}".format(aic))
		if best_lik == []:
			best_lik = [aic]
			results = copy.get_params()
			model.set_params(copy.get_params())
		if aic < best_lik:
			best_lik = [aic]
			results = copy.get_params()
			model.set_params(copy.get_params())
	now = datetime.datetime.now()		
	print("\n\nFinished testing "+name+": "+str(now))
	print(results)
	print("AIC {}".format(best_lik))
	print("Optimizing "+name+" model using these starting params: "+str(now))
	print(model.get_params())
	model.optimize()
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