import momi		## momi2 analysis
import logging		## create log file
import numpy as np
import datetime
import matplotlib as plt

print("-----\n-----\n-----")

print("start logging\n")
logging.basicConfig(level=logging.INFO,
					filename="momi_log_bs_uc2s.txt")

print("load sfs\n")
sfspath = "/home/kprovost/nas1/momi2/cardcard16_sfs_filtered_changelength_monomorphic.txt"
## this is a two-population sfs with monomorphic sites included in "length"
sfs = momi.Sfs.load(sfspath)
#print("Avg pairwise heterozygosity", sfs.avg_pairwise_hets[:5])
#print("populations", sfs.populations)
#print("percent missing data per population", sfs.p_missing)

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
pure_isolation_model.set_params({'tdiv_sc': 807803.1164182945, 'ne_s': 228681.9631920871, 'ne_c': 591672.789892767})
print(pure_isolation_model.get_params())


print("\nUni Chi2Son model")
uni_chi_to_son_model = pure_isolation_model.copy() ## copies the Pure Isolation model 
## add additional parameters 
uni_chi_to_son_model.add_pulse_param("mig_c2s",lower=0,upper=0.1)
uni_chi_to_son_model.add_time_param("tmig_c2s",lower=0,upper=1000000,upper_constraints=["tdiv_sc"])
uni_chi_to_son_model.move_lineages("Son","Chi",t="tmig_c2s",p="mig_c2s")
## randomize parameters and check them
uni_chi_to_son_model.set_params({'tdiv_sc': 1288768.9782706778, 'ne_s': 182530.65255371892, 'ne_c': 572108.2654352491, 'mig_c2s': 0.08839532726042339, 'tmig_c2s': 127212.47479790074})
print(uni_chi_to_son_model.get_params())


## Make bootstrappies

n_bootstraps = 100
# make copies of the original models to avoid changing them
uni_chi_to_son_model_copy = uni_chi_to_son_model.copy()
#add_pulse_copy = add_pulse_model.copy()

bootstrap_results = []
for i in range(n_bootstraps):
	print("Fitting bootstrap "+str(i+1)+" out of "+str(n_bootstraps))

	# resample the data
	resampled_sfs = sfs.resample()
	# tell models to use the new dataset
	uni_chi_to_son_model_copy.set_data(resampled_sfs)
	#add_pulse_copy.set_data(resampled_sfs)

	# choose new random parameters for submodel, optimize
	uni_chi_to_son_model_copy.set_params(randomize=True)
	uni_chi_to_son_model_copy.optimize()
	# initialize parameters from submodel, randomizing the new parameters
	#add_pulse_copy.set_params(pulse_copy.get_params(),
							  #randomize=True)
	#add_pulse_copy.optimize()

	bootstrap_results.append(uni_chi_to_son_model_copy.get_params())

## Save bootstrappies
with open("cardinalis-bootstrap-results-uc2s.txt", 'w') as outfile:
	for bootstrap in bootstrap_results:
		outfile.write("{}\n".format(bootstrap))
