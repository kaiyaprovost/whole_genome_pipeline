# coding: utf-8
import momi        ## momi2 analysis
import numpy as np
import datetime
import matplotlib as plt
import sys
import os
import glob
print("-----\n-----\n-----")

now = datetime.datetime.now()
print("")
print("Current date and time using str method of datetime object:")
print(str(now))

print("load sfs\n")

## set up a try-except block to take in the sfs
try:
    sfspath = sys.argv[1]
    ## this is a two-population sfs with monomorphic sites included in "length"
except:
    print("\nNo 2POP sfs given, quitting")
    sys.exit()
    
print("load priors\n")
try:
    mu=sys.argv[2]
except:
    print("\tMu not given, omitting mu")
    mu=0

try:
    ne_anc=sys.argv[3]
except:
    print("\tAncestral Ne not given, defaulting to 10,000")
    ne_anc=10000
    
try:
    gentime=sys.argv[4]
except:
    print("\tGen time not given, defaulting to 1")
    gentime=1
    
try:
    divrange_bot=sys.argv[5]
    divrange_top=sys.argv[6]
except:
    print("\tBoth divergence time ranges not given, defaulting to 1 to 1e11")
    divrange_bot=1
    divrange_top=1e11

try:
    ne_bot=sys.argv[7]
    ne_top=sys.argv[8]
except:
    print("\tBoth Ne ranges not given, defaulting to 1 to 1e11")
    ne_bot=1
    ne_top=1e11

try:
    mig_bot=sys.argv[9]
    mig_top=sys.argv[10]
except:
    print("\tBoth Migration rate ranges not given, defaulting to 0 to 1")
    mig_bot=0
    mig_top=1

try:
    migdate_bot=sys.argv[11]
    migdate_top=sys.argv[12]
except:
    print("\tBoth Migration time ranges not given, defaulting to 25,000 to 1e11 (not sec con only)")
    migdate_bot=25000
    migdate_top=1e11
    
try:
    seccon_bot=sys.argv[13]
    seccon_top=sys.argv[14]
except:
    print("\tSecondary contact bounds not given, defaulting to 1 to 25,000")
    seccon_bot=1
    seccon_top=25000

  
sfs = momi.Sfs.load(sfspath)
print("Avg pairwise heterozygosity 2D", sfs.avg_pairwise_hets[:5])
print("populations 2D", sfs.populations)
print("percent missing data per population 2D", sfs.p_missing)

sfs1Dpath = sfs.replace("geno","1D")
sfs1D = momi.Sfs.load(sfs1Dpath)
print("Avg pairwise heterozygosity 1D", sfs1D.avg_pairwise_hets[:5])
print("populations 1D", sfs1D.populations)
print("percent missing data per population 1D", sfs1D.p_missing)

print("begin setting up models\n")

##### NO SPLIT MODEL #####
print("\nOne pop model")
if (mu == 0):
     onepop_model = momi.DemographicModel(N_e=ne_anc, gen_time=gentime) 
else:
     onepop_model = momi.DemographicModel(N_e=ne_anc, muts_per_gen=mu, gen_time=gentime) 

onepop_model.set_data(sfs)
## set up effective population size
onepop_model.add_size_param("ne_a",lower=ne_bot,upper=ne_top*2) 
## set up populations and phylogeny
onepop_model.add_leaf("BOTH",N="ne_a")
## randomize parameters and check them
onepop_model.set_params(randomize=True)

##### PURE ISOLATION MODEL #####
print("\nPure Isolation model (base model)")
if (mu == 0):
	pure_isolation_model = momi.DemographicModel(N_e=ne_anc, gen_time=gentime) 
else:
	pure_isolation_model = momi.DemographicModel(N_e=ne_anc, muts_per_gen=mu, gen_time=gentime) 
	
pure_isolation_model.set_data(sfs)
## set up divergence times
pure_isolation_model.add_time_param("tdiv_sc",lower=divrange_bot,upper=divrange_top)
## set up effective population size
pure_isolation_model.add_size_param("ne_s",lower=ne_bot,upper=ne_top) ## this is from Brian's paper on cardinals
pure_isolation_model.add_size_param("ne_c",lower=ne_bot,upper=ne_top) 
## set up populations and phylogeny
pure_isolation_model.add_leaf("SON",N="ne_s")
pure_isolation_model.add_leaf("CHI",N="ne_c")
pure_isolation_model.move_lineages("SON", "CHI", t="tdiv_sc")
## randomize parameters and check them
pure_isolation_model.set_params(randomize=True)

##### CHI TO SON MODEL #####
print("\nC2S model (iso as base)")
uni_chi_to_son_model = pure_isolation_model.copy()
uni_chi_to_son_model.add_pulse_param("mig_c2s",lower=mig_bot,upper=mig_top)
uni_chi_to_son_model.add_time_param("tmig_c2s",lower=migdate_bot,upper=migdate_top,upper_constraints=["tdiv_sc"])
uni_chi_to_son_model.move_lineages("SON","CHI",t="tmig_c2s",p="mig_c2s")
uni_chi_to_son_model.set_params(randomize=True)

##### SON TO CHI MODEL #####

## start unidirectional models
print("\nS2C model (iso as base)")
uni_son_to_chi_model = pure_isolation_model.copy()
uni_son_to_chi_model.add_pulse_param("mig_s2c",lower=mig_bot,upper=mig_top)
uni_son_to_chi_model.add_time_param("tmig_s2c",lower=migdate_bot,upper=migdate_top,upper_constraints=["tdiv_sc"])
uni_son_to_chi_model.move_lineages("CHI","SON",t="tmig_s2c",p="mig_s2c")
uni_son_to_chi_model.set_params(randomize=True)

##### ASYMMETRIC MIGRATION #####
print("\nAsymmetric model (iso as base)")
asym_model = pure_isolation_model.copy() ## copy isol
asym_model.add_pulse_param("mig_s2c",lower=mig_bot,upper=mig_top)
asym_model.add_pulse_param("mig_c2s",lower=mig_bot,upper=mig_top)
asym_model.add_time_param("tmig_asym",lower=migdate_bot,upper=migdate_top,upper_constraints=["tdiv_sc"])
asym_model.move_lineages("CHI","SON",t="tmig_asym",p="mig_s2c")
asym_model.move_lineages("SON","CHI",t="tmig_asym",p="mig_c2s")
asym_model.set_params(randomize=True)

##### SYMMETRIC MIGRATION #####
print("\nSymmetric model (iso as base)")
symm_model = pure_isolation_model.copy()
symm_model.add_pulse_param("mig_symm",lower=mig_bot,upper=mig_top)
symm_model.add_time_param("tmig_symm",lower=migdate_bot,upper=migdate_top,upper_constraints=["tdiv_sc"])
symm_model.move_lineages("CHI","SON",t="tmig_symm",p="mig_symm")
symm_model.move_lineages("SON","CHI",t="tmig_symm",p="mig_symm")
symm_model.set_params(randomize=True)

##### SEC CONTACT #####
print("\nSec contact model (iso as base)")
## constrain time of migration to be after LGM, ~25000 years or sooner
sec_contact_model = pure_isolation_model.copy()
sec_contact_model.add_pulse_param("mig_secc",lower=mig_bot,upper=mig_top)
sec_contact_model.add_time_param("tmig_secc",lower=seccon_bot,upper=seccon_top,upper_constraints=["tdiv_sc"])
sec_contact_model.move_lineages("CHI","SON",t="tmig_secc",p="mig_secc")
sec_contact_model.move_lineages("SON","CHI",t="tmig_secc",p="mig_secc")
sec_contact_model.set_params(randomize=True)

##### OPTIMIZE MODELS #####
models = [onepop_model,pure_isolation_model,asym_model,uni_chi_to_son_model,uni_son_to_chi_model,sec_contact_model,symm_model]
model_names = ["1POP","ISOL","ASYM","UC2S","US2C","SECC","SYMM"]
AICs = []
count = 0
for model in models:
    now = datetime.datetime.now()
    name = str(model_names[count])
    print("Stochastic optimizing "+name+" model: "+str(now))
    model.stochastic_optimize(num_iters=1000, snps_per_minibatch=1000, save_to_checkpoint="momi_checkpoint_isol.txt", svrg_epoch=-1)
    now = datetime.datetime.now()
    print("Finished stochastic optimizing "+name+": "+str(now))
    print("Fully optimizing "+name+" model: "+str(now))
    model.optimize(method="L-BFGS-B")
    now = datetime.datetime.now()
    print("Finished fully optimizing "+name+": "+str(now))
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
print(model_names)
print("Delta AIC per model:\n", delta_aic)
print(model_names)
print("AIC weight per model:\n", np.exp(-0.5 * delta_aic))
print(model_names)

## gotta add in those bootstraps

## Make bootstrappies

#n_bootstraps = 10
# make copies of the original models to avoid changing them
#pure_isolation_model_copy = pure_isolation_model.copy()
#add_pulse_copy = add_pulse_model.copy()

# bootstrap_results = []
# for i in range(n_bootstraps):
#     print(f"Fitting {i+1} out of {n_bootstraps}")
# 
#     # resample the data
#     resampled_sfs = sfs.resample()
#     # tell models to use the new dataset
#     pure_isolation_model_copy.set_data(resampled_sfs)
#     #add_pulse_copy.set_data(resampled_sfs)
# 
#     # choose new random parameters for submodel, optimize
#     pure_isolation_model_copy.set_params(randomize=True)
#     pure_isolation_model_copy.optimize()
#     # initialize parameters from submodel, randomizing the new parameters
#     #add_pulse_copy.set_params(pulse_copy.get_params(),
#                               #randomize=True)
#     #add_pulse_copy.optimize()
# 
#     bootstrap_results.append(pure_isolation_model_copy.get_params())
# 
# ## Save bootstrappies
# with open("test-bootstrap-results-isol.txt", 'w') as outfile:
#     for bootstrap in bootstrap_results:
#         outfile.write("{}\n".format(bootstrap))

# Visualize bootstrappies
#make canvas, but delay plotting the demography (draw=False)
#yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5, 2.5e6, 7.5e6, 1e7]
#
#fig = momi.DemographyPlot(
#    pure_isolation_model, ["SON", "CHI"],
#    #linthreshy=1e5, 
#    figsize=(6,8),
#    major_yticks=yticks,
#    draw=False)
#plot bootstraps onto the canvas in transparency
#for params in bootstrap_results:
#    fig.add_bootstrap(
#        params,
#        #alpha=0: totally transparent. alpha=1: totally opaque
#        alpha=1/10)
#now draw the inferred demography on top of the bootstraps
#fig.draw()
#fig.draw_N_legend(loc="upper left")

