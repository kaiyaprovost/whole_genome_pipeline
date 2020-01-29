import momi		## momi2 analysis
import logging		## create log file
import numpy as np
import datetime

print("-----\n-----\n-----")

print("start logging\n")
logging.basicConfig(level=logging.INFO,
					filename="momi_log.txt")

print("load sfs\n")
sfspath = "/Users/kprovost/Downloads/outfiles_allgroup/cardcard_files2/cardcard16_sfs_filtered_changelength_monomorphic.txt"
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
print("DIV TIME RANGE: 100,000 to 1,000,000")
print("NE RANGE: 1,000 to 2,000,000")
print("MIGRATION RANGE: 0 to 0.1")
print("MIGRATION DATE RANGE: 0 to 1,000,000 [25,000 for SECC model]\n\n")

print("begin setting up models\n")
##### PURE ISOLATION MODEL #####
print("Pure Isolation model (base model)")
pure_isolation_model = momi.DemographicModel(N_e=1000000,muts_per_gen=2.21e-9,gen_time=1) ## why tho -- can you give it something?
pure_isolation_model.set_data(sfs)
## set up divergence times
pure_isolation_model.add_time_param("tdiv_sc",lower=100000,upper=1000000)
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
## check that the model is what you want 
yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5,1e6,3e6]
fig = momi.DemographyPlot(
	pure_isolation_model, 
	["Chi","Son"],
	figsize=(6,8),
	major_yticks=yticks,
	#linthreshy=2.5e6
)

## set up the rest of the models 

##### UNIDIRECTIONAL CHI TO SON MIGRATION #####
print("Uni Chi2Son model")
uni_chi_to_son_model = pure_isolation_model.copy() ## copies the Pure Isolation model 
## add additional parameters 
uni_chi_to_son_model.add_pulse_param("mig_c2s",lower=0,upper=0.1)
uni_chi_to_son_model.add_time_param("tmig_c2s",lower=0,upper=1000000,upper_constraints=["tdiv_sc"])
uni_chi_to_son_model.move_lineages("Son","Chi",t="tmig_c2s",p="mig_c2s")
## randomize parameters and check them
uni_chi_to_son_model.set_params(randomize=True)
print(uni_chi_to_son_model.get_params())
## check that the model is what you want 
yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5,1e6]
fig = momi.DemographyPlot(
	uni_chi_to_son_model, 
	["Chi","Son"],
	#["Out","Chi","Son"],
	figsize=(6,8),
	major_yticks=yticks,
	#linthreshy=2.5e6
)

##### UNIDIRECTIONAL SON TO CHI MIGRATION #####
print("Uni Son2Chi model")
uni_son_to_chi_model = pure_isolation_model.copy() ## copies the Pure Isolation model
## add additional parameters
uni_son_to_chi_model.add_pulse_param("mig_s2c",lower=0,upper=0.1)
uni_son_to_chi_model.add_time_param("tmig_s2c",lower=0,upper=1000000,upper_constraints=["tdiv_sc"])
uni_son_to_chi_model.move_lineages("Chi","Son",t="tmig_s2c",p="mig_s2c")
## randomize parameters and check them
uni_son_to_chi_model.set_params(randomize=True)
print(uni_son_to_chi_model.get_params())
## check that the model is what you want 
yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5,1e6]
fig = momi.DemographyPlot(
	uni_son_to_chi_model, 
	["Chi","Son"],
	#["Out","Chi","Son"],
	figsize=(6,8),
	major_yticks=yticks,
	#linthreshy=2.5e6
)

##### ASYMMETRIC MIGRATION #####
print("Asymmetric model (c2s as base)")
asym_model = uni_chi_to_son_model.copy() ## copy c2s model 
## add s2c parameters
asym_model.add_pulse_param("mig_s2c",lower=0,upper=0.1)
asym_model.add_time_param("tmig_s2c",lower=0,upper=1000000,upper_constraints=["tdiv_sc"])
asym_model.move_lineages("Chi","Son",t="tmig_s2c",p="mig_s2c")
## randomize and check parameters
asym_model.set_params(randomize=True)
print(asym_model.get_params())
## check demography figure
yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5,1e6]
fig = momi.DemographyPlot(
	asym_model, 
	["Chi","Son"],
	#["Out","Chi","Son"],
	figsize=(6,8),
	major_yticks=yticks,
	#linthreshy=2.5e6
)

##### AYMMETRIC MIGRATION #####
print("Symmetric model")
symm_model = pure_isolation_model.copy() ## copy Pure Isolation model
## add parameters
symm_model.add_pulse_param("mig_symm",lower=0,upper=0.1)
symm_model.add_time_param("tmig_symm",lower=0,upper=1000000,upper_constraints=["tdiv_sc"])
symm_model.move_lineages("Chi","Son",t="tmig_symm",p="mig_symm")
symm_model.move_lineages("Son","Chi",t="tmig_symm",p="mig_symm")
## randomize and get params
symm_model.set_params(randomize=True)
print(symm_model.get_params())
## print demgography
yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5,1e6]
fig = momi.DemographyPlot(
	symm_model, 
	["Chi","Son"],
	#["Out","Chi","Son"],
	figsize=(6,8),
	major_yticks=yticks,
	#linthreshy=2.5e6
)

##### SECONDARY CONTACT #####
print("Secondary contact")
## constrain time of migration to be after LGM, ~25000 years or sooner
sec_contact_model = pure_isolation_model.copy() ## copy the Pure Isolation model
## add parameters
sec_contact_model.add_pulse_param("mig_secc",lower=0,upper=0.1)
sec_contact_model.add_time_param("tmig_secc",lower=0,upper=25000,upper_constraints=["tdiv_sc"])
sec_contact_model.move_lineages("Chi","Son",t="tmig_secc",p="mig_secc")
sec_contact_model.move_lineages("Son","Chi",t="tmig_secc",p="mig_secc")
## randomize and get parameters
sec_contact_model.set_params(randomize=True)
print(sec_contact_model.get_params())
## check demography
yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5,1e6]
fig = momi.DemographyPlot(
	sec_contact_model, 
	["Chi","Son"],
	#["Out","Chi","Son"],
	figsize=(6,8),
	major_yticks=yticks,
	#linthreshy=2.5e6
)

## optimize each model once
print("#####")

models = [pure_isolation_model,uni_chi_to_son_model,uni_son_to_chi_model,asym_model,symm_model,sec_contact_model]
model_names = ["ISOL","UC2S","US2C","ASYM","SYMM","SECC"]

AICs = []
count = 0
for model in models:
	now = datetime.datetime.now()
	name = str(model_names[count])
	print("Optimizing "+name+" model: "+str(now))
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
	fig = momi.DemographyPlot(
		model, 
		["Chi","Son"],
		#["Out","Chi","Son"],
		figsize=(6,8),
		major_yticks=yticks,
		#linthreshy=2.5e6
	)
	count += 1
	print("-----")

minv = np.min(AICs)
delta_aic = np.array(AICs) - minv
print("Delta AIC per model: ", delta_aic)
print("AIC weight per model: ", np.exp(-0.5 * delta_aic))
	
## TODO: add searching multiple times
## TODO: add bootstrapping

