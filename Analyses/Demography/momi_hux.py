
# coding: utf-8

# In[1]:


import momi        ## momi2 analysis
import logging        ## create log file
import numpy as np
import datetime
import matplotlib as plt


# In[2]:


print("-----\n-----\n-----")

print("start logging\n")
logging.basicConfig(level=logging.INFO,
                    filename="momi_log_test.txt")


# In[3]:


get_ipython().run_cell_magic('sh', '', '\necho "start"\ncd /Users/kprovost/Dropbox\\ \\(AMNH\\)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/\n\nls *vcf')


# In[4]:


get_ipython().run_cell_magic('bash', '', '\ncd /Users/kprovost/Dropbox\\ \\(AMNH\\)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/\n\nhead -n 2 Amphispiza-bilineata-called.geno.vcf\n## bgzip performs a blockwise compression\n## The -c flag directs bgzip to leave the original vcf file \n##   untouched and create a new file for the vcf.gz\n')


# In[5]:


get_ipython().run_cell_magic('bash', '', 'cd /Users/kprovost/Dropbox\\ \\(AMNH\\)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/\n\n\nbgzip -c Amphispiza-bilineata-called.geno.vcf > Amphispiza-bilineata-called.geno.vcf.gz\n\n## tabix indexes the file for searching\ntabix Amphispiza-bilineata-called.geno.vcf.gz')


# In[6]:


get_ipython().run_cell_magic('bash', '', '\ncd /Users/kprovost/Dropbox\\ \\(AMNH\\)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/\n\nvcf2bed --do-not-split \\\n--do-not-sort < Amphispiza-bilineata-called.geno.vcf \\\n> Amphispiza-bilineata-called.geno.bed\n\n## Print the first 1 lines of this file\nhead -n 1 Amphispiza-bilineata-called.geno.bed\n\n')


# In[ ]:


get_ipython().run_cell_magic('bash', '', '\ncd /Users/kprovost/Dropbox\\ \\(AMNH\\)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/\n\npython -m momi.read_vcf \\\n--no_aa \\\n--verbose \\\nAmphispiza-bilineata-called.geno.vcf.gz \\\nAmphispiza-bilineata-called.popassignments.csv \\\nAmphispiza-bilineata-called.allelecounts.txt \\\n--bed Amphispiza-bilineata-called.geno.bed\n\n#for i in *vcf.gz; do echo $i; j=${i%.geno.vcf.gz}; echo $j; python -m momi.read_vcf --no_aa --verbose $i $j.popassignments.csv $j.allelecounts.txt --bed $j.geno.bed; done')


# In[8]:


get_ipython().run_cell_magic('bash', '', '\ncd /Users/kprovost/Dropbox\\ \\(AMNH\\)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/\n\npython -m momi.extract_sfs Amphispiza-bilineata-called.geno.sfs.txt 20 Amphispiza-bilineata-called.geno.vcf')


# In[9]:


print("load sfs\n")
## sfspath = "/home/kprovost/nas1/momi2/cardcard16_sfs_filtered_changelength_monomorphic.txt"
sfspath = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/Toxostoma-crissale-called_sfs.txt"
## this is a two-population sfs with monomorphic sites included in "length"
sfs = momi.Sfs.load(sfspath)
#print("Avg pairwise heterozygosity", sfs.avg_pairwise_hets[:5])
#print("populations", sfs.populations)
#print("percent missing data per population", sfs.p_missing)


# In[ ]:


print("\nPRIORS")
print("MUT RATE: 2.21e-9")
print("ANCESTRAL NE: 300,000")
print("GEN TIME: 1")
print("DIV TIME RANGE: 5,000 to 2,500,000")
print("NE RANGE: 1,000 to 20,000,000")
print("MIGRATION RANGE: 0 to 0.5")
print("MIGRATION DATE RANGE: 0 to 1,000,000 [25,000 for SECC model]\n\n")


print("begin setting up models\n")
##### PURE ISOLATION MODEL #####
print("\nPure Isolation model (base model)")
pure_isolation_model = momi.DemographicModel(N_e=300000,muts_per_gen=2.21e-9,gen_time=1) ## why tho -- can you give it something?
pure_isolation_model.set_data(sfs)
## set up divergence times
pure_isolation_model.add_time_param("tdiv_sc",lower=5000,upper=2500000)
## set up effective population size
pure_isolation_model.add_size_param("ne_s",lower=1000,upper=20000000) ## this is from Brian's paper on cardinals
pure_isolation_model.add_size_param("ne_c",lower=1000,upper=20000000) 
## set up populations and phylogeny
pure_isolation_model.add_leaf("SON",N="ne_s")
pure_isolation_model.add_leaf("CHI",N="ne_c")
pure_isolation_model.move_lineages("SON", "CHI", t="tdiv_sc")
## randomize parameters and check them

pure_isolation_model.set_params(randomize=True)
print(pure_isolation_model.get_params())


# In[ ]:


## check that the model is what you want 
yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5,1e6,3e6]
fig = momi.DemographyPlot(
	pure_isolation_model, 
	["CHI","SON"],
	figsize=(6,8),
	major_yticks=yticks,
	#linthreshy=2.5e6
)


# In[ ]:


## STOCH OPT WITH SNPS
## SNPS IS FASTER?

pure_isolation_model.stochastic_optimize(num_iters=10, snps_per_minibatch=1000, save_to_checkpoint="momi_checkpoint_isol.txt", svrg_epoch=-1)
print(pure_isolation_model.get_params())



## TRY STOCHASTIC OPTIMIZE
## STILL FAST BUT NOT SO FAST

pure_isolation_model.stochastic_optimize(num_iters=10, n_minibatches=5, save_to_checkpoint="momi_checkpoint_isol.txt", svrg_epoch=-1)
print(pure_isolation_model.get_params())


## TRY THE OPTIMIZE WITH THE OTHER MODEL FROM HERE


pure_isolation_model.optimize(method="L-BFGS-B")
print(pure_isolation_model.get_params())





# In[ ]:


## start unidirectional models
## chi to son migration 

uni_chi_to_son_model = pure_isolation_model.copy()
uni_chi_to_son_model.set_params(randomize=True)

uni_chi_to_son_model.add_pulse_param("mig_c2s",lower=0,upper=0.1)
uni_chi_to_son_model.add_time_param("tmig_c2s",lower=0,upper=3000000,upper_constraints=["tdiv_sc"])
uni_chi_to_son_model.move_lineages("SON","CHI",t="tmig_c2s",p="mig_c2s")
uni_chi_to_son_model.get_params()

yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5,1e6]
fig = momi.DemographyPlot(
    uni_chi_to_son_model, 
    ["CHI","SON"],
    #["Out","CHI","SON"],
    figsize=(6,8),
    major_yticks=yticks,
    #linthreshy=2.5e6
)


# In[ ]:


## optimize model, get params

import datetime
now = datetime.datetime.now()
print("")
print("Current date and time using str method of datetime object:")
print(str(now))

uni_chi_to_son_model.optimize()
uni_chi_to_son_model.get_params()

yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5,1e6]
fig = momi.DemographyPlot(
    uni_chi_to_son_model, 
    ["CHI","SON"],
    #["Out","CHI","SON"],
    figsize=(6,8),
    major_yticks=yticks,
    #linthreshy=2.5e6
)


# In[ ]:


## son to chi migration

## start unidirectional models

uni_son_to_chi_model = pure_isolation_model.copy()
uni_son_to_chi_model.set_params(randomize=True)

uni_son_to_chi_model.add_pulse_param("mig_s2c",lower=0,upper=0.1)
uni_son_to_chi_model.add_time_param("tmig_s2c",lower=0,upper=3000000,upper_constraints=["tdiv_sc"])
uni_son_to_chi_model.move_lineages("CHI","SON",t="tmig_s2c",p="mig_s2c")
uni_son_to_chi_model.get_params()

yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5,1e6]
fig = momi.DemographyPlot(
    uni_son_to_chi_model, 
    ["CHI","SON"],
    #["Out","CHI","SON"],
    figsize=(6,8),
    major_yticks=yticks,
    #linthreshy=2.5e6
)


# In[ ]:


##### ASYMMETRIC MIGRATION #####
print("\nAsymmetric model (iso as base)")
asym_model = pure_isolation_model.copy() ## copy isol
asym_model.add_pulse_param("mig_s2c",lower=0,upper=0.1)
asym_model.add_pulse_param("mig_c2s",lower=0,upper=0.1)
asym_model.add_time_param("tmig_asym",lower=0,upper=1000000,upper_constraints=["tdiv_sc"])
asym_model.move_lineages("CHI","SON",t="tmig_asym",p="mig_s2c")
asym_model.move_lineages("SON","CHI",t="tmig_asym",p="mig_c2s")
## randomize and check parameters
#asym_model.set_params(randomize=True)
#print(asym_model.get_params())


# In[ ]:


## symmetric migration

symm_model = pure_isolation_model.copy()
symm_model.set_params(randomize=True)

symm_model.add_pulse_param("mig_symm",lower=0,upper=0.1)
symm_model.add_time_param("tmig_symm",lower=0,upper=3000000,upper_constraints=["tdiv_sc"])
symm_model.move_lineages("CHI","SON",t="tmig_symm",p="mig_symm")
symm_model.move_lineages("SON","CHI",t="tmig_symm",p="mig_symm")
symm_model.get_params()

yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5,1e6]
fig = momi.DemographyPlot(
    symm_model, 
    ["CHI","SON"],
    #["Out","CHI","SON"],
    figsize=(6,8),
    major_yticks=yticks,
    #linthreshy=2.5e6
)


# In[ ]:


## secondary contact model
## constrain time of migration to be after LGM, ~25000 years or sooner

sec_contact_model = pure_isolation_model.copy()
sec_contact_model.set_params(randomize=True)
sec_contact_model.add_pulse_param("mig_secc",lower=0,upper=0.1)
sec_contact_model.add_time_param("tmig_secc",lower=0,upper=25000,upper_constraints=["tdiv_sc"])
sec_contact_model.move_lineages("CHI","SON",t="tmig_secc",p="mig_secc")
sec_contact_model.move_lineages("SON","CHI",t="tmig_secc",p="mig_secc")
print(sec_contact_model.get_params())

yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5,1e6]
fig = momi.DemographyPlot(
    sec_contact_model, 
    ["CHI","SON"],
    #["Out","CHI","SON"],
    figsize=(6,8),
    major_yticks=yticks,
    #linthreshy=2.5e6
)

#fig.savefig("sec_contact_model.png")


# In[ ]:


models = [pure_isolation_model,asym_model,uni_chi_to_son_model,uni_son_to_chi_model,sec_contact_model,symm_model]
model_names = ["ISOL","ASYM","UC2S","US2C","SECC","SYMM"]
AICs = []
count = 0
for model in models:
	now = datetime.datetime.now()
	name = str(model_names[count])
	print("Stochastic optimizing "+name+" model: "+str(now))
	model.stochastic_optimize(num_iters=10, snps_per_minibatch=1000, save_to_checkpoint="momi_checkpoint_isol.txt", svrg_epoch=-1)
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
print("Delta AIC per model: ", delta_aic)
print("AIC weight per model: ", np.exp(-0.5 * delta_aic))


# In[ ]:


## gotta add in those bootstraps

## Make bootstrappies

n_bootstraps = 10
# make copies of the original models to avoid changing them
pure_isolation_model_copy = pure_isolation_model.copy()
#add_pulse_copy = add_pulse_model.copy()

bootstrap_results = []
for i in range(n_bootstraps):
	print(f"Fitting {i+1} out of {n_bootstraps}")

	# resample the data
	resampled_sfs = sfs.resample()
	# tell models to use the new dataset
	pure_isolation_model_copy.set_data(resampled_sfs)
	#add_pulse_copy.set_data(resampled_sfs)

	# choose new random parameters for submodel, optimize
	pure_isolation_model_copy.set_params(randomize=True)
	pure_isolation_model_copy.optimize()
	# initialize parameters from submodel, randomizing the new parameters
	#add_pulse_copy.set_params(pulse_copy.get_params(),
							  #randomize=True)
	#add_pulse_copy.optimize()

	bootstrap_results.append(pure_isolation_model_copy.get_params())

## Save bootstrappies
with open("test-bootstrap-results-isol.txt", 'w') as outfile:
	for bootstrap in bootstrap_results:
		outfile.write("{}\n".format(bootstrap))
		

#
#
#
#
		
		


# In[ ]:


# Visualize bootstrappies
#make canvas, but delay plotting the demography (draw=False)
yticks = [1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5, 2.5e6, 7.5e6, 1e7]

fig = momi.DemographyPlot(
	pure_isolation_model, ["SON", "CHI"],
	#linthreshy=1e5, 
    figsize=(6,8),
	major_yticks=yticks,
	draw=False)

#plot bootstraps onto the canvas in transparency
for params in bootstrap_results:
	fig.add_bootstrap(
		params,
		#alpha=0: totally transparent. alpha=1: totally opaque
		alpha=1/10)

#now draw the inferred demography on top of the bootstraps
fig.draw()
fig.draw_N_legend(loc="upper left")

