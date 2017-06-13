library(rstan)
library(rstansensitivity)

library(ggplot2)
library(dplyr)
library(reshape2)

rstan_options(auto_write=TRUE)

# Set this to be the appropriate location of the repo on your computer.
example_directory <- file.path(
  Sys.getenv("GIT_REPO_LOC"), "StanSensitivity/examples/example_models")

base_model_name <- file.path(example_directory, "schools/schools-1.stan")
python_script <- file.path(Sys.getenv("GIT_REPO_LOC"), "StanSensitivity/python/generate_models.py")
model_name <- GenerateSensitivityFromModel(base_model_name, python_script=python_script)


##################################
# Compile and run the base model.

model <- stan_model(paste(model_name, "_generated.stan", sep=""))

# Load the data and hyperparameters.
stan_data <- new.env()
source(paste(model_name, "data.R", sep="."), local=stan_data)
stan_data <- as.list(stan_data)

# For now, you must use chains=1 for now to avoid confusion around get_inits.
# The script currently assumes the same number of warm-up draws as final samples.
num_warmup_samples <- 1000
num_samples <- 1000
sampling_result <- sampling(model, data=stan_data, chains=1, iter=(num_samples + num_warmup_samples))
print(summary(sampling_result))


##################################
# Get the sensitivity model and sensitivity.

draws_mat <- extract(sampling_result, permute=FALSE)[,1,]
stan_sensitivity_list <- GetStanSensitivityModel(sampling_result, model_name, stan_data)

sens_result <-
  GetStanSensitivityFromModelFit(sampling_result, draws_mat, stan_sensitivity_list, num_warmup_samples)


##################################
# Inspect the results.

sens_mat <- sens_result$sens_mat
sens_mat <- sens_mat[, colnames(sens_mat) != "lp__"]
sens_mat_normalized <- sens_result$sens_mat_normalized
sens_mat_normalized <- sens_mat_normalized[, colnames(sens_mat_normalized) != "lp__"]

# Look at the sensitivity to other hyperparameters.
print(sens_mat_normalized)
rownames(sens_mat_normalized)

sens_mat_normalized_df <-
  data.frame(sens_mat_normalized) %>%
  mutate(hyperparameter=rownames(sens_mat_normalized)) %>%
  melt(id.var="hyperparameter") %>%
  mutate(base_variable=sub("\\..*", "", variable))

ggplot(sens_mat_normalized_df) +
  geom_bar(stat="identity", position="dodge",
           aes(x=hyperparameter, y=value, fill=variable)) +
  facet_grid(base_variable ~ .) +
  theme(legend.position="none")

unique(sens_mat_normalized_df$base_variable)
ggplot(filter(sens_mat_normalized_df, base_variable == "beta")) +
  geom_bar(stat="identity", position="dodge",
           aes(x=hyperparameter, y=value, fill=variable)) +
  facet_grid(base_variable ~ .) +
  theme(legend.position="none")



############################
# Demonstrate get_inits is screwed up.

examples_loc <- "/home/rgiordan/Documents/git_repos/example-models/bugs_examples/vol2"
stan_data <- new.env()
source(file.path(examples_loc, "schools/schools.data.R"), local=stan_data)
stan_data <- as.list(stan_data)

school_model <- stan_model(file.path(examples_loc, "schools/schools-1.stan"))
sampling_result <- sampling(school_model, data=stan_data, chains=1, iter=200)

par_list <- get_inits(sampling_result, iter=1)[[1]]
pars_free <- unconstrain_pars(sampling_result, par_list) # Fails
# Fails:
# Error in object@.MISC$stan_fit_instance$unconstrain_pars(pars) :
# variable beta missing

# It fails for every variable since they're all malformed.
# But the matrix parameter is especially strange:

class(par_list$alpha)
# [1] "matrix"
dim(par_list$alpha)
# [1] 38  3
class(par_list$alpha[1])
# [1] "list"
par_list$alpha[1, 1]
# [[1]]
# [1] 0.7710187

array(unlist(par_list$alpha), dim(par_list$alpha))



class(par_list$beta)
par_list$beta <- as.numeric(par_list$beta)
par_list$theta <- as.numeric(par_list$theta)
par_list$phi <- as.numeric(par_list$phi)
par_list$alpha <- as.numeric(par_list$alpha)

pars_free <- unconstrain_pars(sampling_result, par_list)

class(par_list$alpha[0])



