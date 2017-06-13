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
num_warmup_samples <- 200
num_samples <- 200
sampling_result <- sampling(model, data=stan_data, chains=1, iter=(num_samples + num_warmup_samples))
print(summary(sampling_result))


##################################
# Get the sensitivity model and sensitivity.

#model_sens <- stan_model(paste(model_name, "_sensitivity.stan", sep=""))

# model_sens_params <- stan(paste(model_name, "_sensitivity_parameters.stan", sep=""),
#                           data=stan_data, algorithm="Fixed_param",
#                           iter=1, chains=1)
# sens_par_list <- get_inits(model_sens_params)[[1]]
# model_par_list <- get_inits(sampling_result)[[1]]
#
# # Get the sensitivity parameters in list form.
# for (par in names(model_par_list)) {
#   cat("Copying parameter '", par, "' from the sampler\n", sep="")
#   sens_par_list[[par]] <- model_par_list[[par]]
# }
# for (par in names(sens_par_list)) {
#   if (par %in% names(stan_data)) {
#     cat("Copying hyperparameter '", par, "' from the data.\n", sep="")
#     sens_par_list[[par]] <- stan_data[[par]]
#   }
# }
#
# model_sens_fit <- stan(paste(model_name, "_sensitivity.stan", sep=""),
#                        data=stan_data, algorithm="Fixed_param",
#                        iter=1, chains=1, init=list(sens_par_list))

draws_mat <- extract(sampling_result, permute=FALSE)[,1,]
stan_sensitivity_list <- GetStanSensitivityModel(sampling_result, model_name, stan_data)

# # These names help sort through the vectors of sensitivity.
# param_names <-
#   sampling_result@.MISC$stan_fit_instance$unconstrained_param_names(FALSE, FALSE)
# sens_param_names <-
#   model_sens_fit@.MISC$stan_fit_instance$unconstrained_param_names(FALSE, FALSE)
#
# stan_sensitivity_list <-
#   list(model_sens_fit=model_sens_fit,
#      param_names=param_names,
#      sens_param_names=sens_param_names,
#      sens_par_list=sens_par_list)
#
# model_sens_fit <- stan_sensitivity_list$model_sens_fit
# param_names <- stan_sensitivity_list$param_names
# sens_param_names <- stan_sensitivity_list$sens_param_names
# sens_par_list_init <- stan_sensitivity_list$sens_par_list
# sens_par_list <- sens_par_list_init
#
# # Get the model gradients with respect to the hyperparameters (and parameters).
# num_samples <- nrow(draws_mat)
# grad_mat <- matrix(NA, num_samples, length(sens_param_names))
# lp_vec <- rep(NA, num_samples)
# cat("Evaluating log gradients at the MCMC draws.\n")
# prog_bar <- txtProgressBar(min=1, max=num_samples, style=3)
#
# for (n in 1:num_samples) {
#   setTxtProgressBar(prog_bar, value=n)
#   par_list <- get_inits(sampling_result, iter=n + num_warmup_samples)[[1]]
#   for (par in ls(par_list)) {
#     # Note that get_inits is currently broken:
#     # https://github.com/stan-dev/rstan/issues/417
#     # ...but this has seemed to fix it (so far):
#     if (length(dim(sens_par_list[[par]])) >= 2) {
#       sens_par_list[[par]] <- array(unlist(par_list[[par]]), dim(par_list[[par]]))
#     } else {
#       sens_par_list[[par]] <- as.numeric(par_list[[par]])
#     }
#   }
#   pars_free <- unconstrain_pars(model_sens_fit, sens_par_list)
#   glp <- grad_log_prob(model_sens_fit, pars_free)
#   grad_mat[n, ] <- glp
#   lp_vec[n] <- attr(glp, "log_prob")
# }
#
#
# close(prog_bar)

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



