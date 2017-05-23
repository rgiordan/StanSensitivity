# StanSensitivity

This repo contains a (very early) draft of some tools that allow users
to automatically generate local sensitivity measures to hyperparameters in
[Stan](http://mc-stan.org/), specifically in [RStan](http://mc-stan.org/interfaces/rstan.html).

A Stan model's data block typically consists of both data and hyperparameters.
Local sensitivity of a posterior mean to the hyperparameters can be
computed as the posterior covariance of a log derivative with with parameter[1].
  These tools allow you to
easily specify which hyperparameters you're intersted in and automatically
calculate the necessary local sensitivity.

Here are the current steps to use it.  Some examples can be found in
```example_models```, which are based on their namesakes in the
[Stan examples](https://github.com/stan-dev/example-models).
1. Start with an existing Stan model and dataset containing hyperparameters (e.g. ```example_models/negative_binomial/negative_binomial.stan``` and ```example_models/negative_binomial.data.R```).
2. Break your model up into chunks using ```python/split_model.py```.  You'll need to
specify the "base name", which is everything in your stan model file except the
extension.  This base name will be used to generate derived scripts for calculating
sensitivity.
```
python/split_model.py --model_name=example_models/negative_binomial/negative_binomial
```
3. This will make a bunch of ```*.stanblock``` files.  You'll need to manually
edit two of them to specify your hyperparameters.  Look in the file ```BASENAME_data_block.stanblock```.  You'll see the contents of your original
datablock, including both data and hyperparameters.  Move (don't copy) the
hyperparameters to the file ```BASENAME_hyperparameters_block.stanblock```.
For example, I moved the parameters ```cauchy_*``` and ```weights```
from the data block to the hyperparameters block in the negative binomial model.
4. Ren ```generate_models.py``` using the same syntax as above:
```
python/generate_models.py --model_name=example_models/negative_binomial/negative_binomial
```
5. Run the R script ```R/run_examples.R```, setting the ```model_name``` variable
to the value of the Python ```--model_name``` flag, e.g.
```
model_name <- "example_models/negative_binomial//negative_binomial"
```
6. Hopefully, if everything works, you will have your local sensitivity measures
in the matrix ```sens_mat```, which you can graph or interpret as you see fit.

This is currently just a sketch, so no guarantees that
it's working yet.  Please contact me with questions or comments.

[1]: [Fast Measurements of Robustness to Changing Priors in Variational Bayes](https://arxiv.org/abs/1611.07469)
