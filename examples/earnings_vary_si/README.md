This example is based on the earnings_vary_si.stan model from the ARM/Ch.13 directory of the Stan examples.

First, run the notebook Sampling_earnings_vary_si.ipynb to run the sampler, evaluate the sensitivity,
and re-run the model with a range of prior parameters to validate the sensitivity results.  Because
the sampling is a little time-consuming, this will also save an Rdata file so the results
can be analyzed without re-running the sampler.

In this case, the sensitivity analysis uncovered a real problem with the Stan model.  For detailed
analysis of what's going on, run the notebook Analysis_earnings_vary_si.ipynb, which takes the
Rdata file from the first notebook as input.  At the end, we run a sampler on a fixed model, which
gives a considerably different posterior than the original.