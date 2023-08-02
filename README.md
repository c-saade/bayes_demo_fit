# bayes_demo_fit
📈 Learn to fit demographic models to time series with Rstan.

In this project, I show how to fit a two-species demographic model onto real experimental data from an ecology [research paper](https://doi.org/10.1098/rspb.2022.0543) for a post on my personal website.

Contents:
  - __Bayes_demo_fit.Rmd__ and __Bayes_demo_fit.html__: a Rmarkdown document (and is html output) explaining the procedure step by step.
  - __Draft_model.R__: a standalone R script to run the model fit without going through the Rmd file.
  - __prior_pred_viz.R__: an R script to vizualize prior predictions.
  - __Data__: a folder with both the data and data description.
  - __out__: a folder with the ouput from Draft_model.R (a stanfit object and figures).
