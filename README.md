# Summary Measures of Performance for Predictive Algorithms

This repository contains the code, data and an R markdown document for recreating the figures and exploring the examples in the paper.

## Getting Started

Download the repository; you should find three directories:

  1. `Code` -- contains the tutorial/walkthrough as an `html` file, the corresponding R markdown notebook and supplementary R code to use the notebook
  2. `Data` -- sample simulation data (including output from Stan / MCMC fitting of the model)
  3. `Figure` -- the individual figures produced by the example code

The main document containing the details is : `Code/revised_predictive_decisions.html` (open in a web browser) -- this is standalone and can be read without installing / using R or Stan

To experiment / run the examples yourself, open the corresponding markdown file `revised_predictive_decisions.Rmd` in RStudio

### Prerequisites

To read the extended example / walkthrough, simply open the html file (described above) in a browser.

To run the example markdown notebook, you will need:
  * the [R](https://www.r-project.org/) environment
  * [RStudio](https://rstudio.com/)

and to install the following R packages:

  * knitr
  * kableExtra
  * ggplot2
  * png
  * grid
  * gridExtra
  * Hmisc
  * dplyr
  * reshape2
  * caret
  * pROC
 
To run the Stan code, you will also need to install and configure [Stan](https://mc-stan.org/)

Note that output from 
 
