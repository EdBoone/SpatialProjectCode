#!/bin/bash

# Load anaconda3/personal module and install if not installed already
module load anaconda3/personal
if [ -d $HOME/anaconda3 ]; then
echo -e "\nanaconda3/personal already installed\n\n";
else
  echo -e "\ninstalling anaconda3/personal\n\n";
anaconda-setup
fi


# Create new conda environment called "env_R405"
if [ -d $HOME/anaconda3/envs/env_R405 ]; then
echo "###############################################"
echo -e "\nenv_R405 conda environment is already present"
echo -e "\nIf you wish you re-install please remove the conda environment first with:"
echo -e "\tconda remove -n env_R405 --all -y"
echo -e "\n\n###############################################"
exit 1
else
  echo -e "\nCreating Conda environment for INLA : env_R405"
echo -e "\nR version: 4.0.5"
conda create -n env_R405 r-base=4.0.5 -c conda-forge -y
source activate env_R405
fi

# Install initial dependencies
echo -e "\nInstalling R packages via conda"
conda install -c conda-forge r-tidyverse -y
conda install r-deriv r-devtools r-doparallel r-fields r-foreach r-gridextra r-matrixmodels r-matrixstats r-mvtnorm r-numderiv r-pixmap r-rgl  r-sf r-sn r-splancs r-spdep r-sp r-polynom r-rgdal bioconductor-rgraphviz r-tidybayes  -c conda-forge -y
conda install -c conda-forge/label/gcc7 r-gmp -y
conda install -c conda-forge r-mpoly -y


# Install additional dependencies
echo -e "\nInstalling additional R packages"
R -e 'options(repos = c(CRAN = "https://cran.rstudio.com"));install.packages("geobr")' 
R -e 'options(repos = c(CRAN = "https://cran.rstudio.com"));install.packages("fastDummies")'
R -e 'options(repos = c(CRAN = "https://cran.rstudio.com"));install.packages("tsModel")' 

# INLA installed manually CRAN 60
# R -e 'options(repos = c(CRAN = "https://cran.rstudio.com"));install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)


# INSTRUCTIONS OF USE:
echo -e "\t\tINSTRUCTIONS FOR USE"
echo -e "\n################################################\n\n"
echo -e "Ensure you add the following lines in your jobscript:\n\n"
echo -e "\tmodule load anaconda3/personal"
echo -e "\tsource activate env_R405"
echo -e "\n\n"
echo -e "\n################################################\n\n"

