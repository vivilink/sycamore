argwas
=======

Creating conda environment
-------------------------

Build environment

    conda env create -f argwas_environment.yml

Install tskit

    conda install -c conda-forge tskit

Install limix_lmm

    conda install -c conda-forge limix

Download and unpack gcta

    https://cnsgenomics.com/software/gcta/#Download, provide correct path with parameter --GCTA

Install statsmodel

    conda install -c conda-forge statsmodels

Install log indenter (I think there is no conda package)

    pip install python_log_indenter

Install tsinfer

    conda install -c conda-forge tsinfer 

Install matplotlib

    conda install -c conda-forge matplotlib
 
Install tsdate

    python3 -m pip install tsdate --user


Install plinkFile library in R if you want to write covariance matrices in .grm.bin format for GCTA

    install.packages("plinkFile")

