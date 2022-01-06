# argwas


Install tskit:
conda install -c conda-forge tskit

Install limix_lmm:
conda install -c conda-forge limix

Download and unpack gcta:
https://cnsgenomics.com/software/gcta/#Download, provide correct path to executable in run_gcta_HE.sh

Build environment with:
conda env create -f argwas_environment.yml

Install statsmodel package from within environment:
conda install -c conda-forge statsmodels

Install tskit package from within environment:
conda install -c conda-forge tskit

Install plinkFile library in R

Install log indenter (I think there is no conda package):
pip install python_log_indenter
