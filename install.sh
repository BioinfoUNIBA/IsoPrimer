#!/bin/bash

# >>> Conda setup >>>
eval "$(conda shell.bash hook)"
# <<< Conda setup <<<

conda create -y --name isoprimer

conda activate isoprimer && \
 conda install -y conda-forge::r-base=4.3.3 \
 bioconda::emboss \
 bioconda::primer3 \
 bioconda::kallisto=0.46.2 \
 conda-forge::r-dplyr \
 conda-forge::r-biocmanager \
 bioconda::bioconductor-biostrings \
 conda-forge::r-doparallel \
 conda-forge::r-openxlsx \
 conda-forge::r-stringr \
 bioconda::bioconductor-msa && \
 echo 'IsoPrimer dependencies installed successfully'
