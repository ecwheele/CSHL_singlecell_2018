#!/usr/bin/env bash

conda create -n single-cell python=3.6;
source activate single-cell;

conda install -y -c bioconda -c anaconda -c conda-forge \
  pandas \
  numpy \
  scipy \
  scikit-learn \
  seaborn \
  matplotlib \
  networkx \
  fastcluster \
  jupyter;
	
pip install polo macosko2015;

pip install git+https://github.com/jacoblevine/phenograph.git;

conda install -y -c r r-essentials
conda install -y -c bioconda r-seurat
