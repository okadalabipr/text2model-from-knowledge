#! /bin/bash

# This script is used to prepare the modified biomass model and the KEGG2Model script used in the supplementary case studies.

# git clone biomass v0.11.0
git clone https://github.com/biomass-dev/biomass -b v0.11.0 biomass_models/modified_biomass

# copy the modified files
cp biomass_modification/temporal_dynamics.py biomass_models/modified_biomass/biomass/dynamics/temporal_dynamics.py

# create symbolic link of figure2/KEGG2Model
ln -s ../figure2/KEGG2Model ./
