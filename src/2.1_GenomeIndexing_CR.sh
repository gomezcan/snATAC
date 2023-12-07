#!/bin/bash

###################################
#######       Modules     #########
###################################

ml Bioinformatics cellrangerATAC/2.1.0

###################################
###    Command lines to Run     ###
###################################

# Requare config file
cellranger-atac mkref --config=cellrangerATAC_B63_v5.config
###################################
