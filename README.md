##############################################################################################################################
# Development of prediction models for outcome of surgery for CMC1 OA
##############################################################################################################################
#
# This is the code used in: 
#
# "Using a machine learning approach to predict clinical improvement following surgery for thumb carpometacarpal osteoarthritis"
#
# Authors: NL Loos; L Hoogendam; JS Souer; HP Slijper; ER Andrinopoulou; MW Coppieters; the Hand-Wrist Study Group‡, RW Selles
#
##############################################################################################################################
# This code was written by N.L. Loos & L. Hoogendam for the Hand Wrist Study Group, Erasmus Medical Center, Rotterdam
##############################################################################################################################
#
# This repository contains two scripts: one for the development and validation of the prediction model for pain and one for the
# model for function
#
# These scripts is not an Extraction, Transformation and Load script. Some processing will probably be required for use on different
# data sets. This script is specific to our data. For the scripts to be used on different data set, some adaptation will be needed.
#
# The script prepares data for analysis. It splits the data into a training, validation and testdata set. Performs resamplin when 
# necessary (for the pain data set) and fits several models to the data: Generalized Linear Model (GLM), Gradient Boosting 
# Machine (GBM) and Random Forest (RF). The models are trained and validated on separate data sets. Hopefully this script might provide 
# those who are interested in developing prediction models with some tips on how to perform this on their own data. 
# 
# For questions or additional contact the corresponding author of the paper.
#
##############################################################################################################################
