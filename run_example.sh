#!/bin/bash

#================================================================================================
#
# Signal Extraction
#-------------------
#  0: probe counting
#  1: Breit-Wigner convolved with Crystal Ball function
#  2: MC template convolved with Gaussian
#  3: Phil's Crystal Ball based "Voigtian" shape
#  4: Unbinned MC data convolved with Gaussian
#  5: BW convolved with Crystal Ball plus Voigtian
#
# Background Model
#------------------
#  0: no background
#  1: exponential model
#  2: erfc*exp model
#  3: double exponential model
#  4: linear*exp model
#  5: quadratic*exp model
#  6: erfc*exp model with background params from <resultsdir>
#________________________________________________________________________________________________
#
#================================================================================================
#
# effZFit usage
#---------------
### Need to do: export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`
# ./effZFit <conf> <sigpass> <bkgpass> <sigfail> <bkgfail> <infname> <outdir> <doPU> <charge> <temfile> <resultsdir>
#________________________________________________________________________________________________

./effZFit IsoMu24_eta2p1.bins 0 0 0 0 Summer12_DYJetsToLL_M-50_TuneZ2Star_smubits.root IsoMu24_eta2p1_mc   1 0 none none
./effZFit IsoMu24_eta2p1.bins 0 0 0 0 SingleMu_2012-22Jan2013_smubits.root             IsoMu24_eta2p1_data 0 0 none none

./effZFit muon_selection.bins 0 0 0 0 Summer12_DYJetsToLL_M-50_TuneZ2Star_muselbits.root muon_selection_mc   1 0 none none
./effZFit muon_selection.bins 2 1 2 2 SingleMu_2012-22Jan2013_muselbits.root             muon_selection_data 0 0 Summer12_DYJetsToLL_M-50_TuneZ2Star_muselbits.root none
