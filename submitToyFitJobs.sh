#!/bin/bash
# needs comments

tarballDir="/scratch5/dhsu/tnpToyData"

#template_erfcexp electron_scalefactors.bins MCTemplateConvGaussian ErfcExpo MCTemplateConvGaussian ErfcExpo 2016-02-26_74x_template_erfcexp/SingleElectron_BaselineToTight_electronTnP 0 1 ~/leptonScaleFactors/root/DYJetsToLL_74x_aMCatNLO_genMatching_BaselineToTight_electronTnP.root - - 0 999

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`
./fitToy $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} -1 

to_run="/home/dhsu/bin/condor-run toyFitJob.sh"
#libTarball="/tmp/${1}_tnpLibs.tgz"
#echo "Making tarball for TnP executable, binning config, and libraries: \"$libTarball\""

mkdir -p ${tarballDir}/${7}

libs="fitToy ${2} CEffUser1D.o  CEffUser2D.o  CEffZFitter.o  CPlot.o  KStyle.o  RooCMSShape.o  RooVoigtianShape.o"
for i in `seq ${13} ${14}`
do
    dataTarball=`printf "${tarballDir}/${7}/${1}_toy%06d_tnpData.tgz" $i`
    echo "Making tarball for the toy input files, fit result reference directories, and any existing MC templates: \"$dataTarball\""
    tar -zcf $dataTarball `printf "toys/${7}/toy%06d/*.dat" $i` `ls templates/${7}/*.root 2> /dev/null` `ls ${11}/plots/fitres*.txt 2> /dev/null` `ls ${12}/plots/fitres*.txt 2> /dev/null`
    args="--auxiliary-input $dataTarball $libs"
    args="${args} --task-name $1"
    args="${args} --job-args \"$1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} $i\" "
    echo "$to_run $args"
    eval "$to_run $args"
    #rm $dataTarball
done
#rm $libTarball
