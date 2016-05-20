#!/bin/bash
# needs comments
function splitPath {
  local IFS=/
  local pieces
  # Disable glob expansion
  set -f
  pieces=( $@ ) 
  set +f
  #printf '%d\n' "${#pieces[@]}"
  #printf '%s\n' "${pieces[@]}"
  echo ${pieces[${#pieces[@]}-1]}
}

work_dir=/home/dhsu/TagAndProbe/toys
repo_dir=/home/dhsu/TagAndProbe
tarballDir="/scratch5/dhsu/tnpToyData"

#template_erfcexp electron_scalefactors.bins MCTemplateConvGaussian ErfcExpo MCTemplateConvGaussian ErfcExpo 2016-02-26_74x_template_erfcexp/SingleElectron_BaselineToTight_electronTnP 0 1 ~/leptonScaleFactors/root/DYJetsToLL_74x_aMCatNLO_genMatching_BaselineToTight_electronTnP.root - - 0 999
cd $repo_dir
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`
./fitToy $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} -1 
cd $work_dir

to_run="/home/dhsu/bin/condor-run ./toyFitJob.sh"
#libTarball="/tmp/${1}_tnpLibs.tgz"
#echo "Making tarball for TnP executable, binning config, and libraries: \"$libTarball\""

mkdir -p ${tarballDir}/${7}

libs="${repo_dir}/toyFitJob.sh ${repo_dir}/fitToy ${repo_dir}/${2} ${repo_dir}/CEffUser1D.o  ${repo_dir}/CEffUser2D.o  ${repo_dir}/CEffZFitter.o  ${repo_dir}/CPlot.o  ${repo_dir}/KStyle.o  ${repo_dir}/RooCMSShape.o  ${repo_dir}/RooVoigtianShape.o ${repo_dir}/libCEffUser1D.so  ${repo_dir}/libCEffUser2D.so  ${repo_dir}/libCEffZFitter.so  ${repo_dir}/libCPlot.so  ${repo_dir}/libKStyle.so  ${repo_dir}/libRooCMSShape.so  ${repo_dir}/libRooVoigtianShape.so"
for i in `seq ${13} ${14}`
do
    dataTarball=`printf "${tarballDir}/${7}/${1}_toy%06d_tnpData.tgz" $i`
    selection_basename=`splitPath $7`
    outputTarball=`printf "${1}_${selection_basename}_toy%06d_tnpOutput.tgz" $i`
    echo "Making tarball for the toy input files, fit result reference directories, and any existing MC templates: \"$dataTarball\""
    echo "Output tarball will arrive at: \"${work_dir}/$outputTarball\""
    cd $repo_dir
    tar -zcf $dataTarball `printf "toys/${7}/toy%06d/*.dat" $i` `ls templates/${7}/*.root 2> /dev/null` `ls ${11}/plots/fitres*.txt 2> /dev/null` `ls ${12}/plots/fitres*.txt 2> /dev/null`
    cd $work_dir

    args="--auxiliary-input $dataTarball $libs --output toys/$outputTarball "
    args="${args} --task-name $1"
    args="${args} --job-args \"$1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} $i\" "
    echo "$to_run $args"
    eval "$to_run $args"
done
