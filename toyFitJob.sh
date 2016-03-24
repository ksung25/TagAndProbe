remoteToysDir="/home/dhsu/TagAndProbe/toys"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`
mkdir -p toys
#tar -xvzf "${1}_tnpLibs.tgz"
tar -xvzf `printf "${1}_toy%06d_tnpData.tgz" ${13}`
#echo "The following files exist in the exec dir:"
#ls -l
echo "Doing the fits now: ./fitToy $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13}"
./fitToy $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13}
#{7} is the output dir from TNP

#rm fitToy ${2} CEffUser1D.o  CEffUser2D.o  CEffZFitter.o  CPlot.o  KStyle.o  RooCMSShape.o  RooVoigtianShape.o
outputTarball=`printf "${1}_toy%06d_tnpOutput.tgz" ${13}`
cd toys
tar -zcf $outputTarball `printf "${7}/toy%06d/*.txt" ${13}`
mv $outputTarball $remoteToysDir
cd $remoteToysDir
tar -xzf $outputTarball
rm $outputTarball
