remoteToysDir="/home/dhsu/TagAndProbe/toys"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`
pwd
mkdir -p toys
tar -xvzf "${1}_tnpLibs.tgz"
tar -xvzf `printf "${1}_toy%06d_tnpData.tgz" ${13}`
echo "The following files exist in the exec dir:"
ls -l
echo "Doing the fits now"
./fitToy $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13}
#{7} is the output dir from TNP
outputTarball=`printf "${1}_toy%06d_tnpOutput.tgz" ${13}`
cd toys
ls -l
tar -zcf $outputTarball `printf "${7}/toy%06d/*.txt" ${13}`
ls -l
mv $outputTarball $remoteToysDir
ls -l
cd $remoteToysDir
ls -l
tar -xzf $outputTarball
ls -l
rm $outputTarball
ls -l
