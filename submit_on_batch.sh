#! /bin/sh
#SBATCH -J Forest 
#SBATCH -D ${rundir}
#SBATCH -o ${rundir}/log_%j.txt
export X509_USER_PROXY='/afs/hephy.at/user/m/mspanring/proxy/x509_proxy'
echo "---------------------"
echo "Grid certificate"
voms-proxy-info --all
echo "---------------------"

eval `scramv1 runtime -sh`
./convertNanoParallel.py ${file} ${channel} ${systShift} ${svfit} ${recoil} ${nevents}

mv -f ${channel}*root ${outdir}
