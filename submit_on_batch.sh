#! /bin/sh
#SBATCH -J Forrest 
#SBATCH -D ${rundir}
#SBATCH -o ${rundir}/log_%j.txt
export X509_USER_PROXY='/afs/hephy.at/user/m/mspanring/proxy/x509_proxy'
echo "---------------------"
echo "Grid certificate"
voms-proxy-info --all
echo "---------------------"

echo "Good night..."
date
sleep ${sleeping}
echo "Good morning..."
date
eval `scramv1 runtime -sh`
./convertNanoParallel.py ${file} ${channel} ${systShift} ${svfit} ${recoil} ${nevents}
date
mv -f ${channel}*root ${outdir}
