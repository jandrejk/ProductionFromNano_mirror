#! /bin/sh
#SBATCH -J ${samplename} 
#SBATCH -D ${rundir}
#SBATCH -o ${rundir}/log_%j.txt
export X509_USER_PROXY=${rundir}/proxy/x509_proxy
echo "---------------------"
echo "Grid certificate"
voms-proxy-info --all
echo "---------------------"

echo "Dummy"
echo ${cell}

echo "Good night..."
date
sleep ${sleeping}
echo "Good morning..."
date
eval `scramv1 runtime -sh`
./convertNanoParallel.py
date
echo ${outdir}
chmod 777 ${channel}*root
mv -f ${channel}*root ${outdir}
