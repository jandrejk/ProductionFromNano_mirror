#! /bin/sh
#SBATCH -J ZoomZoom 
#SBATCH -D ${rundir}
#SBATCH -o ${rundir}/%j.txt
export X509_USER_PROXY='/afs/hephy.at/user/m/mspanring/proxy/x509_proxy'
echo "---------------------"
echo "Grid certificate"
voms-proxy-info --all
echo "---------------------"

eval `scramv1 runtime -sh`
./convertNanoParallel.py ${channel} ${file} ${svfit} ${recoil}

mv -f HTT*root ${outdir}
rm *.h *.cxx *.C *.cc *.py *.pyc *.root *.d *.so *.pcm
