
#qsub -cwd -V -N DG -l highp,h_data=5G,h_rt=00:30:00 run_real.sh 
qsub -cwd -V -N DG -l highp,h_data=5G,h_rt=00:30:00 run.sh 
