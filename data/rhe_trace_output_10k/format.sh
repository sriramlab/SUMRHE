for i in {0..99};
do
    echo 'NSAMPLE,NSNPS,NBLKS' > run_${i}.MN
    echo '10059,454207,1000' >> run_${i}.MN
    echo 'TRACE,NSNPS_JACKKNIFE' >> run_${i}.trace
    awk 'NR % 4 == 1 {print $1}' singleout_10k_run_${i}.txt.normaleq.txt | head -n 1000 > tmp.txt
    paste -d ',' tmp.txt nsnps.txt >> run_${i}.trace
    rm tmp.txt
done
