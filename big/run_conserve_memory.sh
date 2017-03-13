#!/bin/bash -e

#
# Memory usage
#

# parameters
arraysize=9

# submit the array job
#export jobid=$(sbatch --array=1-$arraysize run_conserve_memory.sl | awk '{print $NF}')
#echo "Job ID is $jobid"

# wait for it to finish
#echo "Waiting for job $jobid to complete"
#sleep 5
#while squeue -u ${USER} | grep ${jobid}; do sleep 30; done
#sleep 5
jobid=61101555


# get the results
echo "Processing results"
echo "Dst grid points, MaxRSS (bytes)" > results2.csv
for ((i=1; i<=$arraysize;i++)); do
    echo "Getting MaxRSS for: '${jobid}_${i}'"
    maxrss=$(sacct --format='maxrss' -j ${jobid}_${i} | tail -1 | sed 's/K//')
    gridsz=$(awk '{print $NF}' memory_${i}.csv)
    echo "MaxRSS for $i: $maxrss"
    echo "$gridsz, $maxrss" >> results2.csv
done


