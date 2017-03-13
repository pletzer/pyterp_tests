#!/bin/bash -le
#SBATCH -J pyterp
#SBATCH -A nesi99999          # Project Account
#SBATCH --time=0:59:00        # Walltime HH:MM:SS
#SBATCH --mem-per-cpu=64G     # Memory
#SBATCH --ntasks=1            # number of tasks
#SBATCH --cpus-per-task=1     # number of threads
#SBATCH --partition=debug
#SBATCH --constraint=ib
##SBATCH --array=1-1

echo "Running"

# load modules
ml ESMF Python-Iris
source /projects/nesi99999/ANTS-regridding/local/pythonpath.sh

# some file names
srcfn=/projects/nesi99999/ANTS-regridding/interp_data/coords_CF_ORCA12_GO6.nc
dstfn=dst_$SLURM_ARRAY_TASK_ID.nc
proffn=profile_$SLURM_JOBID_$SLURM_ARRAY_TASK_ID.h5
resultfn=memory_${SLURM_ARRAY_TASK_ID}.csv
rm -f $resultfn

# size of grid
src_nj=3606
src_ni=4322
dst_nj=$(python -c "print 2**(${SLURM_ARRAY_TASK_ID} - 1) * 10")
dst_ni=$(python -c "print 2 * ${dst_nj}")
gridsize=$(python -c "print ${src_nj} * ${src_ni} * ${dst_nj} * ${dst_ni}")
nptsdst=$(python -c "print ${dst_nj} * ${dst_ni}")
echo "Source grid: $src_nj x $src_ni"
echo "Dest grid: $dst_nj x $dst_ni"
echo "Grid n*m: $gridsize"
echo "Number of points in dst grid: $nptsdst"

# generate grid/field
echo "Generating field"
srun --ntasks=1 python generate_field.py --dst_nj $dst_nj --dst_ni $dst_ni --dst_file $dstfn

# run ESMF
echo "Running ESMF"
#srun --profile=task --acctg-freq=5 python esmf_conserve.py --dst_file $dstfn --src_file $srcfn
srun python esmf_conserve.py --dst_file $dstfn --src_file $srcfn
echo "Finished ESMF"

# create hdf5 profile file
#sh5util -j $SLURM_JOBID -o $proffn

# get memory usage
#ml Python/2.7.11-intel-2015a  # need python with h5py
#maxrss=$(srun --ntasks=1 python get_maxrss.py ${proffn})
#echo "$gridsize, $nptsdst, $maxrss" > $resultfn
#echo "Result: $gridsize, $nptsdst, $maxrss"

# write out grid size
echo "$gridsize, $nptsdst" > $resultfn
echo "Result: $gridsize $nptsdst"
