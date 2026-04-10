#!/bin/bash
#SBATCH --job-name=wrappertest
#SBATCH --mem=20G
#SBATCH --partition=luna-cpu-long
#SBATCH --qos=anw-cpu
#SBATCH --cpus-per-task=8
#SBATCH --time=00-12:00:00
#SBATCH --nice=2000
#SBATCH --output=%x_%A.log


# # check eddycpu_container.sh in container
# apptainer exec --bind /scratch/anw/cvriend /scratch/anw/cvriend/TractoFriend_nofast2.sif bash -c '/home/anw/cvriend/my-scratch/dwi-tractography-pipeline/dwi-tract4mega/dwi-02b-eddyCPU_container.sh \
# -i /scratch/anw/cvriend/dwi-test/bids -o /scratch/anw/cvriend/dwi-test/derivatives \
# -w /scratch/anw/cvriend/dwi-test/work -s sub-916015 -z ses-T0 -m default'

apptainer exec --bind /scratch/anw/cvriend,/data/anw/anw-work,/data/anw/anw-archive /scratch/anw/cvriend/TractoFriend_nofast2.sif bash -c '/home/anw/cvriend/my-scratch/dwi-tractography-pipeline/dwi-tract4mega/dwi-preproc_wrapper.sh /scratch/anw/cvriend/spec.json'
