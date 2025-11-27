#!/bin/bash

# (C) C. Vriend - Aumc ANW/Psy - June '23
# c.vriend@amsterdamumc.nl

# Combined wrapper and pipeline script for DWI preprocessing and tractography pipeline
# Automatically determines the number of subjects and submits jobs with the correct dependencies using SLURM arrays.


# input variables and paths
scriptdir=/scratch/anw/cvriend/dwi-tractography-pipeline/dwi-tract4mega
bidsdir=/data/anw/anw-archive/NP/imaging-samples/MDD_MOTAR
workdir=~/my-scratch/dwi-preproc
outputdir=/data/anw/anw-archive/NP/projects/archive_MOTAR/derivatives
freesurferdir=/data/anw/anw-archive/NP/projects/archive_MOTAR/derivatives/freesurfer


# How many in parallel?
simul=2
# Run NODDI? 1/0 = yes/no
noddi=0

# Determine subjects
cd ${bidsdir}
ls -d sub-motar5*/ | sed 's:/.*::' > ${scriptdir}/subjects.txt
nsubj=$(wc -l < ${scriptdir}/subjects.txt)

cd ${scriptdir}

# Submit the SLURM array
#sbatch --array=1-${nsubj}%${simul} \
sbatch --array=7-70%2 \
 --export=ALL,scriptdir=${scriptdir},bidsdir=${bidsdir},workdir=${workdir},outputdir=${outputdir},freesurferdir=${freesurferdir},noddi=${noddi},subjects=${scriptdir}/subjects.txt \
  ${scriptdir}/dwi-01-pipeline_sarray.sh

# Clean up
#rm ${scriptdir}/subjects.txt
echo "SLURM array submitted for ${nsubj} subjects with a maximum of ${simul} running in parallel."
