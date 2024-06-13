#!/bin/bash

#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=PolymerMD_test
#################  
#a file for job output, you can check job progress
#SBATCH --output=PolymerMD.out
#################
# a file for errors from the job
#SBATCH --error=PolymerMD.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH --time=15-00:00:00
#################
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32

# Check if the environment variable 'lmp_dir' exists
if [ -z "$lmp_dir" ]; then
    echo "cannot find lmp_dir and search for LAMMPS_DIR"

    # Check if the environment variable 'LAMMPS_DIR' exists
    if [ -z "$LAMMPS_DIR" ]; then
        echo "assign path to lammps installation as lmp_dir or go through the README to install LAMMPS on your system"
    else
        export lmp_dir=$LAMMPS_DIR
    fi
else
    echo "lmp_dir exists"
fi

n_bead=500
n_chain=500
n_irrev=13600
n_rev=8000

export workdir=$(pwd) ;
export rundir=$workdir/runs/test_run ;
mkdir $rundir ; 
cd $rundir ;

echo "Move the polymer.dat data file generate by LAMMPS chain tool for the remaining lines to execute successfully"

mpirun -n 96 $lmp_dir/src/lmp_mpi -in $workdir/lammps_scripts/in.md ;

python3 $workdir/python/create_irreversible_bonds.py --n_chain $n_chain --n_bead $n_bead --n_quartic $n_irrev ;
mpirun -n 96 $lmp_dir/src/lmp_mpi -in $workdir/lammps_scripts/irreversible_elastomer_loading.md2 ;


python3 $workdir/python/create_reversible_bonds.py --n_chain $n_chain --n_bead $n_bead --n_crosslink $n_rev ;
mpirun -n 96 $lmp_dir/src/lmp_mpi -in $workdir/lammps_scripts/equil_and_activate_reversible.md2 ;
mpirun -n 96 $lmp_dir/src/lmp_mpi -in $workdir/lammps_scripts/reversible_elastomer_loading.md2 ;
