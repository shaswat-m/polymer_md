# POLYMER_MD
Helper scripts to create an equilibrated polymer melt to which irreversible or reversible crosslinks can be added.

## Cite
* Modeling shortest paths in polymeric networks using spatial branching processes (https://www-sciencedirect-com.stanford.idm.oclc.org/science/article/pii/S0022509624001029)
```
@article{zhang2024modeling,
  title={Modeling shortest paths in polymeric networks using spatial branching processes},
  author={Zhang, Zhenyuan and Mohanty, Shaswat and Blanchet, Jose and Cai, Wei},
  journal={Journal of the Mechanics and Physics of Solids},
  volume={187},
  pages={105636},
  year={2024},
  publisher={Elsevier}
}
```

* Network evolution controlling strain-induced damage and self-healing of elastomers with dynamic bonds (https://arxiv.org/html/2401.11087v1)
```
@article{yin2024network,
  title={Network evolution controlling strain-induced damage and self-healing of elastomers with dynamic bonds},
  author={Yin, Yikai and Mohanty, Shaswat and Cooper, Christopher B and Bao, Zhenan and Cai, Wei},
  journal={arXiv preprint arXiv:2401.11087},
  year={2024}
}
```

## Preparing the polymer network
Set work directory
```commandline
export workdir=$(pwd)
```
### Install LAMMPS
Here are the steps to follow:
* Clone the `release` version of LAMMPS from GitHub
* Install the necessary packages for CGMD simulations: `MANBODY`, `MOLECULE`, `MC`, and `DPD-BASIC`
```commandline
git clone -b release https://github.com/lammps/lammps.git mylammps ;
cd mylammps/src ;
make yes-manybody ;
make yes-molecule ;
make yes-mc ; 
make yes-dpd-basic
make mpi ;
```
* Add the cloned repository to the `lmp_dir` environment variable.
```commandline
export lmp_dir=$workdir/mylammps
```
* Setup runs directory
```commandline
export rundir=$workdir/runs/test_run ;
mkdir $rundir ; 
```

### Create initial configuration

* Use the LAMMPS `chain` tool to create the initial configuration. Edit the `$lmp_dir/tools/def.chain` to choose the number of polymer chains (`number of chains`) and the length (degree of polymerization - `monomers/chain`) of the polymer chains before executing the following cell to create the `chain` executable. For production level runs, we use 500 chains with 500 beads.
```commandline
cd $lmp_dir/tools ;
gfortran -o chain chain.f90  ;
```
* Run the `chain` executable to generate the initial input configuration.
```commandline
cd $rundir ;
$lmp_dir/tools/chain < $lmp_dir/tools/def.chain > polymer.dat ;
```
* Run the lammps script `in.md` to create an equilibrated polymer melt (run on more CPUs if running on a cluster)
```commandline
mpirun -n 1 $lmp_dir/src/lmp_mpi -in $workdir/lammps_scripts/in.md ;
```
* The equilibrated polymer melt, `polymer_melt.data`, will be created at the end of this simulation for the assignment of the `irreversible` or `reversible` cross-links.


### The irreversible cross-linked polymer network simulation
* We begin by assigning irreversible cross-links to the polymer melt. This is done by identifying close enough beads that lie on different chains (avoids self-loops) as candiddate cross-linking sites. Depedning on the number of cross-links to be assigned the script chooses from these candidate sites either randomly or through a bunch of other options available in `create_irreversible_bonds.py`
```commandline
python3 $workdir/python/create_irreversible_bonds.py ;
```
Note, the file I/O might appear rudimentary and can be replaced with the Atomic Simulation Environment (ASE) for ease of usage. Current functions have been retained in their legacy versions.
* We can now equilibrate the system for a small duration before loading it uniaxially (if that is the desired simulation)
```commandline
mpirun -n 1 $lmp_dir/src/lmp_mpi -in $workdir/lammps_scripts/irreversible_elastomer_loading.md2 ;
```

### The reversible cross-linked polymer network simulation
* We begin by assigning reversible cross-links to the polymer melt. This is done by dissolving chains into monomers to serve as the emtal ions that will form the reversible cross-links. Depending on the ligand-metal ion ratio, we assign beads along the remaining polymer backbones as the ligands. These ligands can either be spaced uniformly (default) or randomly using the options available in `create_reversible_bonds.py`
```commandline
python3 $workdir/python/create_reversible_bonds.py ;
```
* We can now start activating the reversible cross-link potential gradually to not disturb the system aggressively.
```commandline
mpirun -n 1 $lmp_dir/src/lmp_mpi -in $workdir/lammps_scripts/equil_and_activate_reversible.md2 ;
```
* We can now equilibrate the system for a small duration before loading it uniaxially (if that is the desired simulation)
```commandline
mpirun -n 1 $lmp_dir/src/lmp_mpi -in $workdir/lammps_scripts/reversible_elastomer_loading.md2 ;
```

## Authors and acknowledgment
Authored by Shaswat Mohanty, and Wei Cai. Development of CGMD scripts for the simulations of DPNs were derived from earlier works of the group carried out by Yikai Yin.

## Project status
Development on the project will be sporadically updated by Shaswat Mohanty, Myng Chul Kim and Wei Cai.
