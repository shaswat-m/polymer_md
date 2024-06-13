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

## Create initial configuration

* Use the LAMMPS `chain` tool to create the initial configuration. Edit the `$lmp_dir/tools/def.chain` to choose the number of polymer chains (`number of chains`) and the length (degree of polymerization - `monomers/chain`) of the polymer chains before executing the following cell to create the `chain` executable.
```commandline
cd $lmp_dir/tools ;
gfortran -o chain chain.f90  ;
```
* Run the `chain` executable to generate the initial input configuration.
```commandline
cd $rundir ;
$lmp_dir/tools/chain < $lmp_dir/tools/def.chain > polymer.dat ;
```
* Run the lammps script `in.md` to create an equilibrated polymer melt
```commandline
mpirun -n 1 $lmp_dir/src/lmp_mpi -in $workdir/lammps_scripts/in.md ;
```

## Authors and acknowledgment
Authored by Shaswat Mohanty, and Wei Cai. Development of CGMD scripts for the simulations of DPNs were derived from earlier works of the group carried out by Yikai Yin.

## Project status
Development is currently ongoing and is intended to be the dissertation project of Shaswat Mohanty.
