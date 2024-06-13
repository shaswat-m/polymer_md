# Set the number of processors for parallel execution
#processors 2 2 6

# Choose LJ units and specify bond atom style
units lj
atom_style bond

# Set special LJ/Coulombic bonds
special_bonds lj/coul 0 1 1

# Read data file containing polymer information
read_data polymer.dat 

# Set up neighbor lists for pairwise interactions
neighbor 0.4 bin
neigh_modify every 1 delay 1
comm_modify vel yes

# Define FENE bond style and coefficients
bond_style fene
bond_coeff 1 30.0 1.5 1.0 1.0

# Compute local properties and dump chain configurations
compute 1 all property/local batom1 batom2 btype
dump mydump all custom 1000 chain.xyz id type xu yu zu
dump bonddump all local 1000 bond.xyz index c_1[1] c_1[2] c_1[3]

# Set simulation timestep and output thermo data every 1000 steps
timestep 0.01
thermo 10000

# Set up DPD pairwise interactions
pair_style dpd 1.0 1.0 122347
pair_coeff 1 1 25 4.5 1.0

# Set initial velocities and perform NVE dynamics for 10000 steps
velocity all create 1.0 17786140
fix 1 all nve
run 100000

# Adjust pair coefficients and velocities for subsequent runs
pair_coeff * * 50.0 4.5 1.0
velocity all create 1.0 15086120
run 50000
pair_coeff * * 100.0 4.5 1.0
velocity all create 1.0 15786120 
run 50000
pair_coeff * * 150.0 4.5 1.0
velocity all create 1.0 15486120
run 50000
pair_coeff * * 200.0 4.5 1.0
velocity all create 1.0 17986120
run 100000
pair_coeff * * 250.0 4.5 1.0
velocity all create 1.0 15006120
run 100000
pair_coeff * * 500.0 4.5 1.0
velocity all create 1.0 15087720
run 100000
pair_coeff * * 1000.0 4.5 1.0
velocity all create 1.0 15086189
run 100000

# Set up hybrid pair style with LJ/Cut and DPD/Tstat
pair_style hybrid/overlay lj/cut 1.122462 dpd/tstat 1.0 1.0 1.122462 122347
pair_modify shift yes
pair_coeff * * lj/cut 1.0 1.0 1.122462
pair_coeff * * dpd/tstat 4.5 1.122462

# Perform runs with different velocity initializations
velocity all create 1.0 1508612013
run 10000
velocity all create 1.0 15086121
run 10000
velocity all create 1.0 15086111
run 10000
velocity all create 1.0 15086125
run 10000

# Write final data and restart files
write_data polymer_melt.data
write_restart my_restart_lj
