# gro file
gro=$1
# pdb file
pdb=$2

# rename old files
mv $gro old_$gro
mv $pdb old_$pdb

# create a fake tpr file
gmx_mpi grompp -f ../2-Equilibration/0-em-steep.mdp -c old_$gro -p topol.top -o fake.tpr

# convert old gro to new gro and pdb
echo 0 | gmx_mpi trjconv -f old_$gro -o $gro -s fake.tpr
echo 0 | gmx_mpi trjconv -f old_$gro -o $pdb -s fake.tpr

# clean
rm fake.tpr
