# library to read/write mrc files
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import MDAnalysis as mda
import argparse

## create a string of supported atoms for MDAnalysis selection
selt_="type C O N S P F"

# initialize parser
parser = argparse.ArgumentParser(prog='python check_PDB.py', description='Check PDB')
parser.add_argument('pdb', type=str, help='PDB filename')

args = parser.parse_args()

#### INPUT
# input PDB file
PDB_=vars(args)["pdb"]

# create MDA universe
u = mda.Universe(PDB_)

# do heavy atoms selection (exclude hydrogens)
atoms = u.select_atoms('not type H*')
# save PDB without hydrogens
atoms.write(PDB_.split(".")[0]+"-noH.pdb")

# number of heavy atoms
nat = len(atoms)
# number of hydrogen atoms
nH = len(u.select_atoms('type H*'))

# bfactor and occupancy lists
bfactors=[]; occupancy=[]

# loop on supported atom types 
for at in atoms.select_atoms(selt_):
    # bfactor
    bfactors.append(at.tempfactor)
    # occupancy
    occupancy.append(at.occupancy)

# fraction of C O N S P F atoms
tf = float(len(bfactors)) / float(nat)

# fraction of atoms with bfactor > 10
bf = 0.0
for b in bfactors:
    if(b>10.): bf += 1.0
bf /= float(len(bfactors))

# fraction of atoms with occupancy = 1.0
of = 0.0
for o in occupancy: 
    if(o>=1.0): of += 1.0
of /= float(len(occupancy))

# number of water molecules
nwat = len(u.select_atoms("resname HOH TIP TIP3 SOL").residues)

## PRINT OUT SOME INFO
log=open("log.PDB", "w")
log.write("PDB info\n")
log.write("%37s %s\n" % ("PDB filename :", PDB_))
log.write("%37s %d\n" % ("# heavy atoms :", nat))
log.write("%37s %d\n" % ("# hydrogen atoms :", nH))
log.write("%37s %d\n" % ("# water molecules :", nwat))
log.write("%37s %4.2lf\n" % ("fraction of C O N S P F atoms :", tf))
log.write("%37s %4.2lf\n" % ("   with bfactor > 10 :", bf))
log.write("%37s %4.2lf\n" % ("   with occupancy = 1 :", of))
