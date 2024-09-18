#####################################################################
#
#####################################################################

import os, sys
import MDAnalysis as mda
import warnings 
warnings.filterwarnings('ignore')
from MDAnalysis.analysis.rms import RMSD 
import subprocess
import pandas as pd 



def SystermPrepare(inpfile):

    pdb = mda.Universe(inpfile)
    chain1 = pdb.select_atoms('protein and segid R and not element H')
    chain1.write('chainR.pdb')
    chain2 = pdb.select_atoms('protein and segid P and not element H')
    chain2.write('chainP.pdb')

    system = """
    source leaprc.protein.ff14SB
source leaprc.RNA.OL3
source leaprc.DNA.bsc1
source leaprc.gaff
source leaprc.water.tip3p
loadamberparams frcmod.ionsjc_tip3p

COMP1 = loadpdb chainR.pdb
COMP2 = loadpdb chainP.pdb
COMP = combine {COMP1 COMP2}
list

saveamberparm COMP comp.prmtop comp.inpcrd
saveamberparm COMP1 comp1.prmtop comp1.inpcrd
saveamberparm COMP2 comp2.prmtop comp2.inpcrd
savepdb COMP comp.pdb

source leaprc.water.tip3p
solvatebox COMP TIP3PBOX 12.0
addionsrand COMP Na+ 0
saveamberparm COMP comp-solvate.prmtop comp-solvate.inpcrd
savepdb COMP comp-solvate.pdb
quit
    """
    with open('prepare.in', 'w') as f:
        f.write(system)
    subprocess.run("conda run -n AmberTools23 tleap -f prepare.in ", shell=True)
    return 


def Minimization():
    minim = """
 Minimization of everything excluding backbone
 &cntrl
  imin = 1, maxcyc = 1000,
  ncyc = 30, ntx = 1,
  ntwe = 0, ntwr = 500, ntpr = 50,
  ntc = 2, ntf = 2, ntb = 1, ntp = 0,
  cut = 8.0,
  ntr=1, restraintmask="@CA,N,C", restraint_wt=10.
  ioutfm=1, ntxo=2,
 /
    """
    with open("01min.in", 'w') as m:
        m.write(minim)
    subprocess.run("conda run -n AmberTools23 sander -O -i 01min.in -o min.out -p comp-solvate.prmtop -c comp-solvate.inpcrd -r min.rst7 -inf min.info -ref comp-solvate.inpcrd -x min.nc ", shell=True)
    return 

def MDSimulation():
    mdsim = """
10-mer DNA MD in-vacuo, 12 angstrom cut off
 &cntrl
  imin = 0, ntb = 0,
  igb = 0, ntpr = 100, ntwx = 100,
  ntt = 3, gamma_ln = 1.0,
  tempi = 300.0, temp0 = 300.0
  nstlim = 1000, dt = 0.001,
  cut = 12.0,
  ntr=1, restraintmask="@CA,N,C", restraint_wt=10.
 /
    """
    with open('03mdsim.in', 'w') as mds:
        mds.write(mdsim)
    subprocess.run("conda run -n AmberTools23 sander -O -i 03mdsim.in -o mdsim.out -p comp-solvate.prmtop -c min.rst7 -r md.rst7 -inf md.info -ref min.rst7 -x mdcrd.nc ", shell=True)
    return 

def LowestRMSDStructure(reference, prmtop, trajectory):

    inp = mda.Universe(prmtop, trajectory)
    mobile = inp.atoms

    ref_pdb = mda.Universe(reference).atoms 

    R = RMSD(mobile, ref_pdb, select="backbone").run()

    ref_rmsd = min(R.rmsd.T[2])

    thisdict = {}
    for ts, rmsd in zip(R.rmsd.T[1], R.rmsd.T[2]):
        thisdict[rmsd] = ts

    lowestRmsdStr = list(thisdict).index(ref_rmsd)

    receptor = inp.select_atoms("protein")
    receptor.write('lowestrmsdstructure.pdb', frames=inp.trajectory[[lowestRmsdStr]])

    return 


def MMGBSAForSingle():
    mmgbsa = """
Input file for running PB and GB
&general
   endframe=99, verbose=1,
#   entropy=1,
/
&gb
  igb=5, saltcon=0.100
/
&decomp
  idecomp=1,
  dec_verbose=1,
/
    """
    with open('mmgbsa0.in', 'w') as m:
        m.write(mmgbsa)
    subprocess.run("conda run -n AmberTools23 MMPBSA.py -O -i mmgbsa0.in -o FINAL_RESULTS_MMPBSA.dat -sp comp-solvate.prmtop -cp comp.prmtop -rp comp1.prmtop -lp comp2.prmtop -y min.rst7 ", shell=True)
    return



def MMGBSAwithDecomp():
    mmgbsa = """
Input file for running PB and GB
&general
   endframe=99, verbose=1,
#   entropy=1,
/
&gb
  igb=5, saltcon=0.100
/
&decomp
  idecomp=1,
  dec_verbose=1,
/
    """
    with open('mmgbsa.in', 'w') as m:
        m.write(mmgbsa)
    subprocess.run("conda run -n AmberTools23 MMPBSA.py -O -i mmgbsa.in -o FINAL_RESULTS_MMPBSA.dat -sp comp-solvate.prmtop -cp comp.prmtop -rp comp1.prmtop -lp comp2.prmtop -y mdcrd.nc ", shell=True)
    return



# Total:370:407, backbone:779:816, sidechain:1188:1234 
def DecompostionProcessing(inpfile, firstline, lastline, specific):
    df = pd.read_csv(inpfile, skiprows=6)
    df1 = df.iloc[firstline:lastline, 1:10]
    df2 = df1.drop(['Internal', 'Unnamed: 3', 'Unnamed: 4'], axis=1)
    df3 = df2.rename(columns={'van der Waals': 'E_vdW', 'Unnamed: 6': 'E_Elec', 'Unnamed: 7': 'E_Polar','Electrostatic': 'E_nonPolar', 'Unnamed: 9': 'E_Total' })
    df4a = df3.T
    df4a.columns = df4a.iloc[0]
    df4a = df4a[1:]
    df4a.to_csv(str(inpfile.split('.')[0])+"_"+str(specific)+".csv")
    return



def main():
    SystermPrepare(inpfile)
    Minimization()
    MDSimulation()
    LowestRMSDStructure('comp-solvate.pdb', 'comp-solvate.prmtop', 'mdcrd.nc')
    MMGBSAForSingle()
    rename = "mv FINAL_DECOMP_MMPBSA.dat "+ str(inpfile.strip().split(".")[0])+"_FINAL_DECOMP_MMPBSA_min.dat"
    os.system(rename)
    ddg_for_min = subprocess.check_output("grep DELTA FINAL_RESULTS_MMPBSA.dat | awk '{print $3}' | tail -n1", shell=True)
    rename = "mv FINAL_RESULTS_MMPBSA.dat "+ str(inpfile.strip().split(".")[0])+"_FINAL_RESULTS_MMPBSA_min.dat"
    os.system(rename)

    DecompostionProcessing(str(inpfile.strip().split(".")[0])+"_FINAL_DECOMP_MMPBSA_min.dat", 370, 407, "total")
    DecompostionProcessing(str(inpfile.strip().split(".")[0])+"_FINAL_DECOMP_MMPBSA_min.dat", 779, 816, "bkb")
    DecompostionProcessing(str(inpfile.strip().split(".")[0])+"_FINAL_DECOMP_MMPBSA_min.dat", 1188, 1234, "sidechain")
#
    MMGBSAwithDecomp()
    rename = "mv FINAL_DECOMP_MMPBSA.dat "+ str(inpfile.strip().split(".")[0])+"_FINAL_DECOMP_MMPBSA.dat"
    os.system(rename)
#
    DecompostionProcessing(str(inpfile.strip().split(".")[0])+"_FINAL_DECOMP_MMPBSA.dat", 370, 407, "total")
    DecompostionProcessing(str(inpfile.strip().split(".")[0])+"_FINAL_DECOMP_MMPBSA.dat", 779, 816, "bkb")
    DecompostionProcessing(str(inpfile.strip().split(".")[0])+"_FINAL_DECOMP_MMPBSA.dat", 1188, 1234, "sidechain")
#
    rename = "mv mdcrd.nc "+ str(inpfile.strip().split(".")[0])+"_mdcrd.nc"
    os.system(rename)
    rename = "mv comp-solvate.prmtop "+ str(inpfile.strip().split(".")[0])+"_comp-solvate.prmtop"
    os.system(rename)
    rename = "mv comp-solvate.pdb "+ str(inpfile.strip().split(".")[0])+"_comp-solvate.pdb"
    os.system(rename)
    rename = "mv lowestrmsdstructure.pdb "+ str(inpfile.strip().split(".")[0])+"_lowestrmsdstructure.pdb"
    os.system(rename)
#
    subprocess.run("rm comp* chain* leap.log _MMPBSA* min.* mdsim.* *.in md.*", shell=True)
    outfile = open('mmgbsa_freeenergies.txt', 'a')
    ddg = subprocess.check_output("grep DELTA FINAL_RESULTS_MMPBSA.dat | awk '{print $3}' | tail -n1", shell=True)
    rename = "mv FINAL_RESULTS_MMPBSA.dat "+ str(inpfile.strip().split(".")[0])+"_FINAL_RESULTS_MMPBSA.dat"
    os.system(rename)
    outfile.write(str(inpfile)+"\t"+str(ddg_for_min.decode().strip())+"\t"+str(ddg.decode().strip()))
    outfile.close()
#
    os.system("mv *min* MMGBSA-min")
    os.system("mv *rmsdstructure.pdb RMSD")
    os.system("mv *csv *dat MMGBSA-decomp")
    os.system("mv "+str(inpfile.strip().split(".")[0]) +"_* INP")
    return 

if __name__ == "__main__":
    if len(sys.argv) > 1:
        inpfile = sys.argv[1]
        dirt_list = ["INP", "RMSD", "MMGBSA-decomp", "MMGBSA-min"]
        for dirt in dirt_list:
            if os.path.exists(dirt):
                os.mkdir(dirt)
        main()
    else:
        print("Usage: python3 PeptideOptimizationWithAnalysis.py {inputpdb}")
