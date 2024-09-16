import pymol
from pymol import cmd, CmdException

def MUTATION(CHAINID, RESID, MUTNAME):
    cmd.wizard("mutagenesis")
    cmd.load('8f2a.pdb', '8f2a')
    cmd.get_wizard().set_mode(MUTNAME)
    cmd.get_wizard().do_select("chain " + str(CHAINID) + " and resid "+str(RESID))
    cmd.get_wizard().apply()
    cmd.save(str(CHAINID)+"_"+str(RESID)+"_structure.pdb")
    return 

MUTATION("P", 3, "ARG")
