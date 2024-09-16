import pymol
from pymol import cmd, CmdException
import pandas as pd 
import os, sys
import subprocess

def MUTATION(NAME, CHAINID, RESID, MUTNAME):
    cmd.wizard("mutagenesis")
    cmd.load('8f2a.pdb', '8f2a')
    cmd.get_wizard().set_mode(MUTNAME)
    cmd.get_wizard().do_select("chain " + str(CHAINID) + " and resid "+str(RESID))
    cmd.get_wizard().apply()
    cmd.h_add()
    cmd.save("ms_"+str(NAME)+"_"+str(RESID)+"_"+str(MUTNAME)+".pdb")
    cmd.delete("8f2a")
    return


def CompareSeq(s1, s2):
    k1 = list(s1)
    k2 = list(s2)
    diff_list = []
    count=0
    for x, y in zip(k1, k2):
        count+=1
        if x == y:
            pass
        else:
            diff_list.append(str(y) +":"+ str(count))
    return diff_list



def ALLMUTATION(NAME, CHAINID, diff_list):
    cmd.wizard("mutagenesis")
    cmd.load('8f2a.pdb', '8f2a')
    cmd.h_add()
    cmd.save("8f2a_h.pdb")
    for elem in range(0, len(diff_list)):
        ele = diff_list[elem]
        mut, mut_pos = ele.split(":")[0], ele.split(":")[1]
        MUTNAME = AAlist[mut]
        RESID = mut_pos
        cmd.get_wizard().set_mode(MUTNAME)
        cmd.get_wizard().do_select("chain " + str(CHAINID) + " and resid "+str(RESID))
        cmd.get_wizard().apply()
        cmd.h_add()
        cmd.save("multi_ms_"+str(NAME)+"_"+str(RESID)+"_"+str(MUTNAME)+".pdb")
    cmd.delete("8f2a")
    return


AAlist = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS', 'E':'GLU', 'Q':'GLN', 'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS', 'M':'MET', 'F':'PHE', 'P':'PRO', 'S':'SER', 'T':'THR', 'W':'TRP', 'Y':'TYR', 'V':'VAL'}


def PROGRESS(name, chainid, ref_seq, seq):
    s1 = ref_seq
    s2 = seq
    mutations_list = CompareSeq(s1, s2)
    print(mutations_list)
    if len(mutations_list) > 0:
        for elem in range(0, len(mutations_list)):
            ele = CompareSeq(s1, s2)[elem]
            mut, mut_pos = ele.split(":")[0], ele.split(":")[1]
            print(mut, mut_pos)
            MUTATION(name, chainid, mut_pos, AAlist[mut])
        ALLMUTATION(name, chainid, mutations_list)
    else:
        pass 


df = pd.read_csv('20240903__Amylin_peptide_design_iteration_ref_sequences.csv')
df['Sequence'] = df.iloc[:,1:38].apply(lambda x: ''.join(x), axis=1)
updated_df = df[['Name', 'Sequence']][1:]

reference = updated_df['Sequence'][1]

for idx, row in updated_df[1:].iterrows():
    name = row['Name']
    sequence = row['Sequence']
    PROGRESS(name, "P", reference, sequence)
