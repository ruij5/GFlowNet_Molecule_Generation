from rdkit.Chem import BRICS
from rdkit import Chem
from rdkit.Chem import Draw
import re
import csv

training_set_list = []
with open("./mols/data/training_set.csv", 'r') as file:
  csvreader = csv.reader(file)
  for row in csvreader:
    training_set_list.append(row)

# Store all the molecules from the training_set.csv
smi_list = []
for i in training_set_list:
    smi_list.append(i[0])
# print(smi_list[0]) # 'CSCCC(NCc1ccc2c(c1)CN(C)C2)c1nnc2ccccn12'
# print(len(smi_list)) # 450000



def identify_free_bond_pos(smi,matched_atom_idxs,link_value,bond_lst):
    if len(matched_atom_idxs)>1:
        for atom_idxs in matched_atom_idxs:
            if sorted(list(set(atom_idxs).intersection(set(bond_lst))))==sorted(list(bond_lst)):
                unique_atom_idx=atom_idxs
    elif len(matched_atom_idxs)==1:
        unique_atom_idx=matched_atom_idxs[0]
    else:
        print("no match:{len(matched_atom_idxs)}")
        return None
    
    print(f"INPUT\nfrag smi:{smi};unique_atom_idx:{unique_atom_idx};link_value:{link_value};bond_lst:{bond_lst}")
    smi_char_lst=re.split('\[|\]|\(|\)',smi)
    new_smi_char_lst=[]
    
    num_lst=[str(i) for i in list(range(9))]
    for smi_char in smi_char_lst:
        
        if not smi_char=='':
            if '*' in smi_char:
                new_smi_char_lst.append(smi_char)
            else:
                for a in smi_char:
                    if not a in num_lst:
                        new_smi_char_lst.append(a)

    mol=Chem.MolFromSmiles(smi)
    atoms_num=len(mol.GetAtoms())
    bond_pos=-1

    print(f"atom num==len(new_smi_char_lst):{atoms_num==len(new_smi_char_lst)}")
    n=0
    for cdx,atom in enumerate(mol.GetAtoms()):
        if (atom.GetSymbol() == "*"):
            if new_smi_char_lst[cdx]==link_value+"*":
                bond=unique_atom_idx[cdx]
                if bond in bond_lst:
                    bond_pos=atom.GetIdx()
                    if not bond_pos==0:
                        bond_pos=bond_pos-n-1
                    print(f"matched bond_pos:{bond_pos}")
            n+=1
    if bond_pos==-1:
        print(f"\nplease check the order of atom index...")
        print(f"atom list:{new_smi_char_lst}; link index list:{unique_atom_idx}; link_bond:{bond_lst}\n")
    return bond_pos