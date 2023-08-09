import pandas as pd

import pickle
from rdkit.Chem import BRICS

# Used to read info.pkl (This is the output of gflownet_activelearning.py)
# with open("mols/data/info.pkl", 'rb') as file:
#     object = pickle.load(file)
# print(object['train_metrics'])

# #print(len(object['smis'][0]))


#df = pd.read_pickle('mols/info.pkl')
#df.to_csv('mols/info.pkl', index=False)


from rdkit import Chem
from rdkit.Chem import Draw
# import matplotlib.pyplot as plt
# mols=[]
# for smi in object['smis'][0]:
#     mol=Chem.MolFromSmiles(smi)
#     mols.append(mol)
#     if len(mols)==10:
#         break

# fig=Draw.MolsToGridImage(mols,molsPerRow=5,subImgSize=(500,500))




# read h5
# columns = ["smiles", "dockscore", "blockidxs", "slices", "jbonds", "stems"]
# store = pd.HDFStore("mols/data/docked_mols.h5", 'r')
# df = store.select('df')
# df.to_csv('file_name.csv')
# print(df)

# smis=["CC(=O)C1CCN(c2c[nH]c3ccccc23)CC1","C1CCNCC1","CC=O","c1ccc2[nH]ccc2c1"]#"C1=CCC(C2(C34CC(C3)OC4c3nccs3)=CCNCCC2)C=C1"
# mols=[Chem.MolFromSmiles(smi) for smi in smis]
# for mol in mols:
#     for i, atom in enumerate(mol.GetAtoms()):
#         atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))

# fig=Draw.MolsToGridImage(mols,subImgSize=(500,500))
# fig.save("./mols.png")


# smis=["O=c1nc(Cl)ccn1-c1[nH]ncc1C1CC=CCC1","O=c1nccc[nH]1","Cl","c1cn[nH]c1", "C1=CCCCC1"]#"C1=CCC(C2(C34CC(C3)OC4c3nccs3)=CCNCCC2)C=C1"
# mols=[Chem.MolFromSmiles(smi) for smi in smis]
# for mol in mols:
#     for i, atom in enumerate(mol.GetAtoms()):
#         atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))

# fig=Draw.MolsToGridImage(mols,subImgSize=(500,500))
# fig.save("./mols_new.png")



# smis=["CC(=O)C1CCN(c2c[nH]c3ccccc23)CC1","C1CCNCC1", "CC=O", "c1ccc2[nH]ccc2c1"]#"C1=CCC(C2(C34CC(C3)OC4c3nccs3)=CCNCCC2)C=C1"
# mols=[Chem.MolFromSmiles(smi) for smi in smis]
# for mol in mols:
#     for i, atom in enumerate(mol.GetAtoms()):
#         atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))

# fig=Draw.MolsToGridImage(mols,subImgSize=(500,500))
# fig.save("./mols_new.png")



smi="*c1nc2c(ncn2*)c(=O)[nH]1"
mol=Chem.MolFromSmiles(smi)
for i, atom in enumerate(mol.GetAtoms()):
    # print(atom.GetSymbol())
    # print(atom.GetIdx())
    atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))

fig=Draw.MolToImage(mol,subImgSize=(500,500))
fig.save("./mol.png")



# smi="CC(=O)C1CCN(c2c[nH]c3ccccc23)CC1"#"C1=CCC(C2(C34CC(C3)OC4c3nccs3)=CCNCCC2)C=C1"
# m = Chem.MolFromSmiles(smi)
# res = list(BRICS.BRICSDecompose(m))
# smis = []
# smis.append(smi)
# for i in sorted(res):
#     print(i)
#     smis.append(i)
# mols=[Chem.MolFromSmiles(smi) for smi in smis]
# for mol in mols:
#     for i, atom in enumerate(mol.GetAtoms()):
#         atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
# res = list(BRICS.FindBRICSBonds(m))
# print(res)
# fig=Draw.MolsToGridImage(mols,subImgSize=(500,500))
# fig.save("./mols.png")


# new_smis=[]
# for smi in smis:
#     mol=Chem.MolFromSmiles(smi)
#     # new_smi=Chem.MolToSmiles(mol)
#     new_smi=Chem.MolToSmiles(mol,isomericSmiles=False)
#     # new_smi=new_smi.replace("(*)", "")
#     # new_smi=new_smi.replace("*", "")
#     new_smis.append(new_smi)
# print(new_smis)

# mols=[Chem.MolFromSmiles(smi) for smi in new_smis]
# for mol in mols:
#     for i, atom in enumerate(mol.GetAtoms()):
#         atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
# fig=Draw.MolsToGridImage(mols,subImgSize=(500,500))
# fig.save("./mols_canonical.png")


# smi="CSCCC(NCc1ccc2c(c1)CN(C)C2)c1nnc2ccccn12"#"C1=CCC(C2(C34CC(C3)OC4c3nccs3)=CCNCCC2)C=C1"
# m = Chem.MolFromSmiles(smi)
# res = list(BRICS.BRICSDecompose(m))
# fragms = [Chem.MolFromSmiles(x) for x in sorted(res)]
# ms = BRICS.BRICSBuild(fragms)
# print(type(ms))
# prods = [next(ms) for x in range(20)]
# fig=Chem.Draw.MolsToGridImage(prods, molsPerRow=4, subImgSize=(200, 200))
# fig.save("./BRICS_BUILD.png")



