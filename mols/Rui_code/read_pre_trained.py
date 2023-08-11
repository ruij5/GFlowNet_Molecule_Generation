import pandas as pd

import json
import pickle



# with open("mols/data/Pre_trained_data/info2.pkl", 'rb') as file:
#     object = pickle.load(file)
# #print(object.keys())


# with open("mols/data/Pre_trained_data/best_params.pkl", 'rb') as file2:
#     object2 = pickle.load(file2)
# #print(object2)


# with open("mols/data/Pre_trained_data/params.pkl", 'rb') as file3:
#     object3 = pickle.load(file3)
# #print(object3)

blocks = pd.read_json('mols/data/blocks_PDB_105.json')
block_rs = blocks["block_r"].to_list()
print(block_rs[73][0])