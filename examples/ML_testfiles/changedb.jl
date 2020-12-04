using PyCall

py"""
import ase
path="H2_Cu.db"
db2=ase.db.connect((path))
metadata={}
metadata["ref"]="H2 on different cupper facets and cupper (bulk), data from ref. Zhu et al., Phys.Chem.Chem.Phys.,2020, 22, 13958"
metadata["force_mask"]=True
metadata["force_mask_index_1"]=[1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0]
metadata["force_mask_index_2"]=[1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0]
metadata["force_mask_string_1"]=["H","H","Cu","Cu","Cu","Cu","Cu","Cu","Cu","Cu","Cu","Cu","Cu","Cu","Cu","Cu","Cu","Cu"]
metadata["force_mask_string_2"]=["Cu","Cu","Cu","Cu","Cu","Cu","Cu","Cu","Cu","Cu","Cu","Cu","Cu","Cu","Cu","Cu"]
metadata["force_mask_string_3"]=["H","H"]
metadata["force_mask_index_3"]= [1,1]
print(metadata)
db2.metadata=metadata
"""