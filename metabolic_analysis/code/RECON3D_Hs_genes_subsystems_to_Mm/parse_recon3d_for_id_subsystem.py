import re
import json


infile=open("Recon3D.json")
out=open("Recon3D_id_subsystem.txt","w")
inf=infile.read()

r3d=json.loads(inf)

#Example field input:
#{'id': '25HVITD3t', 'name': '25-hydroxyvitamin D3 transport from cytoplasm', 'metabolites': {'25hvitd3_c': -1.0, '25hvitd3_e': 1.0}, 'lower_bound': 0.0, 'upper_bound': 1000.0, 'gene_reaction_rule': '', 'subsystem': 'Transport, extracellular', 'notes': {'original_bigg_ids': ['25HVITD3t']}, 'annotation': {'bigg.reaction': ['25HVITD3t'], 'metanetx.reaction': ['MNXR94739'], 'sbo': 'SBO:0000185'}}
#>>> r3d['reactions'][1]['subsystem']
#'Transport, extracellular'

r3drx=r3d['reactions']
nlen=len(r3drx)

for rxn in range(0,nlen):
   rx=r3drx[rxn]

   print(rx['id'], rx['name'], rx['subsystem'],'\n')
   out.write(rx['id'] + "\t" +  rx['name'] + "\t" + rx['subsystem'] + "\t" + '\n')

out.close()
infile.close()
