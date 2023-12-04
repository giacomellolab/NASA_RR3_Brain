

infile=open("Recon3D_pathway_Mm_gene_symbol_withATP_tbl.tab", "r")
outfile=open("Recon3D_pathway_Mm_gene_symbol_withATP_tbl.gmt ", "w")

tab="\t"

pathwaydb={}

fs=infile.readline()

while 1:
    fs=infile.readline()
    if fs=="":
        break
    fsl=fs[:-1].split(tab)
    path=fsl[0]
    gene=fsl[1]
    if path in pathwaydb.keys():
        templist=pathwaydb[path]
        templist.append(gene)
        pathwaydb[path]=templist
    else:
        pathwaydb[path]=[gene]

keys=pathwaydb.keys()
keys=list(keys)
keys.sort()

for key in keys:

    genelist=pathwaydb[key]
    outfile.write(key + tab + key + tab + tab.join(genelist) + "\n")




