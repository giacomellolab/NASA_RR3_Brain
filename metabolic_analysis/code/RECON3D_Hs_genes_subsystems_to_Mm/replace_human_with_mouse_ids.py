
infile=open("Recon3D_pathway_gene_symbol_tbl.csv", "r")
hcopfile=open("HCOP_mouse_human_15col.txt", "r")
outfile=open("Recon3D_pathway_Mm_gene_symbol_tbl.tab", "w") 

tab="\t"
hcopdb={}

fs=hcopfile.readline()
while 1:
    fs=hcopfile.readline()
    if fs=="":
        break
    fsl=fs.split(tab)
    hs=fsl[4]
    mm=fsl[11]
    if hs !="-" and mm !="-":
       try:
           templs=hcopdb[hs]
           templs.append(mm)
           hcopdb[hs]=templs
       except KeyError:
           hcopdb[hs]=[mm]

hcopfile.close()

# ,pathway,gene_symbol
#0,Alanine and aspartate metabolism,ACY3
#1,Alanine and aspartate metabolism,AGXT
geneseen=[]

fs=infile.readline()
while 1:
    fs=infile.readline()
    if fs=="":
        break
    fsl=fs[:-1].split(",")
    path=fsl[1]
    hsgene=fsl[2]
    try:
        mmgene=hcopdb[hsgene]
    except KeyError:
        mmgene=[]
    for mmg in mmgene:
        if mmg not in geneseen:

            outfile.write(path + tab + mmg  + "\n")
            geneseen.append(mmg)





