import cobra
import itertools as it
import pandas as pd

def get_bigg_mmu_nopt_model():
    bigg_mmu_json = cobra.io.load_json_model('Recon3D.json')

    bigg_rxn_subsys_lut = {}
    for ir in bigg_mmu_json.reactions:
        if ir._id in bigg_rxn_subsys_lut:
            print(
                'dup {} {} {}'.format(
                    ir._id, ir._subsystem, bigg_rxn_subsys_lut[ir._id]))
        bigg_rxn_subsys_lut[ir._id] = ir.subsystem

    bigg_mmu_xml = cobra.io.read_sbml_model('Recon3D.xml')

    for ir in bigg_mmu_xml.reactions:
        ir.subsystem = bigg_rxn_subsys_lut[ir._id]

    return bigg_mmu_xml

Recon3D_model = get_bigg_mmu_nopt_model()

entrez_id_symbol_lut = {}
#entrez_id_symbol_lut['514'] = 'Atp5e' # hu
#entrez_id_symbol_lut['10476'] = 'Atp5h' # hu
#entrez_id_symbol_lut['4710'] = 'Ndufb4' # hu
#entrez_id_symbol_lut['4541'] = 'mt-Nd6' #hu
#entrez_id_symbol_lut['4707'] = 'Ndufb1-ps' #hu
#entrez_id_symbol_lut['222'] = 'Aldh3b2' #hu
#entrez_id_symbol_lut['26237'] = 'Eno1b' # hu
#entrez_id_symbol_lut['348477'] = 'Mia3' # hu


all_subsystems = sorted(set(x.subsystem for x in Recon3D_model.reactions))

pathway_gene_list = []
for i_pathway in all_subsystems:
    i_rxns = [
        x for x in Recon3D_model.reactions
        if x.subsystem == i_pathway]

    i_genes = sorted(set(list(
        it.chain(*[[xx.name for xx in x.genes] for x in i_rxns]))))

    for j in i_genes:
        if j in entrez_id_symbol_lut:
            pathway_gene_list.append((i_pathway, entrez_id_symbol_lut[j]))
        else:
            pathway_gene_list.append((i_pathway, j))

del i_pathway, i_rxns, i_genes, j

pathway_gene_df = pd.DataFrame(
    pathway_gene_list, columns = ['pathway', 'gene_symbol'])

pathway_gene_df.to_csv('Recon3D_pathway_gene_symbol_tbl.csv')
