
# coding: utf-8

# # Create a hetnet of genesets for automatic gene expression compression interpretation
# 
# This script was modified from https://github.com/dhimmel/integrate.
# 
# The script creates a hetnet as described in the eLIFE publication _"Systematic integration of biomedical knowledge prioritizes drugs for repurposing"_ by [Himmelstein et al. 2017](https://doi.org/10.7554/eLife.26726)
# 
# ## Datasets
# 
# 1. [MSigDb](https://doi.org/10.1073/pnas.0506580102 "Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles") - Curated genesets that represent various biological processes
# 2. [xCell](https://doi.org/10.1186/s13059-017-1349-1, "xCell: digitally portraying the tissue cellular heterogeneity landscape") - Curated genesets that describe profiles of different cell-types

# In[1]:


import os
import csv
import pandas as pd
import seaborn as sns

import hetio.hetnet
import hetio.readwrite
import hetio.stats


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# ## Define the metagraph and instantiate the graph

# In[3]:


kind_to_abbev = {
    
    # metanodes
    'Gene': 'G',

    # MSigDB Nodes
    'Cancer-Hallmarks': 'H',
    'Positional-Gene-Sets': 'C1',
    'Curated-Gene-Sets-CPG': 'C2CPG',
    'Curated-Gene-Sets-REACTOME': 'C2CPREACTOME',
    'Motif-Gene-Sets-MIR': 'C3MIR',
    'Motif-Gene-Sets-TFT': 'C3TFT',
    'Computational-Gene-Sets-CGN': 'C4CGN',
    'Computational-Gene-Sets-CM': 'C4CM',
    'GO-Gene-Sets-BP': 'C5BP',
    'GO-Gene-Sets-CC': 'C5CC',
    'GO-Gene-Sets-MF': 'C5MF',
    'Oncogenic-Gene-Sets': 'C6',
    'Immunologic-Gene-Sets': 'C7',
    
    # xCell Nodes
    'xCell-Cell-Type': 'XCELL',
    
    # metaedges
    'participates': 'p',
}

metaedge_tuples = [
    # MSigDB metaedges
    ('Gene', 'Cancer-Hallmarks', 'participates', 'both'),
    ('Gene', 'Positional-Gene-Sets', 'participates', 'both'),
    ('Gene', 'Curated-Gene-Sets-CPG', 'participates', 'both'),
    ('Gene', 'Curated-Gene-Sets-REACTOME', 'participates', 'both'),
    ('Gene', 'Motif-Gene-Sets-MIR', 'participates', 'both'),
    ('Gene', 'Motif-Gene-Sets-TFT', 'participates', 'both'),
    ('Gene', 'Computational-Gene-Sets-CGN', 'participates', 'both'),
    ('Gene', 'Computational-Gene-Sets-CM', 'participates', 'both'),
    ('Gene', 'GO-Gene-Sets-BP', 'participates', 'both'),
    ('Gene', 'GO-Gene-Sets-CC', 'participates', 'both'),
    ('Gene', 'GO-Gene-Sets-MF', 'participates', 'both'),
    ('Gene', 'Oncogenic-Gene-Sets', 'participates', 'both'),
    ('Gene', 'Immunologic-Gene-Sets', 'participates', 'both'),
    
    # xCell metaedges
    ('Gene', 'xCell-Cell-Type', 'participates', 'both'),
]


# In[4]:


# Initialize the graph
metagraph = hetio.hetnet.MetaGraph.from_edge_tuples(metaedge_tuples, kind_to_abbev)
graph = hetio.hetnet.Graph(metagraph)


# ## Gene Nodes

# In[5]:


# Load curated gene names from versioned resource 
commit = '721204091a96e55de6dcad165d6d8265e67e2a48'
url = 'https://raw.githubusercontent.com/cognoma/genes/{}/data/genes.tsv'.format(commit)
gene_df = pd.read_table(url)

# Only consider protein-coding genes
gene_df = (
    gene_df.query("gene_type == 'protein-coding'")
)

coding_genes = set(gene_df['entrez_gene_id'].astype(int))
print(gene_df.shape)
gene_df.head(2)


# ## Add genes as nodes to the graph
# 
# Use the gene-symbol identifier for easier interpretation

# In[6]:


get_ipython().run_cell_magic('time', '', "for i, row in gene_df.iterrows():\n    # Build dictionary of descriptive elements for each gene\n    meta_data = {\n        'description': row['description'],\n        'source': 'Entrez Gene',\n        'url': 'http://identifiers.org/ncbigene/{}'.format(row['entrez_gene_id']),\n        'license': 'CC0 1.0',\n    }\n    \n    if pd.notnull(row['chromosome']):\n        meta_data['chromosome'] = row['chromosome']\n\n    # Add genes to graph\n    graph.add_node(kind='Gene', identifier=int(row['entrez_gene_id']), name=row['symbol'],\n                   data=meta_data)")


# In[7]:


# Load gene updater
url = 'https://raw.githubusercontent.com/cognoma/genes/{}/data/updater.tsv'.format(commit)
updater_df = pd.read_table(url)
old_to_new_entrez = dict(zip(updater_df.old_entrez_gene_id,
                             updater_df.new_entrez_gene_id))


# ## Add gene set nodes and associated genes as edges
# 
# Add each MSigDB collection as distinct nodes with a `participates` edge for representative gene sets and corresponding membership.

# In[8]:


def add_node_to_graph(current_graph, collection_file, collection_kind,
                      collection_source, gene_list, min_geneset_size=4,
                      max_geneset_size=1000, license='CC BY 4.0'):
    """
    Add nodes and edges to current graph based on geneset memembership of collection
    
    Arguments:
    current_graph - a hetnet object to add node-edge info to
    collection_file - location of msigdb file
    collection_kind - the kind of node already initialized in the graph
    collection_source - alternative ID for collection
    gene_list - a list of genes to consider when building the graph
    min_geneset_size - filter out a given gene set if it has fewer genes
    max_geneset_size - filter out a given gene set if it has more genes
    license - given license associated with node
    
    Output:
    Adds to current graph; Returns the amount of filtered genesets
    """
    
    # Build meta data dictionary to store node info
    meta_data = {'license': license, 'source': collection_source}
    
    # Open the .gmt file and process each geneset
    filtered_genesets = []
    with open(collection_file, 'r') as collection_fh:
        collection_reader = csv.reader(collection_fh, delimiter='\t')

        for row in collection_reader:
            # Get geneset and and metadata info
            geneset_name = row[0]
            meta_data['url'] = row[1]
            
            # Process geneset membership
            genes = row[2:]

            # Update entrez_gene_id
            genes = set(old_to_new_entrez[x] if x in old_to_new_entrez else x for x in genes)
            
            # The genes must exist in curated resource
            genes = [int(x) for x in genes if int(x) in gene_list]

            # Filter geneset if its too big or small
            if min_geneset_size > len(genes) or len(genes) > max_geneset_size:
                filtered_genesets.append(geneset_name)
                continue
                
            # Add the genesetname as a node (based on collection) to the graph
            current_graph.add_node(kind=collection_kind,
                                   identifier=geneset_name,
                                   data=meta_data)
            # Loop through all genes and add to the graph it should be considered
            for gene in genes:
                source_id = ('Gene', gene)
                target_id = (collection_kind, geneset_name)
                edge_data = meta_data.copy()
                current_graph.add_edge(source_id, target_id, 'participates',
                                       'both', edge_data)

    return filtered_genesets


# In[9]:


hetnet_build = {
    # Format: `Collection Source`: [`Collection File`, `Collection Kind`]

    # MSigDB
    'MSigDB-H': ['h.all.v6.1.entrez.gmt', 'Cancer-Hallmarks'],
    'MSigDB-C1': ['c1.all.v6.1.entrez.gmt', 'Positional-Gene-Sets'],
    'MSigDB-C2-CPG': ['c2.cgp.v6.1.entrez.gmt', 'Curated-Gene-Sets-CPG'],
    'MSigDB-C2-Reactome': ['c2.cp.reactome.v6.1.entrez.gmt', 'Curated-Gene-Sets-REACTOME'],
    'MSigDB-C3-MIR': ['c3.mir.v6.1.entrez.gmt', 'Motif-Gene-Sets-MIR'],
    'MSigDB-C3-TFT': ['c3.tft.v6.1.entrez.gmt', 'Motif-Gene-Sets-TFT'],
    'MSigDB-C4-CGN': ['c4.cgn.v6.1.entrez.gmt', 'Computational-Gene-Sets-CGN'],
    'MSigDB-C4-CM': ['c4.cm.v6.1.entrez.gmt', 'Computational-Gene-Sets-CM'],
    'MSigDB-C5-BP': ['c5.bp.v6.1.entrez.gmt', 'GO-Gene-Sets-BP'],
    'MSigDB-C5-CC': ['c5.cc.v6.1.entrez.gmt', 'GO-Gene-Sets-CC'],
    'MSigDB-C5-MF': ['c5.mf.v6.1.entrez.gmt', 'GO-Gene-Sets-MF'],
    'MSigDB-C6': ['c6.all.v6.1.entrez.gmt', 'Oncogenic-Gene-Sets'],
    'MSigDB-C7': ['c7.all.v6.1.entrez.gmt', 'Immunologic-Gene-Sets'],
    
    # xCell
    'xCell-X': ['xcell_all_entrez.gmt', 'xCell-Cell-Type'],

}


# In[10]:


get_ipython().run_cell_magic('time', '', "\n# Add all collections genesets to hetnet\nfiltered = {}\nfor collection_source, collection_info in hetnet_build.items():\n    path, collection_kind = collection_info\n    collection_file = os.path.join('data', path)\n    filtered[collection_kind] = add_node_to_graph(current_graph=graph,\n                                                  collection_file=collection_file,\n                                                  collection_kind=collection_kind,\n                                                  collection_source=collection_source,\n                                                  gene_list=coding_genes)")


# ## Network visualizations and stats

# In[11]:


# Export node degree tables
node_degree_file = os.path.join('results', 'interpret_node_degrees.xlsx')
hetio.stats.degrees_to_excel(graph, node_degree_file)


# In[12]:


# Summary of metanodes and cooresponding nodes
metanode_df = hetio.stats.get_metanode_df(graph)

metanode_file = os.path.join('results', 'interpret_metanode_summary.tsv')
metanode_df.to_csv(metanode_file, sep='\t', index=False)
metanode_df


# In[13]:


# Summary of metaedges and cooresponding edges
metaedge_df = hetio.stats.get_metaedge_df(graph)

rows = list()
for metaedge, edges in graph.get_metaedge_to_edges(exclude_inverts=True).items():
    rows.append({'metaedge': str(metaedge)})

metaedge_file = os.path.join('results', 'interpret_metaedges.tsv')
metaedge_df = metaedge_df.merge(pd.DataFrame(rows))

sum_total = metaedge_df.sum()
sum_total.metaedge = 'Total'
sum_total.abbreviation = ''

metaedge_df = (
    pd.concat([metaedge_df.T, sum_total], axis='columns')
    .transpose()
    .reset_index(drop=True)
)
# Number of edges in the network
metaedge_df.edges.sum()

metaedge_df.to_csv(metaedge_file, sep='\t', index=False)
metaedge_df


# In[14]:


# Summary of different styles for representing each metaedge
metaedge_style_file = os.path.join('results', 'interpret_metaedge_styles.tsv')
metaedge_style_df = hetio.stats.get_metaedge_style_df(metagraph)
metaedge_style_df.to_csv(metaedge_style_file, sep='\t', index=False)
metaedge_style_df


# In[15]:


# How many genesets were filtered per collection?
{x: len(y) for x, y in filtered.items()}


# ## Save graph

# In[16]:


get_ipython().run_cell_magic('time', '', "# Write nodes to a table\nnodes_file = os.path.join('hetnets', 'interpret_nodes.tsv')\nhetio.readwrite.write_nodetable(graph, nodes_file)\n\n# Write edges to a table\nedges_file = os.path.join('hetnets', 'interpret_edges.sif.gz')\nhetio.readwrite.write_sif(graph, edges_file)")


# In[17]:


get_ipython().run_cell_magic('time', '', "# Write metagraph as json\nmetagraph_file = os.path.join('hetnets', 'interpret_metagraph.json')\nhetio.readwrite.write_metagraph(metagraph, metagraph_file)")


# In[18]:


get_ipython().run_cell_magic('time', '', "# Write graph as json\nhetnet_json_path = os.path.join('hetnets', 'interpret_hetnet.json.bz2')\nhetio.readwrite.write_graph(graph, hetnet_json_path)")


# In[19]:


get_ipython().system(" sha256sum 'hetnets/interpret_hetnet.json.bz2'")


# ## Visualize hetnet node and edge counts

# In[20]:


ax = sns.barplot(x='metanode', y='nodes', data=metanode_df.sort_values('nodes'))
for tick in ax.get_xticklabels():
    tick.set_rotation(90)
ax.set_xlabel(''); ax.set_ylabel('nodes');


# In[21]:


ax = sns.barplot(x='metaedge', y='edges', data=metaedge_df.sort_values('edges'))
for tick in ax.get_xticklabels():
    tick.set_rotation(90)
ax.set_xlabel(''); ax.set_ylabel('edges');

