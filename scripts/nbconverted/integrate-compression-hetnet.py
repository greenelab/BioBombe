
# coding: utf-8

# # Create a hetnet of MSigDB genesets
# 
# This script was modified from https://github.com/dhimmel/integrate.
# 
# ## For the automatic interpretation of gene expression compression features
# 
# The script creates an MSigDB hetnet as described in the eLIFE publication _"Systematic integration of biomedical knowledge prioritizes drugs for repurposing"_ by [Himmelstein et al. 2017](https://doi.org/10.7554/eLife.26726)

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

    'Cancer-Hallmarks': 'H',
    'Positional-Gene-Sets': 'C1',
    'Curated-Gene-Sets-CPG': 'C2-CPG',
    'Curated-Gene-Sets-REACTOME': 'C2-CP-REACTOME',
    'Motif-Gene-Sets-MIR': 'C3-MIR',
    'Motif-Gene-Sets-TFT': 'C3-TFT',
    'Computational-Gene-Sets-CGN': 'C4-CGN',
    'Computational-Gene-Sets-CM': 'C4-CM',
    'GO-Gene-Sets-BP': 'C5-BP',
    'GO-Gene-Sets-CC': 'C5-CC',
    'GO-Gene-Sets-MF': 'C5-MF',
    'Oncogenic-Gene-Sets': 'C6',
    'Immunologic-Gene-Sets': 'C7',
    
    # metaedges
    'participates': 'p',
}

metaedge_tuples = [
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
]


# In[4]:


# Initialize the graph
metagraph = hetio.hetnet.MetaGraph.from_edge_tuples(metaedge_tuples, kind_to_abbev)
graph = hetio.hetnet.Graph(metagraph)


# ## Get all genes found in MSigDB

# In[5]:


full_msigdb_file = os.path.join('data', 'full_msigdb_binary_matrix.tsv.bz2')
msigdb_df = pd.read_table(full_msigdb_file, index_col=0)

print(msigdb_df.shape)
msigdb_df.head()


# In[6]:


# Ask the distribution of pathway sizes
msigdb_df.sum(axis=1).sort_values().hist(bins=30);


# In[7]:


# What about the distribution of gene participation?
msigdb_df.sum(axis=0).sort_values().hist(bins=30);


# ## Gene Nodes

# In[8]:


# Load curated gene names from versioned resource 
commit = '721204091a96e55de6dcad165d6d8265e67e2a48'
url = 'https://raw.githubusercontent.com/cognoma/genes/{}/data/genes.tsv'.format(commit)
gene_df = pd.read_table(url)

# Only consider protein-coding genes
gene_df = (
    gene_df.query("gene_type == 'protein-coding'")
    .drop_duplicates('symbol')
)

coding_genes = set(gene_df['symbol'])
print(gene_df.shape)
gene_df.head(2)


# In[9]:


# Subset the msigdb genes to the gene curation above
common_subset_genes = set(gene_df['symbol']).intersection(set(msigdb_df.index))
msigdb_subset_genes_df = msigdb_df.loc[common_subset_genes, :]

# What is the new distribution of pathway sizes?
msigdb_subset_genes_df.sum(axis=1).sort_values().hist(bins=30);


# In[10]:


# What is the new distribution of gene representation?
msigdb_subset_genes_df.sum(axis=0).sort_values().hist(bins=30);


# In[11]:


# The genes that were removed from MSigDB participate mostly in few pathways
diff_genes = set(msigdb_df.index).difference(gene_df['symbol'])
msigdb_other_genes_df = msigdb_df.loc[diff_genes, :]
msigdb_other_genes_df.sum(axis=1).sort_values().hist(bins=30);


# ## Add genes as nodes to the graph
# 
# Use the gene-symbol identifier for easier interpretation

# In[12]:


get_ipython().run_cell_magic('time', '', "for i, row in gene_df.iterrows():\n    # Build dictionary of descriptive elements for each gene\n    meta_data = {\n        'description': row['description'],\n        'source': 'Entrez Gene',\n        'url': 'http://identifiers.org/ncbigene/{}'.format(row['entrez_gene_id']),\n        'license': 'CC0 1.0',\n    }\n    \n    if pd.notnull(row['chromosome']):\n        meta_data['chromosome'] = row['chromosome']\n\n    # Add genes to graph\n    graph.add_node(kind='Gene', identifier=row['symbol'], name=row['entrez_gene_id'],\n                   data=meta_data)")


# ## Add gene set nodes and associated genes as edges
# 
# Add each MSigDB collection as distinct nodes with a `participates` edge for representative gene sets and corresponding membership.

# In[13]:


def add_msigdb_node_to_graph(current_graph, collection_file, collection_kind,
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
                if gene in gene_list:
                    source_id = ('Gene', gene)
                    target_id = (collection_kind, geneset_name)
                    
                    edge_data = meta_data.copy()

                    current_graph.add_edge(source_id, target_id, 'participates',
                                           'both', edge_data)

    return filtered_genesets


# In[14]:


get_ipython().run_cell_magic('time', '', "hallmarks_file = os.path.join('data', 'h.all.v6.1.symbols.gmt')\nhallmark_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                             collection_file=hallmarks_file,\n                                             collection_kind='Cancer-Hallmarks',\n                                             collection_source='MSigDB-H',\n                                             gene_list=common_subset_genes)")


# In[15]:


get_ipython().run_cell_magic('time', '', "positional_file = os.path.join('data', 'c1.all.v6.1.symbols.gmt')\npositional_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                               collection_file=positional_file,\n                                               collection_kind='Positional-Gene-Sets',\n                                               collection_source='MSigDB-C1',\n                                               gene_list=common_subset_genes)")


# In[16]:


get_ipython().run_cell_magic('time', '', "curated_file = os.path.join('data', 'c2.cgp.v6.1.symbols.gmt')\ncurated_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                            collection_file=curated_file,\n                                            collection_kind='Curated-Gene-Sets-CPG',\n                                            collection_source='MSigDB-C2-CPG',\n                                            gene_list=common_subset_genes)")


# In[17]:


get_ipython().run_cell_magic('time', '', "reactome_file = os.path.join('data', 'c2.cp.reactome.v6.1.symbols.gmt')\nreactome_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                             collection_file=reactome_file,\n                                             collection_kind='Curated-Gene-Sets-REACTOME',\n                                             collection_source='MSigDB-C2-Reactome',\n                                             gene_list=common_subset_genes)")


# In[18]:


get_ipython().run_cell_magic('time', '', "micro_file = os.path.join('data', 'c3.mir.v6.1.symbols.gmt')\nmicro_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                          collection_file=micro_file,\n                                          collection_kind='Motif-Gene-Sets-MIR',\n                                          collection_source='MSigDB-C3-MIR',\n                                          gene_list=common_subset_genes)")


# In[19]:


get_ipython().run_cell_magic('time', '', "tf_file = os.path.join('data', 'c3.tft.v6.1.symbols.gmt')\ntf_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                          collection_file=tf_file,\n                                          collection_kind='Motif-Gene-Sets-TFT',\n                                          collection_source='MSigDB-C3-TFT',\n                                          gene_list=common_subset_genes)")


# In[20]:


get_ipython().run_cell_magic('time', '', "cancer_n_file = os.path.join('data', 'c4.cgn.v6.1.symbols.gmt')\ncancer_n_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                             collection_file=cancer_n_file,\n                                             collection_kind='Computational-Gene-Sets-CGN',\n                                             collection_source='MSigDB-C4-CGN',\n                                             gene_list=common_subset_genes)")


# In[21]:


get_ipython().run_cell_magic('time', '', "cancer_mod_file = os.path.join('data', 'c4.cm.v6.1.symbols.gmt')\ncancer_mod_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                             collection_file=cancer_mod_file,\n                                             collection_kind='Computational-Gene-Sets-CM',\n                                             collection_source='MSigDB-C4-CM',\n                                             gene_list=common_subset_genes)")


# In[22]:


get_ipython().run_cell_magic('time', '', "go_bp_file = os.path.join('data', 'c5.bp.v6.1.symbols.gmt')\ngo_bp_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                          collection_file=go_bp_file,\n                                          collection_kind='GO-Gene-Sets-BP',\n                                          collection_source='MSigDB-C5-BP',\n                                          gene_list=common_subset_genes) ")


# In[23]:


get_ipython().run_cell_magic('time', '', "go_cc_file = os.path.join('data', 'c5.cc.v6.1.symbols.gmt')\ngo_cc_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                          collection_file=go_cc_file,\n                                          collection_kind='GO-Gene-Sets-CC',\n                                          collection_source='MSigDB-C5-CC',\n                                          gene_list=common_subset_genes) ")


# In[24]:


get_ipython().run_cell_magic('time', '', "go_mf_file = os.path.join('data', 'c5.mf.v6.1.symbols.gmt')\ngo_mf_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                          collection_file=go_mf_file,\n                                          collection_kind='GO-Gene-Sets-MF',\n                                          collection_source='MSigDB-C5-MF',\n                                          gene_list=common_subset_genes) ")


# In[25]:


get_ipython().run_cell_magic('time', '', "oncogenic_file = os.path.join('data', 'c6.all.v6.1.symbols.gmt')\noncogenic_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                              collection_file=oncogenic_file,\n                                              collection_kind='Oncogenic-Gene-Sets',\n                                              collection_source='MSigDB-C6',\n                                              gene_list=common_subset_genes)")


# In[26]:


get_ipython().run_cell_magic('time', '', "immunologic_file = os.path.join('data', 'c7.all.v6.1.symbols.gmt')\nimmunologic_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                                collection_file=immunologic_file,\n                                                collection_kind='Immunologic-Gene-Sets',\n                                                collection_source='MSigDB-C7',\n                                                gene_list=common_subset_genes)")


# ## Network visualizations and stats

# In[27]:


# Export node degree tables
node_degree_file = os.path.join('results', 'msigdb_node_degrees.xlsx')
hetio.stats.degrees_to_excel(graph, node_degree_file)


# In[28]:


# Summary of metanodes and cooresponding nodes
metanode_df = hetio.stats.get_metanode_df(graph)

metanode_file = os.path.join('results', 'msigdb_metanode_summary.tsv')
metanode_df.to_csv(metanode_file, sep='\t', index=False)
metanode_df


# In[29]:


# Total number of nodes
metanode_df.nodes.sum()


# In[30]:


# Summary of metaedges and cooresponding edges
metaedge_df = hetio.stats.get_metaedge_df(graph)

rows = list()
for metaedge, edges in graph.get_metaedge_to_edges(exclude_inverts=True).items():
    rows.append({'metaedge': str(metaedge)})

metaedge_file = os.path.join('results', 'msigdb_metaedges.tsv')
metaedge_df = metaedge_df.merge(pd.DataFrame(rows))
metaedge_df.to_csv(metaedge_file, sep='\t', index=False)
metaedge_df


# In[31]:


# Summary of different styles for representing each metaedge
metaedge_style_file = os.path.join('results', 'msigdb_metaedge_styles.tsv')
metaedge_style_df = hetio.stats.get_metaedge_style_df(metagraph)
metaedge_style_df.to_csv(metaedge_style_file, sep='\t', index=False)
metaedge_style_df


# In[32]:


# Number of edges in the network
metaedge_df.edges.sum()


# In[33]:


# How many genesets were filtered?
{'Cancer-Hallmarks': len(hallmark_filtered),
 'Positional-Gene-Sets': len(positional_filtered),
 'Curated-Gene-Sets-CPG': len(curated_filtered),
 'Curated-Gene-Sets-REACTOME': len(reactome_filtered),
 'Motif-Gene-Sets-MIR': len(micro_filtered),
 'Motif-Gene-Sets-TFT': len(tf_filtered),
 'Computational-Gene-Sets-CGN': len(cancer_n_filtered),
 'Computational-Gene-Sets-CM': len(cancer_mod_filtered),
 'GO-Gene-Sets-BP': len(go_bp_filtered),
 'GO-Gene-Sets-CC': len(go_cc_filtered),
 'GO-Gene-Sets-MF': len(go_mf_filtered),
 'Oncogenic-Gene-Sets': len(oncogenic_filtered),
 'Immunologic-Gene-Sets': len(immunologic_filtered)}


# ## Save graph

# In[34]:


get_ipython().run_cell_magic('time', '', "# Write nodes to a table\nnodes_file = os.path.join('hetnets', 'msigdb_nodes.tsv')\nhetio.readwrite.write_nodetable(graph, nodes_file)\n\n# Write edges to a table\nedges_file = os.path.join('hetnets', 'msigdb_edges.sif.gz')\nhetio.readwrite.write_sif(graph, edges_file)")


# In[35]:


get_ipython().run_cell_magic('time', '', "# Write metagraph as json\nmetagraph_file = os.path.join('hetnets', 'msigdb_metagraph.json')\nhetio.readwrite.write_metagraph(metagraph, metagraph_file)")


# In[36]:


get_ipython().run_cell_magic('time', '', "# Write graph as json\nhetnet_json_path = os.path.join('hetnets', 'msigdb_hetnet.json.bz2')\nhetio.readwrite.write_graph(graph, hetnet_json_path)")


# In[37]:


get_ipython().system(" sha256sum 'hetnets/msigdb_hetnet.json.bz2'")


# ## Visualize hetnet node and edge counts

# In[38]:


ax = sns.barplot(x='metanode', y='nodes', data=metanode_df.sort_values('nodes'))
for tick in ax.get_xticklabels():
    tick.set_rotation(90)
ax.set_xlabel(''); ax.set_ylabel('nodes');


# In[39]:


ax = sns.barplot(x='metaedge', y='edges', data=metaedge_df.sort_values('edges'))
for tick in ax.get_xticklabels():
    tick.set_rotation(90)
ax.set_xlabel(''); ax.set_ylabel('edges');

