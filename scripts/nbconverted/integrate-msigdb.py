
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
    'Curated-Gene-Sets': 'C2',
    'Motif-Gene-Sets': 'C3',
    'Computational-Gene-Sets': 'C4',
    'GO-Gene-Sets': 'C5',
    'Oncogenic-Gene-Sets': 'C6',
    'Immunologic-Gene-Sets': 'C7',
    
    # metaedges
    'participates': 'p',
}

metaedge_tuples = [
    ('Gene', 'Cancer-Hallmarks', 'participates', 'both'),
    ('Gene', 'Positional-Gene-Sets', 'participates', 'both'),
    ('Gene', 'Curated-Gene-Sets', 'participates', 'both'),
    ('Gene', 'Motif-Gene-Sets', 'participates', 'both'),
    ('Gene', 'Computational-Gene-Sets', 'participates', 'both'),
    ('Gene', 'GO-Gene-Sets', 'participates', 'both'),
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
commit = 'a7362748a34211e5df6f2d185bb3246279760546'
url = 'https://raw.githubusercontent.com/dhimmel/entrez-gene/{}/data/genes-human.tsv'.format(commit)
gene_df = pd.read_table(url)

# Only consider protein-coding genes
gene_df = (
    gene_df.query("type_of_gene == 'protein-coding'")
    .drop_duplicates('Symbol')
)

coding_genes = set(gene_df['GeneID'])
print(gene_df.shape)
gene_df.head(2)


# In[9]:


# Subset the msigdb genes to the gene curation above
common_subset_genes = set(gene_df['Symbol']).intersection(set(msigdb_df.index))
msigdb_subset_genes_df = msigdb_df.loc[common_subset_genes, :]

# What is the new distribution of pathway sizes?
msigdb_subset_genes_df.sum(axis=1).sort_values().hist(bins=30);


# In[10]:


# What is the new distribution of gene representation?
msigdb_subset_genes_df.sum(axis=0).sort_values().hist(bins=30);


# In[11]:


# The genes that were removed from MSigDB participate mostly in few pathways
diff_genes = set(msigdb_df.index).difference(gene_df['Symbol'])
msigdb_other_genes_df = msigdb_df.loc[diff_genes, :]
msigdb_other_genes_df.sum(axis=1).sort_values().hist(bins=30)


# ## Add genes as nodes to the graph
# 
# Use the gene-symbol identifier for easier interpretation

# In[12]:


get_ipython().run_cell_magic('time', '', "for i, row in gene_df.iterrows():\n    # Build dictionary of descriptive elements for each gene\n    meta_data = {\n        'description': row['description'],\n        'source': 'Entrez Gene',\n        'url': 'http://identifiers.org/ncbigene/{}'.format(row['GeneID']),\n        'license': 'CC0 1.0',\n    }\n    \n    if pd.notnull(row['chromosome']):\n        meta_data['chromosome'] = row['chromosome']\n\n    # Add genes to graph\n    graph.add_node(kind='Gene', identifier=row['Symbol'], name=row['GeneID'],\n                   data=meta_data)")


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
            current_graph.add_node(kind=collection_kind, identifier=geneset_name, data=meta_data)

            # Loop through all genes and add to the graph it should be considered
            for gene in genes:
                if gene in gene_list:
                    source_id = ('Gene', gene)
                    target_id = (collection_kind, geneset_name)
                    
                    edge_data = meta_data.copy()
                    edge_data['unbiased'] = False

                    current_graph.add_edge(source_id, target_id, 'participates', 'both', edge_data)

    return filtered_genesets


# In[14]:


get_ipython().run_cell_magic('time', '', "hallmarks_file = 'data/h.all.v6.1.symbols.gmt'\nhallmark_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                             collection_file=hallmarks_file,\n                                             collection_kind='Cancer-Hallmarks',\n                                             collection_source='MSigDB-H',\n                                             gene_list=common_subset_genes)")


# In[15]:


get_ipython().run_cell_magic('time', '', "positional_file = 'data/c1.all.v6.1.symbols.gmt'\npositional_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                               collection_file=positional_file,\n                                               collection_kind='Positional-Gene-Sets',\n                                               collection_source='MSigDB-C1',\n                                               gene_list=common_subset_genes)")


# In[16]:


get_ipython().run_cell_magic('time', '', "curated_file = 'data/c2.all.v6.1.symbols.gmt'\ncurated_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                            collection_file=curated_file,\n                                            collection_kind='Curated-Gene-Sets',\n                                            collection_source='MSigDB-C2',\n                                            gene_list=common_subset_genes)")


# In[17]:


get_ipython().run_cell_magic('time', '', "motif_file = 'data/c3.all.v6.1.symbols.gmt'\nmotif_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                          collection_file=motif_file,\n                                          collection_kind='Motif-Gene-Sets',\n                                          collection_source='MSigDB-C3',\n                                          gene_list=common_subset_genes)")


# In[18]:


get_ipython().run_cell_magic('time', '', "computational_file = 'data/c4.all.v6.1.symbols.gmt'\ncomputational_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                                  collection_file=computational_file,\n                                                  collection_kind='Computational-Gene-Sets',\n                                                  collection_source='MSigDB-C4',\n                                                  gene_list=common_subset_genes)")


# In[19]:


get_ipython().run_cell_magic('time', '', "go_file = 'data/c5.all.v6.1.symbols.gmt'\ngo_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                       collection_file=go_file,\n                                       collection_kind='GO-Gene-Sets',\n                                       collection_source='MSigDB-C5',\n                                       gene_list=common_subset_genes) ")


# In[20]:


get_ipython().run_cell_magic('time', '', "oncogenic_file = 'data/c6.all.v6.1.symbols.gmt'\noncogenic_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                              collection_file=oncogenic_file,\n                                              collection_kind='Oncogenic-Gene-Sets',\n                                              collection_source='MSigDB-C6',\n                                              gene_list=common_subset_genes)")


# In[21]:


get_ipython().run_cell_magic('time', '', "immunologic_file = 'data/c7.all.v6.1.symbols.gmt'\nimmunologic_filtered = add_msigdb_node_to_graph(current_graph=graph,\n                                                collection_file=immunologic_file,\n                                                collection_kind='Immunologic-Gene-Sets',\n                                                collection_source='MSigDB-C7',\n                                                gene_list=common_subset_genes)")


# ## Network visualizations and stats

# In[22]:


# Export node degree tables
node_degree_file = os.path.join('results', 'msigdb_node_degrees.xlsx')
hetio.stats.degrees_to_excel(graph, node_degree_file)


# In[23]:


# Summary of metanodes and cooresponding nodes
metanode_df = hetio.stats.get_metanode_df(graph)

metanode_file = os.path.join('results', 'msigdb_metanode_summary.tsv')
metanode_df.to_csv(metanode_file, sep='\t', index=False)
metanode_df


# In[24]:


# Total number of nodes
metanode_df.nodes.sum()


# In[25]:


# Summary of metaedges and cooresponding edges
metaedge_df = hetio.stats.get_metaedge_df(graph)

# Calculate number of unbiased edges
rows = list()
for metaedge, edges in graph.get_metaedge_to_edges(exclude_inverts=True).items():
    unbiased = sum(edge.data['unbiased'] for edge in edges)
    rows.append({'metaedge': str(metaedge), 'unbiased': unbiased})

metaedge_file = os.path.join('results', 'msigdb_metaedges.tsv')
metaedge_df = metaedge_df.merge(pd.DataFrame(rows))
metaedge_df.to_csv(metaedge_file, sep='\t', index=False)
metaedge_df


# In[26]:


# Summary of different styles for representing each metaedge
metaedge_style_file = os.path.join('results', 'msigdb_metaedge_styles.tsv')
metaedge_style_df = hetio.stats.get_metaedge_style_df(metagraph)
metaedge_style_df.to_csv(metaedge_style_file, sep='\t', index=False)
metaedge_style_df


# In[27]:


# Number of edges in the network
metaedge_df.edges.sum()


# In[28]:


# How many genesets were filtered?
{'Cancer-Hallmarks': len(hallmark_filtered),
 'Positional-Gene-Sets': len(positional_filtered),
 'Curated-Gene-Sets': len(curated_filtered),
 'Motif-Gene-Sets': len(motif_filtered),
 'Computational-Gene-Sets': len(computational_filtered),
 'GO-Gene-Sets': len(go_filtered),
 'Oncogenic-Gene-Sets': len(oncogenic_filtered),
 'Immunologic-Gene-Sets': len(immunologic_filtered)}


# ## Save graph

# In[29]:


get_ipython().run_cell_magic('time', '', "# Write nodes to a table\nnodes_file = os.path.join('hetnets', 'msigdb_nodes.tsv')\nhetio.readwrite.write_nodetable(graph, nodes_file)\n\n# Write edges to a table\nedges_file = os.path.join('hetnets', 'msigdb_edges.sif.gz')\nhetio.readwrite.write_sif(graph, edges_file)")


# In[30]:


get_ipython().run_cell_magic('time', '', "# Write a subset of edges to a table\nedges_file = os.path.join('hetnets', 'msigdb_edges-10.sif.gz')\nhetio.readwrite.write_sif(graph, edges_file, max_edges=10)\n\nedges_file = os.path.join('hetnets', 'msigdb_edges-5k.sif.gz')\nhetio.readwrite.write_sif(graph, edges_file, max_edges=5000)")


# In[31]:


get_ipython().run_cell_magic('time', '', "# Write metagraph as json\nmetagraph_file = os.path.join('hetnets', 'msigdb_metagraph.json')\nhetio.readwrite.write_metagraph(metagraph, metagraph_file)")


# In[32]:


get_ipython().run_cell_magic('time', '', "# Write graph as json\nhetnet_json_path = os.path.join('hetnets', 'msigdb_hetnet.json.bz2')\nhetio.readwrite.write_graph(graph, hetnet_json_path)")


# In[33]:


get_ipython().system(" sha256sum 'hetnets/msigdb_hetnet.json.bz2'")


# ## Visualize hetnet node and edge counts

# In[34]:


ax = sns.barplot(x='metanode', y='nodes', data=metanode_df.sort_values('nodes'))
for tick in ax.get_xticklabels():
    tick.set_rotation(90)
ax.set_xlabel(''); ax.set_ylabel('nodes');


# In[35]:


ax = sns.barplot(x='metaedge', y='edges', data=metaedge_df.sort_values('edges'))
for tick in ax.get_xticklabels():
    tick.set_rotation(90)
ax.set_xlabel(''); ax.set_ylabel('edges');

