
# coding: utf-8

# # Permute MSigDB Hetnets
# 
# Modified from @dhimmel - https://github.com/dhimmel/integrate/blob/master/permute.ipynb
# 
# Generate several randomly permuted hetnets to serve as a null distribution. The permutations preserve node degree but randomizes connections between nodes. See [Himmelstein et al. 2017](https://doi.org/10.7554/eLife.26726) for more details.

# In[1]:


import os
import pandas as pd

import hetio.readwrite
import hetio.permute


# In[2]:


get_ipython().run_cell_magic('time', '', "msigdb_hetnet_path = os.path.join('hetnets', 'msigdb_hetnet.json.bz2')\ngraph = hetio.readwrite.read_graph(msigdb_hetnet_path)")


# In[3]:


num_permuted_hetnets = 5
num_swaps = 5


# In[4]:


get_ipython().run_cell_magic('time', '', "stat_dfs = list()\npermuted_graph = graph\n \nfor i in range(num_permuted_hetnets):\n    i += 1\n    print('Starting permutation', i)\n    permuted_graph, stats = hetio.permute.permute_graph(permuted_graph,\n                                                        multiplier=num_swaps,\n                                                        seed=i)\n    stat_df = pd.DataFrame(stats)\n    stat_df['permutation'] = i\n    stat_dfs.append(stat_df)\n    perm_path = os.path.join('hetnets', 'permuted',\n                             'msigdb_hetnet_perm-{}.json.bz2'.format(i))\n    hetio.readwrite.write_graph(permuted_graph, perm_path)")


# In[5]:


# Save stats
stat_df = pd.concat(stat_dfs)
stat_path = os.path.join('hetnets', 'permuted', 'stats.tsv')
stat_df.to_csv(stat_path, sep='\t', index=False, float_format='%.5g')

