"""
2018 Gregory Way
5.analyze-components/scripts/latent.py

Functions to perform gene set analyses on compressed gene expression weights

Usage:

    import only

    from scripts.latent import latentModel

Output:

A python class with several attributes and methods to assign biology to
compressed gene expression features.
"""

import os
import glob
import gseapy as gp
import pandas as pd


class latentModel():
    """
    Methods for processing a weight matrix from a compression model

    Usage:
    from scripts.latent import latentModel
    w = latentModel(filename)
    """
    def __init__(self, filename, z_dim=None, dataset_name=None,
                 algorithm_name=None, weight_seed=None, shuffled_true=False):
        """
        The latent model is initialized with a filename that will load a weight
        matrix and set appropriate attributes.

        Arguments:
        filename - will store the filename for the selected weight matrix
        z_dim - the dimensionality of the bottleneck layer
        algorithm_name - the name of the compression algorithm
        weight_seed - the seed used to compile the specific weight matrix
        shuffled_true - boolean if the data were shuffled prior to compression
        """
        # Load gene expression data
        self.filename = filename
        self.w_df = pd.read_table(self.filename, index_col=0)

        # Set attributes
        self.z_dim = z_dim
        self.dataset_name = dataset_name.upper()
        self.weight_seed = int(weight_seed)
        self.shuffled_true = shuffled_true
        self.is_gene_dictionary_loaded = False

        if algorithm_name:
            self.algorithm_name = algorithm_name.lower()
            use_cols = self.w_df.columns.str.contains(self.algorithm_name)
            self.w_df = self.w_df.loc[:, use_cols]

    def _split_genes(self, pos_w_df, neg_w_df):
        """
        Internal method to process positive and negative tailed genes by some
        prior determination

        Arguments:
        pos_w_df - a boolean matrix of genes in the positive tail
        neg_w_df - a boolean matrix of genes in the negative tail

        Output:
        A long dataframe of positive and negative genes
        """

        rename_cols = ['entrez_gene', 'feature', 'feature_weight']

        pos_w_df = pos_w_df.stack().reset_index()
        pos_w_df.columns = rename_cols
        pos_w_df = pos_w_df.dropna()

        neg_w_df = neg_w_df.stack().reset_index()
        neg_w_df.columns = rename_cols
        neg_w_df = neg_w_df.dropna()

        return pos_w_df, neg_w_df

    def _load_gene_dict(self, hash='ad9631bb4e77e2cdc5413b0d77cb8f7e93fc5bee'):
        """
        Internal method to load gene dictionary from
        https://github.com/cognoma/genes

        Arguments:
        hash - a string indicating which commit hash version to load data
        """

        url = 'https://raw.githubusercontent.com/cognoma/genes/' + \
            '{}/data/genes.tsv'.format(hash)
        self.gene_dictionary_df = (
            pd.read_table(url)
            .query("gene_type == 'protein-coding'")
        )

        self.is_gene_dictionary_loaded = True

    def get_high_weight_genes(self, std_dev=2.5, long_format=False):
        """
        Method to obtain high weight genes for all latent space features in a
        given weight matrix

        Arguments:
        std_dev - how many standard deviations above the mean to consider
        long_format - boolean to output data in wide or long format
        """

        feature_mean = self.w_df.mean()
        std_dev_cutoff = std_dev * self.w_df.std()

        self.pos_high_w_df = feature_mean + std_dev_cutoff
        self.neg_high_w_df = feature_mean - std_dev_cutoff

        self.pos_high_w_df = self.w_df.ge(self.pos_high_w_df)
        self.neg_high_w_df = self.w_df.le(self.neg_high_w_df)

        if long_format:
            self.pos_high_w_df, self.neg_high_w_df = (
                self._split_genes(
                    pos_w_df=self.pos_high_w_df,
                    neg_w_df=self.neg_high_w_df
                )
            )

        self.high_weight_std_cutoff = std_dev

    def split_pos_neg_genes(self, long_format=False):
        """
        Method to split positive and negative weight genes into two features

        Arguments:
        long_format - boolean to output data in wide or long format
        """

        self.pos_w_df = self.w_df[self.w_df > 0]
        self.neg_w_df = self.w_df[self.w_df < 0]

        if long_format:
            self.pos_w_df, self.neg_w_df = (
                self._split_genes(
                    pos_w_df=self.pos_w_df,
                    neg_w_df=self.neg_w_df
                )
            )

    def get_squared_distribution(self):
        """
        Method to process distribution of feature scores by squaring
        """

        self.sqr_w_df = self.w_df ** 2

    def translate_gene_ids(self, gene_list, translate_from, translate_to,
                           **kwargs):
        """
        Method to translate from entrez gene ids to other symbols of choice

        Arguments:
        gene_list - a list or pandas series of gene identifiers to translate
        translate_from - string indicating the symbol type to translate from
        translate_to - string indicating the symbol type to translate to

        Return:
        Outputs a translated list of gene IDs
        """

        # unpack kwargs
        hash = kwargs.pop('hash', 'ad9631bb4e77e2cdc5413b0d77cb8f7e93fc5bee')

        # The translation must be limited
        choices = ['entrez_gene_id', 'symbol']
        if translate_from not in choices or translate_to not in choices:
            err_txt = "Translation for 'entrez_gene_id' and 'symbol' only"
            raise ValueError(err_txt)

        if not self.is_gene_dictionary_loaded:
            self._load_gene_dict(hash=hash)

        gene_mapper = dict(zip(self.gene_dictionary_df.loc[:, translate_from],
                               self.gene_dictionary_df.loc[:, translate_to]))

        gene_list = pd.Series(gene_list)
        return gene_list.map(gene_mapper)

    def _get_compressed_feature(self, weight_df, subset_algs, current_z):
        """
        Helper function to isolate the feature of interest

        Arguments:
        weight_df - a gene by compressed feature pandas DataFrame
        subset_algs - a boolean vector indicating which columns correpond to
                      the algorithm of focus
        current_z - int describing which z vector to subset

        Output:
        A pandas Series of the compressed feature with a specific distribution
        """

        compressed_feature_df = (
            weight_df
            .loc[:, subset_algs]
            .iloc[:, current_z]
        )

        return compressed_feature_df

    def _get_gsea_results(self, use_df, algorithm, current_z, distribution,
                          num_perm, shuffled, direction, seed):
        """
        Helper function to perform and process the GSEA results

        Arguments:
        use_df - a compressed feature
        algorithm - a string indicating which algorithm is being analyzed
        current_z - int describing which z vector to subset
        distribution - a string indicating which distribution is being tested
        num_perm - int specifying the number of permutations to use
        shuffled - a boolean indicating if data is shuffled before training
        direction - string indicating the distribution tail being considered
        seed - the random seed used to train the model

        Output:
        A pandas DataFrame of processed results and attributes
        """
        results_df = (
            self.run_gsea_prerank(
                gene_score_df=use_df,
                permutation_num=num_perm)
            )

        # Store various attributes for downstream analysis
        results_df = (
            results_df.assign(
                distrib=distribution,
                algorithm=algorithm,
                current_z=current_z,
                shuffled=shuffled,
                direction=direction,
                seed=seed
                )
            )

        return results_df

    def get_gsea_compressed_matrix(self, algorithms, distrib_methods,
                                   shuffled, seed, num_perm=15):
        """
        Perform GSEA on all compressed features in a weight matrix across all
        algorithms and distribution methods

        Arguments:
        lm - a latent model storing weight matrices and several attributes
        algorithms - a list of algorithms to subset and analyze
        distrib_methods - a list of methods to obtain compressed feature genes
        shuffled - a boolean indicating if data is shuffled before training
        seed - the random seed that the model is trained using
        num_perm - int specifying the number of permutations to use

        Output:
        a gseapy res2d object with significant pathways (fdr) assigned to
        compressed features
        """

        current_z_list = []
        # Loop through all of the algorithms
        for alg in algorithms:

            # Determine the first subsetting logic
            subset_algs = self.w_df.columns.str.contains(alg)

            # Loop through the current z assignment
            # take each column individually
            for current_z in range(0, int(self.z_dim)):

                # Now, loop through each of the distributions
                for distrib in distrib_methods:
                    print(alg, current_z, distrib, seed)

                    # Full distributions do not require additional processing
                    num_distrib = 1
                    if distrib == 'full':
                        use_df = self._get_compressed_feature(
                            weight_df=self.w_df,
                            subset_algs=subset_algs,
                            current_z=current_z
                        )
                    elif distrib == 'full_squared':
                        use_df = self._get_compressed_feature(
                            weight_df=self.sqr_w_df,
                            subset_algs=subset_algs,
                            current_z=current_z
                        )
                    elif distrib == 'pos_neg':
                        num_distrib = 2
                        use_pos_df = self._get_compressed_feature(
                            weight_df=self.pos_w_df,
                            subset_algs=subset_algs,
                            current_z=current_z
                        )
                        use_neg_df = self._get_compressed_feature(
                            weight_df=self.neg_w_df,
                            subset_algs=subset_algs,
                            current_z=current_z
                        )

                    elif distrib == 'pos_neg_high_weight':
                        num_distrib = 2
                        use_pos_df = self._get_compressed_feature(
                            weight_df=self.pos_high_w_df,
                            subset_algs=subset_algs,
                            current_z=current_z
                        )
                        use_neg_df = self._get_compressed_feature(
                            weight_df=self.neg_high_w_df,
                            subset_algs=subset_algs,
                            current_z=current_z
                        )

                    # Perform GSEA on the single distribution
                    if num_distrib == 1:
                        results_df = (
                            self._get_gsea_results(
                                use_df=use_df,
                                algorithm=alg,
                                current_z=current_z,
                                distribution=distrib,
                                direction='both',
                                num_perm=num_perm,
                                shuffled=shuffled,
                                seed=seed)
                        )

                        current_z_list += [results_df]
                    # Or, perform GSEA on the double distributions, if exist
                    elif num_distrib == 2:

                        if len(use_pos_df.dropna()) != 0:
                            results_df = (
                                self._get_gsea_results(
                                    use_df=use_df,
                                    algorithm=alg,
                                    current_z=current_z,
                                    distribution=distrib,
                                    direction='positive',
                                    num_perm=num_perm,
                                    shuffled=shuffled,
                                    seed=seed)
                            )

                            current_z_list += [results_df]

                        if len(use_neg_df.dropna()) != 0:
                            results_df = (
                                self._get_gsea_results(
                                    use_df=use_df,
                                    algorithm=alg,
                                    current_z=current_z,
                                    distribution=distrib,
                                    direction='negative',
                                    num_perm=num_perm,
                                    shuffled=shuffled,
                                    seed=seed)
                            )

                            current_z_list += [results_df]

    def run_gsea_prerank(self, gene_score_df, **kwargs):
        """
        Method to perform ranked GSEA on an input gene list and score

        Arguments:
        gene_score_df - A two column pandas DataFrame (no column names) with
                        the first column storing Hugo gene symbols and the
                        second column storing the scores

        Output:
        A gseapy res2d dataframe
        https://gseapy.readthedocs.io/en/master/run.html#gseapy.prerank

        NOTE: To get list of all available genesets run:

            import gseapy
            names = gseapy.get_library_name()
            print(names)
        """

        gene_score_df.columns = [0, 1]

        # unpack kwargs
        gene_sets = kwargs.pop('gene_sets', 'KEGG_2016')
        processes = kwargs.pop('processes', 1)
        permutation_num = kwargs.pop('permutation_num', 1000)
        verbose = kwargs.pop('verbose', False)
        outdir = kwargs.pop('outdir', None)
        format = kwargs.pop('format', 'png')
        no_plot = kwargs.pop('no_plot', True)
        fdr_cutoff = kwargs.pop('fdr_cutoff', 0.05)

        # Perform the GSEA prerank
        result = gp.prerank(rnk=gene_score_df,
                            gene_sets=gene_sets,
                            processes=processes,
                            permutation_num=permutation_num,
                            outdir=outdir,
                            format=format,
                            no_plot=no_plot,
                            verbose=verbose)

        result = result.res2d
        result = (
            result.loc[:, ['fdr', 'geneset_size']]
            .query("fdr < @fdr_cutoff")
        )

        return result


def run_gsea_pipeline_command(input_weight_dir, z_dim, dataset_name, num_perm,
                              shuffled_true, algorithms, distrib_methods):
    """
    Perform the entire pipeline to obtain all GSEA results built into a single
    function for multiprocessing to act upon

    Arguments:
    input_z_dir - a folder storing how many dimensions are being compressed
    z_dim - the dimensionality of the bottleneck layer
    dataset_name - the name of the dataset
    num_perm - int specifying the number of permutations to use
    shuffled_true - boolean indicating if the data were shuffled
    algorithms - a list of algorithms to analyze
    distrib_methods - a list of methods to obtain compressed feature gene lists

    Output:
    for each input z directory, a results file will be written to disk
    """

    # Step 0 - Compile all of the weight matrices in a dictionary
    weight_matrices = {}
    for file_name in glob.glob('{}/*_weight_matrix*'.format(input_weight_dir)):
        seed = file_name.split('/')[-1].split('_')[1]
        weight_matrices[seed] = file_name

    # Loop over all seeds in the directory
    gsea_results_list = []
    for seed in weight_matrices.keys():

        # Step 1 - instatiate the latent model
        file = weight_matrices[seed]
        lm = latentModel(filename=file,
                         z_dim=z_dim,
                         dataset_name=dataset_name,
                         algorithm_name=None,
                         weight_seed=seed,
                         shuffled_true=shuffled_true)

        # Step 2 - load the gene dictionary that will help with translation
        lm._load_gene_dict()

        # Step 3 - Translate the indeces of the weight matrix
        lm.w_df.index = (
            lm.translate_gene_ids(
                gene_list=lm.w_df.index,
                translate_from='entrez_gene_id',
                translate_to='symbol')
        )

        # Step 4 - Process the distriubtions that will be used for comparisons
        lm.get_high_weight_genes()
        lm.split_pos_neg_genes()
        lm.get_squared_distribution()

        gsea_results_df = (
            lm.get_gsea_compressed_matrix(
                algorithms=algorithms,
                distrib_methods=distrib_methods,
                num_perm=num_perm,
                shuffled=shuffled_true,
                seed=seed
                )
            )

        gsea_results_list += gsea_results_df

    out_file = 'gsea_results_{}_zdim_{}.tsv'.format(dataset_name, z_dim)
    out_file = os.path.join('results', 'gsea', out_file)
    pd.concat(gsea_results_list).to_csv(out_file, sep='\t')
