"""
Gregory Way 2018
Interpret Compression
7.tcga-classify/scripts/tcga_util.py

Usage: For import only
"""


def build_feature_dictionary(dataset="TCGA"):
    """
    Generate a nested dictionary of the directory structure pointing to compressed
    feature matrices for training and testing sets

    Arguments:
    dataset - which dataset to load the z matrices from (default - "TCGA")

    Output: a nested dictionary storing the feature matrices from training and testing
    sets
    """
    import os
    import glob

    z_matrix_dict = {}
    for signal in ["signal", "shuffled"]:
        z_matrix_dict[signal] = {}

        if signal == "signal":
            results_dir = "{}_results".format(dataset)
        else:
            results_dir = "{}_shuffled_results".format(dataset)

        matrix_dir = os.path.join(
            "..", "2.ensemble-z-analysis", "results", results_dir, "ensemble_z_matrices"
        )

        for comp_dir in os.listdir(matrix_dir):
            matrix_comp_dir = os.path.join(matrix_dir, comp_dir)
            z_dim = comp_dir.split("_")[2]
            z_matrix_dict[signal][z_dim] = {}

            for z_file in glob.glob("{}/*_z_*".format(matrix_comp_dir)):
                seed = os.path.basename(z_file).split("_")[1]

                if seed not in z_matrix_dict[signal][z_dim].keys():
                    z_matrix_dict[signal][z_dim][seed] = {}

                if "_test_" in z_file:
                    z_matrix_dict[signal][z_dim][seed]["test"] = z_file
                else:
                    z_matrix_dict[signal][z_dim][seed]["train"] = z_file

    return z_matrix_dict


def build_top_feature_dictionary(algorithms, genes, num_features, load_random=True):
    """
    Generate a nested dictionary of the specific x matrices to train using

    Arguments:
    algorithms - list of algorithms to match ['pca', 'ica', 'nmf', 'dae', 'vae', 'all']
    genes - list of genes to match ['TP53', 'PTEN', 'KRAS', 'PIK3CA', 'TTN']
    num_features - list of the number of features used in prediction
    load_random - boolean if the random features should be loaded as well

    Output: a nested dictionary storing the x matrices for training and testing
    """

    import os
    import pandas as pd

    base_path = os.path.join("results", "top_feature_matrices")
    x_matrix_dict = {}
    for algorithm in algorithms:
        x_matrix_dict[algorithm] = {}

        for gene in genes:
            x_matrix_dict[algorithm][gene] = {}

            for n in num_features:
                x_matrix_dict[algorithm][gene][n] = {}

                base_file = os.path.join(
                    base_path,
                    "top_model_algorithm_{}_gene_{}_numtopfeatures_{}".format(
                        algorithm, gene, n
                    ),
                )

                train_file = "{}_train.tsv.gz".format(base_file)
                test_file = "{}_test.tsv.gz".format(base_file)

                train_df = pd.read_table(train_file, sep="\t", index_col=0)
                test_df = pd.read_table(test_file, sep="\t", index_col=0)

                x_matrix_dict[algorithm][gene][n]["train"] = train_df
                x_matrix_dict[algorithm][gene][n]["test"] = test_df

                if load_random and n == 200:
                    x_matrix_dict[algorithm][gene][n]["randomized"] = {}
                    base_file = "{}_randomized".format(base_file)

                    train_file = "{}_train.tsv.gz".format(base_file)
                    test_file = "{}_test.tsv.gz".format(base_file)

                    train_df = pd.read_table(train_file, sep="\t", index_col=0)
                    test_df = pd.read_table(test_file, sep="\t", index_col=0)

                    x_matrix_dict[algorithm][gene][n]["randomized"]["train"] = train_df
                    x_matrix_dict[algorithm][gene][n]["randomized"]["test"] = test_df

    return x_matrix_dict


def get_threshold_metrics(y_true, y_pred, drop=False):
    """
    Retrieve true/false positive rates and auroc/aupr for class predictions

    Arguments:
    y_true - an array of gold standard mutation status
    y_pred - an array of predicted mutation status
    drop - boolean if intermediate thresholds are dropped

    Output:
    dict of AUROC, AUPR, pandas dataframes of ROC and PR data, and cancer-type
    """
    import pandas as pd
    from sklearn.metrics import roc_auc_score, roc_curve
    from sklearn.metrics import precision_recall_curve, average_precision_score

    roc_columns = ["fpr", "tpr", "threshold"]
    pr_columns = ["precision", "recall", "threshold"]

    roc_results = roc_curve(y_true, y_pred, drop_intermediate=drop)
    roc_items = zip(roc_columns, roc_results)
    roc_df = pd.DataFrame.from_dict(dict(roc_items))

    prec, rec, thresh = precision_recall_curve(y_true, y_pred)
    pr_df = pd.DataFrame.from_records([prec, rec]).T
    pr_df = pd.concat([pr_df, pd.Series(thresh)], ignore_index=True, axis=1)
    pr_df.columns = pr_columns

    auroc = roc_auc_score(y_true, y_pred, average="weighted")
    aupr = average_precision_score(y_true, y_pred, average="weighted")

    return {"auroc": auroc, "aupr": aupr, "roc_df": roc_df, "pr_df": pr_df}


def summarize_results(
    results, gene_or_cancertype, signal, z_dim, seed, algorithm, data_type
):
    """
    Given an input results file, summarize and output all pertinent files

    Arguments:
    results - a results object output from `get_threshold_metrics`
    gene_or_cancertype - the gene or cancertype of interest
    signal - the signal of interest
    z_dim - the internal bottleneck dimension of the compression model
    seed - the seed used to compress the data
    algorithm - the algorithm used to compress the data
    data_type - the type of data (either training, testing, or cv)
    """

    results_append_list = [
        gene_or_cancertype,
        signal,
        z_dim,
        seed,
        algorithm,
        data_type,
    ]

    metrics_out_ = [results["auroc"], results["aupr"]] + results_append_list

    roc_df_ = results["roc_df"]
    pr_df_ = results["pr_df"]

    roc_df_ = roc_df_.assign(
        predictor=gene_or_cancertype,
        signal=signal,
        z_dim=z_dim,
        seed=seed,
        algorithm=algorithm,
        data_type=data_type,
    )

    pr_df_ = pr_df_.assign(
        predictor=gene_or_cancertype,
        signal=signal,
        z_dim=z_dim,
        seed=seed,
        algorithm=algorithm,
        data_type=data_type,
    )

    return metrics_out_, roc_df_, pr_df_


def extract_coefficients(cv_pipeline, feature_names, signal, z_dim, seed, algorithm):
    """
    Pull out the coefficients from the trained classifiers

    Arguments:
    cv_pipeline - the trained sklearn cross validation pipeline
    feature_names - the column names of the x matrix used to train model (features)
    results - a results object output from `get_threshold_metrics`
    gene - the gene of interest
    signal - the signal of interest
    z_dim - the internal bottleneck dimension of the compression model
    seed - the seed used to compress the data
    algorithm - the algorithm used to compress the data
    """
    import pandas as pd

    final_pipeline = cv_pipeline.best_estimator_
    final_classifier = final_pipeline.named_steps["classify"]

    coef_df = pd.DataFrame.from_dict(
        {"feature": feature_names, "weight": final_classifier.coef_[0]}
    )

    coef_df = (
        coef_df.assign(abs=coef_df["weight"].abs())
        .sort_values("abs", ascending=False)
        .reset_index(drop=True)
        .assign(signal=signal, z_dim=z_dim, seed=seed, algorithm=algorithm)
    )

    return coef_df


def process_y_matrix(
    y_mutation,
    y_copy,
    include_copy,
    gene,
    sample_freeze,
    mutation_burden,
    filter_count,
    filter_prop,
    output_directory,
    hyper_filter=5,
):
    """
    Combine copy number and mutation data and filter cancer-types to build y matrix

    Arguments:
    y_mutation - Pandas DataFrame of mutation status
    y_copy - Pandas DataFrame of copy number status
    include_copy - boolean if the copy number data should be included in status calc
    gene - string indicating gene of interest (used for writing proportion file)
    sample_feeze - pandas dataframe storing which samples to use
    mutation_burden - pandas dataframe storing log10 mutation counts
    filter_count - the number of positives or negatives required per cancer-type
    filter_prop - the proportion of positives or negatives required per cancer-type
    output_directory - the name of the directory to store the gene summary
    hyper_filter - the number of std dev above log10 mutation burden to filter

    Output:
    Write file of cancer-type filtering to disk and output a processed y vector
    """
    import os
    import pandas as pd

    if include_copy:
        y_df = y_copy + y_mutation
    else:
        y_df = y_mutation

    y_df.loc[y_df > 1] = 1
    y_df = pd.DataFrame(y_df)
    y_df.columns = ["status"]

    y_df = (
        y_df.merge(
            sample_freeze, how="left", left_index=True, right_on="SAMPLE_BARCODE"
        )
        .set_index("SAMPLE_BARCODE")
        .merge(mutation_burden, left_index=True, right_index=True)
    )

    # Get statistics per gene and disease
    disease_counts_df = pd.DataFrame(y_df.groupby("DISEASE").sum()["status"])

    disease_proportion_df = disease_counts_df.divide(
        y_df["DISEASE"].value_counts(sort=False).sort_index(), axis=0
    )

    # Filter diseases with low counts or proportions for classification balance
    filter_disease_df = (disease_counts_df > filter_count) & (
        disease_proportion_df > filter_prop
    )
    filter_disease_df.columns = ["disease_included"]

    disease_stats_df = disease_counts_df.merge(
        disease_proportion_df,
        left_index=True,
        right_index=True,
        suffixes=("_count", "_proportion"),
    ).merge(filter_disease_df, left_index=True, right_index=True)

    filter_file = "{}_filtered_cancertypes.tsv".format(gene)
    filter_file = os.path.join(output_directory, filter_file)
    disease_stats_df.to_csv(filter_file, sep="\t")

    # Filter
    use_diseases = disease_stats_df.query("disease_included").index.tolist()
    burden_filter = y_df["log10_mut"] < hyper_filter * y_df["log10_mut"].std()
    y_df = y_df.loc[burden_filter, :].query("DISEASE in @use_diseases")

    return y_df


def process_y_matrix_cancertype(
    acronym, sample_freeze, mutation_burden, hyper_filter=5
):
    """
    Build a y vector based on cancer-type membership

    Arguments:
    acronym - the TCGA cancer-type barcode
    sample_freeze - a dataframe storing TCGA barcodes and cancer-types
    mutation_burden - a log10 mutation count per sample (added as covariate)

    Output:
    A y status DataFrame and a status count dataframe
    """

    import pandas as pd

    y_df = sample_freeze.assign(status=0)
    y_df.loc[y_df.DISEASE == acronym, "status"] = 1

    y_df = y_df.set_index("SAMPLE_BARCODE").merge(
        mutation_burden, left_index=True, right_index=True
    )

    burden_filter = y_df["log10_mut"] < hyper_filter * y_df["log10_mut"].std()
    y_df = y_df.loc[burden_filter, :]

    count_df = pd.DataFrame(y_df.status.value_counts()).reset_index()
    count_df.columns = ["status", acronym]

    return y_df, count_df


def align_matrices(x_file_or_df, y, add_cancertype_covariate=True, algorithm=None):
    """
    Process the x matrix for the given input file and align x and y together

    Arguments:
    x_file_or_df - string location of the x matrix or matrix df itself
    y - pandas DataFrame storing status of corresponding samples
    algorithm - a string indicating which algorithm to subset the z matrices

    Output:
    The samples used to subset and the processed X and y matrices
    """
    import pandas as pd
    from sklearn.preprocessing import StandardScaler

    # Load Data
    try:
        x_df = pd.read_table(x_file_or_df, index_col=0)
        if algorithm:
            x_df = x_df.loc[:, x_df.columns.str.contains(algorithm)]
    except:
        x_df = x_file_or_df

    # Subset samples
    use_samples = set(y.index).intersection(set(x_df.index))

    x_df = x_df.reindex(use_samples)
    y = y.reindex(use_samples)

    # Transform features to between zero and one
    x_scaled = StandardScaler().fit_transform(x_df)
    x_df = pd.DataFrame(x_scaled, columns=x_df.columns, index=x_df.index)

    # create covariate info
    mutation_covariate_df = pd.DataFrame(y.loc[:, "log10_mut"], index=y.index)

    # Merge log10 mutation burden covariate
    x_df = x_df.merge(mutation_covariate_df, left_index=True, right_index=True)

    if add_cancertype_covariate:
        # Merge features with covariate data
        covariate_df = pd.get_dummies(y.DISEASE)
        x_df = x_df.merge(covariate_df, left_index=True, right_index=True)

    return use_samples, x_df, y


def train_model(x_train, x_test, y_train, alphas, l1_ratios, n_folds=5, max_iter=1000):
    """
    Build the logic and sklearn pipelines to train x matrix based on input y

    Arguments:
    x_train - pandas DataFrame of feature matrix for training data
    x_test - pandas DataFrame of feature matrix for testing data
    y_train - pandas DataFrame of processed y matrix (output from align_matrices())
    alphas - list of alphas to perform cross validation over
    l1_ratios - list of l1 mixing parameters to perform cross validation over
    n_folds - int of how many folds of cross validation to perform
    max_iter - the maximum number of iterations to test until convergence

    Output:
    The full pipeline sklearn object and y matrix predictions for training, testing,
    and cross validation
    """
    from sklearn.model_selection import cross_val_predict
    from sklearn.pipeline import Pipeline
    from sklearn.linear_model import SGDClassifier
    from dask_searchcv import GridSearchCV

    # Setup the classifier parameters
    clf_parameters = {
        "classify__loss": ["log"],
        "classify__penalty": ["elasticnet"],
        "classify__alpha": alphas,
        "classify__l1_ratio": l1_ratios,
    }

    estimator = Pipeline(
        steps=[
            (
                "classify",
                SGDClassifier(
                    random_state=0,
                    class_weight="balanced",
                    loss="log",
                    max_iter=max_iter,
                    tol=1e-3,
                ),
            )
        ]
    )

    cv_pipeline = GridSearchCV(
        estimator=estimator,
        param_grid=clf_parameters,
        n_jobs=-1,
        cv=n_folds,
        scoring="roc_auc",
        return_train_score=True,
    )

    # Fit the model
    cv_pipeline.fit(X=x_train, y=y_train.status)

    # Obtain cross validation results
    y_cv = cross_val_predict(
        cv_pipeline.best_estimator_,
        X=x_train,
        y=y_train.status,
        cv=n_folds,
        method="decision_function",
    )

    # Get all performance results
    y_predict_train = cv_pipeline.decision_function(x_train)
    y_predict_test = cv_pipeline.decision_function(x_test)

    return cv_pipeline, y_predict_train, y_predict_test, y_cv


def check_status(file):
    """
    Check the status of a gene or cancer-type application

    Arguments:
    file - the file to check if it exists. If exists, then there is no need to rerun

    Output:
    boolean if the file exists or not
    """

    import os

    return os.path.isfile(file)


def get_feature(z_dim, seed, feature):
    """
    Given arguments, find the specific compressed feature for test and train data

    Arguments:
    z_dim - string indicating the bottleneck dimension
    seed - string indicating which seed to search for model
    feature - string indicating algorithm and number of specific feature

    Output:
    Reads compressed z matrices for training and testing and outputs them as a tuple
    of two pandas series.
    """
    import os
    import pandas as pd

    base_file = os.path.join(
        "..",
        "2.ensemble-z-analysis",
        "results",
        "TCGA_results",
        "ensemble_z_matrices",
        "tcga_components_{}".format(z_dim),
    )
    train_file = os.path.join(base_file, "model_{}_z_matrix.tsv.gz".format(seed))
    test_file = os.path.join(base_file, "model_{}_z_test_matrix.tsv.gz".format(seed))

    test_df = pd.read_table(test_file, index_col=0).loc[:, feature].sort_index()
    train_df = pd.read_table(train_file, index_col=0).loc[:, feature].sort_index()

    return test_df, train_df


def build_good_feature_matrix(coef_df):
    """
    Compile training and testing feature matrices (X) given a coefficient DataFrame
    output from the various mutation classification analyses through repeated calls to
    `get_feature()`

    Arguments:
    coef_df - an input dataframe with the following information:
        z_dim, seed, and feature columns

    Output:
    Reads a series of compressed z matrices for training and testing and combines them
    into two dataframes. The dataframes represent the training and testing matrices
    used in downstream analyses.
    """
    import pandas as pd

    all_test_features = []
    all_train_features = []
    for feature_idx, feature_row in coef_df.iterrows():
        z_dim = feature_row.z_dim
        seed = feature_row.seed
        feature = feature_row.feature

        test_feature_df, train_feature_df = get_feature(z_dim, seed, feature)
        all_test_features.append(test_feature_df)
        all_train_features.append(train_feature_df)

    all_test_features_df = pd.concat(all_test_features, axis="columns")
    all_train_features_df = pd.concat(all_train_features, axis="columns")

    return all_test_features_df, all_train_features_df
