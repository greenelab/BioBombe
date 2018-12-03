"""
Gregory Way 2018
Interpret Compression
7.tcga-classify/scripts/tcga_util.py

Usage: For import only
"""


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


def summarize_results(results, gene, signal, z_dim, seed, algorithm, data_type):
    """
    Given an input results file, summarize and output all pertinent files

    Arguments:
    results - a results object output from `get_threshold_metrics`
    gene - the gene of interest
    signal - the signal of interest
    z_dim - the internal bottleneck dimension of the compression model
    seed - the seed used to compress the data
    algorithm - the algorithm used to compress the data
    data_type - the type of data (either training, testing, or cv)
    """

    results_append_list = [gene, signal, z_dim, seed, algorithm, data_type]

    metrics_out_ = [results["auroc"], results["aupr"]] + results_append_list

    roc_df_ = results["roc_df"]
    pr_df_ = results["pr_df"]

    roc_df_ = roc_df_.assign(
        gene=gene,
        signal=signal,
        z_dim=z_dim,
        seed=seed,
        algorithm=algorithm,
        data_type=data_type,
    )

    pr_df_ = pr_df_.assign(
        gene=gene,
        signal=signal,
        z_dim=z_dim,
        seed=seed,
        algorithm=algorithm,
        data_type=data_type,
    )

    return metrics_out_, roc_df_, pr_df_


def extract_coefficients(
    cv_pipeline, feature_names, gene, signal, z_dim, seed, algorithm
):
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
        .assign(gene=gene, signal=signal, z_dim=z_dim, seed=seed, algorithm=algorithm)
    )


def process_y_matrix(
    y_copy,
    y_mutation,
    gene,
    sample_freeze,
    mutation_burden,
    filter_count,
    filter_prop,
    hyper_filter=5,
):
    """
    Combine copy number and mutation data and filter cancer-types to build y matrix

    Arguments:
    y_copy - Pandas DataFrame of copy number status
    y_mutation - Pandas DataFrame of mutation status
    gene - string indicating gene of interest (used for writing proportion file)
    sample_feeze - pandas dataframe storing which samples to use
    mutation_burden - pandas dataframe storing log10 mutation counts
    filter_count - the number of positives or negatives required per cancer-type
    filter_prop - the proportion of positives or negatives required per cancer-type
    hyper_filter - the number of std dev above log10 mutation burden to filter

    Output:
    Write file of cancer-type filtering to disk and output a processed y vector
    """
    import os
    import pandas as pd

    y_df = y_copy + y_mutation
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
    filter_file = os.path.join("results", filter_file)
    disease_stats_df.to_csv(filter_file, sep="\t")

    # Filter
    use_diseases = disease_stats_df.query("disease_included").index.tolist()
    burden_filter = y_df["log10_mut"] < hyper_filter * y_df["log10_mut"].std()
    y_df = y_df.loc[burden_filter, :].query("DISEASE in @use_diseases")

    return y_df


def align_matrices(x_file, y, algorithm=None):
    """
    Process the x matrix for the given input file and align x and y together

    Arguments:
    x_file - string location of the x matrix
    y - pandas DataFrame storing status of corresponding samples
    algorithm - a string indicating which algorithm to subset the z matrices

    Output:
    The samples used to subset and the processed X and y matrices
    """
    import pandas as pd
    from sklearn.preprocessing import StandardScaler

    # Load Data
    x_df = pd.read_table(x_file, index_col=0)
    if algorithm:
        x_df = x_df.loc[:, x_df.columns.str.contains(algorithm)]

    # Subset samples
    use_samples = set(y.index).intersection(set(x_df.index))

    x_df = x_df.reindex(use_samples)
    y = y.reindex(use_samples)

    # Transform features to between zero and one
    x_scaled = StandardScaler().fit_transform(x_df)
    x_df = pd.DataFrame(x_scaled, columns=x_df.columns, index=x_df.index)

    # create covariate info
    covariate_df = pd.get_dummies(y.DISEASE)

    mutation_covariate_df = pd.DataFrame(y.loc[:, "log10_mut"], index=y.index)

    # Merge features with covariate data
    x_df = x_df.merge(covariate_df, left_index=True, right_index=True).merge(
        mutation_covariate_df, left_index=True, right_index=True
    )

    return use_samples, x_df, y


def train_model(x_train, x_test, y_train, alphas, l1_ratios, n_folds=5):
    """
    Build the logic and sklearn pipelines to train x matrix based on input y

    Arguments:
    x_train - pandas DataFrame of feature matrix for training data
    x_test - pandas DataFrame of feature matrix for testing data
    y_train - pandas DataFrame of processed y matrix (output from align_matrices())
    alphas - list of alphas to perform cross validation over
    l1_ratios - list of l1 mixing parameters to perform cross validation over
    n_folds - int of how many folds of cross validation to perform

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
                    max_iter=100,
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
