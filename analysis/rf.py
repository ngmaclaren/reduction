# WARNING: the .fit() and permutation_importance() steps are taking some time. This may be slow.

# the RFR will want an X array (n, m) and a y array (n, )
# So, get a data.frame/csv from R. That will be the node features, X.
# The y will be the freq counts. These can come from R, too, and be in the same csv. Can simply select w/ pandas.
# Then, go through the example pipeline. May need to adjust hyperparameters, but shouldn't be too difficult.

import time
import joblib
N_CORES = joblib.cpu_count(only_physical_cores=True)
import numpy as np
import pandas as pd
from sklearn.preprocessing import OrdinalEncoder # ordinal encoder fine for regression trees, apparently
from sklearn.model_selection import train_test_split # need to check performance somehow?
# from sklearn.pipeline import make_pipeline # not using, not sure why
from sklearn.ensemble import RandomForestRegressor # the classifier itself
from sklearn.inspection import permutation_importance # skl's preferred variable importance measure

### This is a list of RFR arguments I might want to adjust
# n_estimators : default=100, number of trees
# set criterion="poisson", default is "squared_error"
# max_depth : default=None, nodes are expanded until all leaves are pure or contain less than min_samples_split
# min_samples_split : default=2
# max_features : try "sqrt" or "log2", default=1.0 is probably too large. Might want 0.3, but more recent work apparently recommends max_features=1.0
# oob_scores : set True, probably
# n_jobs : set 2, probably; for parallelization
# max_samples : not sure if want to mess with this. Better to provide train and test data? KFold CV?
# It is recommended to adjust max_depth, min_samples_leaf, etc. to avoid overlarge memory usage

networks = [
    "dolphin", "proximity", "celegans", "euroroad", "email", "gkk", "ba", "hk", "lfr", "er",
    "drosophila", "reactome", "route_views", "spanish", "foldoc", "tree_of_life",  "word_assoc",
    "enron", "marker_cafe", "prosper"
]

results = []
for network in networks:
    # network = "prosper"
    csvfile = f"../data/nodefeatures-{network}.csv"

    df = pd.read_csv(csvfile)

    encoder = OrdinalEncoder().set_output(transform="pandas")
    df["dyn"] = encoder.fit_transform(df[["dynamics"]])

    xcol = ["k", "knn", "lcl", "cc", "bc", "kcore", "dyn"]
    ycol = "Freq"
    X = df[xcol]
    y = df[ycol]

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.33)
    model = RandomForestRegressor(criterion="poisson", n_jobs=N_CORES, max_features="log2", oob_score=True)

    t_rf_start = time.time()
    rf = model.fit(X_train, y_train)
    t_rf_stop = time.time()
    t_rf = t_rf_stop - t_rf_start

    t_imp_start = time.time()
    imp = permutation_importance(rf, X_test, y_test)
    t_imp_stop = time.time()
    t_imp = t_imp_stop - t_imp_start

    print(f"Network: {network}")
    print(f"Model wall time: {round(t_rf, 3)} seconds.")
    print(f"Importance computation wall time: {round(t_imp, 3)} seconds")
    # print(f"Coefficient of determination: {round(rf.score(X_test, y_test), 3)}") # default is "r2_score", which is 1 - (squared error)/(variance), where squared error is (predicted y - true y)^2. Docs report that r2_score can be negative but max is 1.0 (the case when predicted y == true y for all i.
    # print(f"Out-of-bag error: {round(rf.oob_score_, 3)}")
    # print("Importances:")
    # print(pd.Series(imp.importances_mean, index = xcol).round(3))

    # test = model.fit(X, y)
    # print(f"Out-of-bag error on the full data: {round(test.oob_score_, 3)}")

    # Ok, this basically works. I want to run this model on each network, saving:
    # - coefficient of determination
    # - importance
    imps = pd.Series(imp.importances_mean, index = xcol).to_dict()
    r2 = {"r2": rf.score(X_test, y_test)}
    info = {"N": int(df["N"][0]), "M": int(df["M"][0])}
    res = pd.DataFrame({network: r2 | imps | info}).transpose() # now it's a row
    results.append(res)


rdf = pd.concat(results)
rdf.to_csv("../data/rf-importances.csv")
