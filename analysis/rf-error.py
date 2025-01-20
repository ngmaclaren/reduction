import time
import joblib
N_CORES = joblib.cpu_count(only_physical_cores=True)
import numpy as np
import pandas as pd
from sklearn.preprocessing import OrdinalEncoder # ordinal encoder fine for regression trees, apparently
from sklearn.model_selection import train_test_split # need to check performance somehow?
from sklearn.ensemble import RandomForestRegressor # the classifier itself
from sklearn.inspection import permutation_importance # skl's preferred variable importance measure

separate_dynamics = True # False

def analyze(X, y, show_summary=False):
    """
    X and y are selected from df and are the predictor and outcome variables, respectively.
    Returns a data frame row with results (R^2, variable importances).
    """
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.1, random_state=123)
    model = RandomForestRegressor(criterion="squared_error",
                                  max_features="sqrt",
                                  oob_score=True,
                                  n_jobs=N_CORES, random_state=456)
    
    t_rf_start = time.time()
    rf = model.fit(X_train, y_train)
    t_rf_stop = time.time()
    t_rf = t_rf_stop - t_rf_start
    
    t_imp_start = time.time()
    ## For the below, note that the two vectors of importances can be somewhat different. I believe this is because the model itself fits poorly.
    # imp = permutation_importance(rf, X_test, y_test) # results in more negative importances
    imp = permutation_importance(rf, X_train, y_train) # results in more positive importances
    t_imp_stop = time.time()
    t_imp = t_imp_stop - t_imp_start

    imps = pd.Series(imp.importances_mean, index = xcol).to_dict()
    r2 = {"r2": rf.score(X_test, y_test)}
    res = pd.DataFrame({network: r2 | imps | info}).transpose() # now it's a row

    if show_summary:
        print(f"Network: {network}")
        print(f"Dynamics: {dyn}")
        print(f"Model wall time: {round(t_rf, 3)} seconds.")
        print(f"Importance computation wall time: {round(t_imp, 3)} seconds")
        print(f"Coefficient of determination: {round(rf.score(X_test, y_test), 3)}") # default is "r2_score", which is 1 - (squared error)/(variance), where squared error is (predicted y - true y)^2. Docs report that r2_score can be negative but max is 1.0 (the case when predicted y == true y for all i.
        print(f"Out-of-bag error: {round(rf.oob_score_, 3)}")
        print("Importances:")
        print(pd.Series(imp.importances_mean, index = xcol).round(3))
    else:
        print(f"{network} network, {dyn} dynamics, complete.")
    
    return res
#results.append(res)    

networks = [
    "dolphin", "proximity", "celegans", "euroroad", "email", "gkk", "ba", "hk", "lfr", "er",
    "drosophila", "reactome", "route_views", "spanish", "foldoc", "tree_of_life",  "word_assoc",
    "enron", "marker_cafe", "prosper"
]
dynamics = ["doublewell", "mutualistic", "SIS", "genereg"]

results = []
for network in networks:
    # network = "route_views"
    csvfile = f"../data/nodefeatures-{network}.csv" # will cause problems on Windows?

    df = pd.read_csv(csvfile)
    df["logerror"] = np.log(df["error"])
    ycol = "logerror"

    info = {"N": int(df["N"][0]), "M": int(df["M"][0])}

    encoder = OrdinalEncoder().set_output(transform="pandas")
    df["dyn"] = encoder.fit_transform(df[["dynamics"]])

    if separate_dynamics:
        for dyn in dynamics:
            # dyn = "genereg"
            sdf = df.loc[df["dynamics"] == dyn]
            info["dynamics"] = dyn
            xcol = ["k", "knn", "lcl", "cc", "bc", "kcore"]
            X = sdf[xcol]
            y = sdf[ycol]

            results.append(analyze(X, y))
    else:
        dyn = "all"
        xcol = ["k", "knn", "lcl", "cc", "bc", "kcore", "dyn"]
        X = df[xcol]
        y = df[ycol]

        results.append(analyze(X, y))

rdf = pd.concat(results)
col_order = ["N", "r2", "k", "cc", "bc", "knn", "lcl", "kcore", "dynamics"]
# rdf.to_csv("../data/rf-error-importances.csv")
rdf[col_order].drop("er").sort_values(["dynamics", "N"]).to_csv("../data/rf-error-importances.csv")

rdf[col_order].drop("er").round(3).sort_values(["dynamics", "N"])

#### Sanity check:
## Use V and w instead of X and y
# Vcol = ["k", "knn", "lcl", "cc", "bc", "Freq"]
# wcol = "kcore"
# V = df[Vcol]
# w = df[wcol]

# V_train, V_test, w_train, w_test = train_test_split(V, w, test_size=0.1)
# VWmodel = RandomForestRegressor(criterion="poisson", n_jobs=N_CORES, max_features="log2", oob_score=True)
# weirdrf = VWmodel.fit(V_train, w_train)
# # weirdimp = permutation_importance(weirdrf, V_train, w_train)
# weirdimp = permutation_importance(weirdrf, V_test, w_test)
# pd.Series(weirdimp.importances_mean, index = ["k", "knn", "lcl", "cc", "bc", "Freq"])
# round(weirdrf.oob_score_, 3)
# weirdrf.score(V_test, w_test)

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

# from sklearn.pipeline import make_pipeline # not using, not sure why
