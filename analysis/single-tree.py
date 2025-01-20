# Need a deep tree (> 5 deep) to get reasonable R^2
# Not easy to figure out what's going on, but the first splits are on k-core, eliminating the worse node sets
# having trouble figuring out how to learn more from the data structures


import time
import joblib
N_CORES = joblib.cpu_count(only_physical_cores=True)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import OrdinalEncoder # ordinal encoder fine for regression trees, apparently
from sklearn.model_selection import train_test_split # need to check performance somehow?
# from sklearn.ensemble import RandomForestRegressor # the classifier itself
from sklearn import tree
# from sklearn.tree import DecisionTreeRegressor
from sklearn.inspection import permutation_importance # skl's preferred variable importance measure

networks = [
    "dolphin", "proximity", "celegans", "euroroad", "email", "gkk", "ba", "hk", "lfr", "er",
    "drosophila", "reactome", "route_views", "spanish", "foldoc", "tree_of_life",  "word_assoc",
    "enron", "marker_cafe", "prosper"
]
dynamics = ["doublewell", "mutualistic", "SIS", "genereg"]

csvfile = "../data/nodeset-features.csv"
df = pd.read_csv(csvfile)
df = df.loc[df["network"] != "er"]
df["logerror"] = np.log(df["error"])

encoder = OrdinalEncoder().set_output(transform="pandas")
df["dyn"] = encoder.fit_transform(df[["dynamics"]])
df["net"] = encoder.fit_transform(df[["network"]])
df["ns"] = encoder.fit_transform(df[["ns.type"]])

ycol = "logerror"

# for dyn in dynamics:
net = "email"
dyn = "doublewell"
sdf = df.loc[(df["dynamics"] == dyn) & (df["network"] == net)]
# sdf = df.loc[df["dynamics"] == dyn]
info = {"dynamics": dyn}
xcol = ["k", "knn", "lcl", "cc", "bc", "kcore", "pairs", "geods"] # "net", # , "kcore" is highly correlated w/ k
X = sdf[xcol]
y = sdf[ycol]

# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.1)#, random_state=123)
model = tree.DecisionTreeRegressor(criterion="squared_error", max_depth=5)

# might want a cross-validation step for this one
res = model.fit(X, y)
# tree.plot_tree(res, filled=True, feature_names=xcol, fontsize=8)
# plt.show()

thetree = model.tree_
children_left = thetree.children_left
children_right = thetree.children_right
feature = thetree.feature
threshold = thetree.threshold
values = thetree.value

print(f"Mean degree of optimized node sets is {round(sdf.loc[sdf['ns.type'] == 'opt']['k'].mean(), 3)}")
print(f"Degree-based thresholds in regression tree are:")
print(np.round([t for t, f in zip(threshold, feature) if f == 0], 3))
print(f"Closeness-based thresholds are:")
print(np.round([t for t, f in zip(threshold, feature) if f == 3], 3))
print(f"kcore-based thresholds are:")
print(np.round([t for t, f in zip(threshold, feature) if f == 5], 3))

sdf.groupby("ns.type")["k"].mean()
sdf.groupby("ns.type")["k"].std()
sdf.groupby("ns.type")["cc"].mean()
sdf.groupby("ns.type")["cc"].std()
sdf.groupby("ns.type")["kcore"].mean()
sdf.groupby("ns.type")["kcore"].std()

print(pd.Series([xcol[f] for f in feature if f >= 0]).value_counts())

fi = res.feature_importances_
print(pd.Series({xcol[i]: fi[i] for i in range(len(fi))}))

results = []
for net in networks:
    if net == "er":
        continue
    
    sdf = df.loc[(df["dynamics"] == dyn) & (df["network"] == net)]
    xcol = ["k", "knn", "lcl", "cc", "bc", "kcore", "pairs", "geods"]
    X = sdf[xcol]
    y = sdf[ycol]
    res = model.fit(X, y)
    fi = res.feature_importances_
    results.append({xcol[i]: fi[i] for i in range(len(fi))})

idf = pd.DataFrame(results)

# Across all networks (except er), k, cc, and kcore are the most important.
# For k, the decision tree is clearly looking for the mean degree
# What about for kcore and cc?
