load("proximity-data.rda")

plot(data.frame(A.e, k, bc, knn, lcl), log = "xy")
plot(data.frame(B.e, k, bc, knn, lcl), log = "xy")

plot(A.e, B.e, log = "xy")


idx <- which(k > 2 & k < 11)
plot(data.frame(A.e, k, bc, knn, lcl)[idx, ], log = "xy")
plot(data.frame(B.e, k, bc, knn, lcl)[idx, ], log = "xy")
