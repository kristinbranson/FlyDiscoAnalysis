function X = ClustergramImputeFun(X)

X(isnan(X)) = 0;