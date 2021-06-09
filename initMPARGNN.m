function filter = initMPARGNN(detection)

filter = initcvekf(detection);
filter.StateCovariance(2,2) = 1e6;
filter.StateCovariance(4,4) = 1e6;
filter.StateCovariance(6,6) = 1e6;