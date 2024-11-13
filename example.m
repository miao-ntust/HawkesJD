clear

S0 = 100;
K = 90;
r = 0.05;
q = 0;
sigma = 0.2;
t = 1;
lambda = 5;
theta = -0.05;
delta = 0.1;

alpha = 1.5;
beta = 2;
ga = lambda * (1 - 1 ./ (alpha * beta));
lambda0 = 2;

N = 20;

[~, probC] = hawkesPMF(ga, alpha, beta, t, lambda0, N);
price = callJD(S0, K, r, q, sigma, t, theta, delta, probC);