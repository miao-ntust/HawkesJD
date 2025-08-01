clear

lambda = 5;
a = 0.3;
b = 5;
ab = a * b;
r = lambda * (1 - 1 / ab);
t = 1;
lambda0 = 2;

N = 20;

[P, Pc] = hawkesPMF(r, a, b, t, lambda0, N);