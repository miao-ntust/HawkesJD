% JD params
S0    = 100;
K     = [90 100 110];   % row vector
r     = 0.05;                  % risk-free
q     = 0.00;                  % div. yield
sigma = 0.20;
t     = 1;                  % maturity (years)
theta = -0.05;                 % mean jump size
delta = 0.10;                  % jump size vol

% Hawkes PMF params
rH      = 2.50;
a       = 1.00;
b       = 2.00;
lambda0 = 2.00;
N       = 50;

% get PMFs from Hawkes
[P, Pc] = hawkesPMF(rH, a, b, t, lambda0, N);

% use Pc as jumpProb
jumpProb = Pc(:);
pricePc  = callJD(S0, K, r, q, sigma, t, theta, delta, jumpProb);

% use P as jumpProb
jumpProb = P(:);
priceP   = callJD(S0, K, r, q, sigma, t, theta, delta, jumpProb);

% show results
format long g
disp('Prices using Pc:'), disp(pricePc)
disp('Prices using P:'),  disp(priceP)

