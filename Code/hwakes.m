function [PMF, Tail, Price, MVSK] = hwakes(S0, K, r, q, sigma, t, gamma, delta, lambda, v, a, b, N)

na = numel(a);

PMF = zeros(N + 1, na + 1);
PMF(:, 1) = poisPMF(lambda, t, N);

Price = cell(10, 18);
call = zeros(na + 2, 3);
iv = zeros(na + 2, 3);
call(1, :) = callMJD(S0, K, r, q, sigma, t, gamma, delta, PMF(:, 1));
iv(1, :) = calIV(call(1, :), S0, K, r, q, t);
sigBS = sqrt((gamma ^ 2 + delta ^ 2) * lambda + sigma ^ 2);
call(2, :) = callBS(S0, K, r, q, sigBS, t);
iv(2, :) = sigBS;

MVSK = zeros(na + 1, 4);
MVSK(1, :) = mvskMJD(r, q, sigma, t, 0, delta, PMF(:, 1));
for ia = 1 : na
    PMF(:, ia + 1) = hawkesPMF(v(ia), a(ia), b(ia), t, N);
    call(ia + 2, :)= callMJD(S0, K, r, q, sigma, t, gamma, delta, PMF(:, ia + 1));
    iv(ia + 2, :) = calIV(call(ia + 2, :), S0, K, r, q, t);
    MVSK(ia + 1, :) = mvskMJD(r, q, sigma, t, 0, delta, PMF(:, ia + 1));
end
for iK = 1 : 3
    Price{1, 6 * iK - 5} = num2str(call(1, iK), '%.3f');
    Price{1, 6 * iK - 2} = ['$', num2str(iv(1, iK) * 100, '%.3f'), '%$'];
    Price{2, 6 * iK - 5} = ['(', num2str(call(2, iK), '%.3f'), ')'];
    Price{2, 6 * iK - 2} = ['($', num2str(sigBS * 100, '%.3f'), '%$)'];
    for ia = 1 : 8
        Price{ia + 2, 6 * iK - 5} =num2str(call(ia + 2, iK), '%.3f');
        Price{ia + 2, 6 * iK - 4} =['$', num2str((call(ia + 2, iK) - call(1, iK)) ./ call(1, iK) * 100, '%.3f'), '%$'];
        Price{ia + 2, 6 * iK - 2} =['$', num2str(iv(ia + 2, iK) * 100, '%.3f'), '%$'];
        Price{ia + 2, 6 * iK - 1} =['$', num2str((iv(ia + 2, iK) - iv(1, iK)) ./ iv(1, iK) * 100, '%.3f'), '%$'];
    end
end

Tail = 1 - sum(PMF, 1);

end

function P = hawkesPMF(r, a, b, t, N)

n = (1 : N)';

[nn, kk] = ndgrid(1 : N);
nn = reshape(nn, [], 1);
kk = reshape(kk, [], 1);
ind = find(kk > nn);
nn(ind) = [];
kk(ind) = [];
nkB = N * (kk - 1) + nn;

[nnn, kkk, jjj] = ndgrid(1 : (N - 1));
nnn = reshape(nnn, [], 1);
kkk = reshape(kkk, [], 1);
jjj = reshape(jjj, [], 1);
ind = find((kkk > nnn) | (jjj > kkk));
nnn(ind) = [];
kkk(ind) = [];
jjj(ind) = [];
Binom = gamma(nnn + 1) ./ gamma(kkk + 1) ./ gamma(nnn - kkk + 1);
kjB = N * (jjj - 1) + kkk;


w = [1 / b; -round(2 .^ (n - 1) .* gamma(n - 0.5) / sqrt(pi)) .* (2 * b) .^ (n - 1) .* a .^ n ./ (a * b + 1) .^ (2 * n - 1)];
ww = -w;
ww(1) = -a;

B2w = bellPoly(2 * w(2 : end));
Bw = bellPoly(w(2 : end));
Bww = bellPoly(ww(2 : end));

www = [b / (a * b + 1); accumarray(nn, (-1) .^ kk .* gamma(kk + 1) .* (b / (a * b + 1)) .^ (kk + 1) .* B2w(nkB))];

O_w_ww = zeros(N, 1);
O_w_ww(1) = 1;
O_w_ww(2 : end) = b * w(2 : (end - 1)) / (a * b + 1) + accumarray(nnn, Binom .* (-1) .^ jjj .* gamma(jjj + 1) .* b .^ (jjj + 1) .* (w(nnn - kkk + 1) + a * (kkk == nnn)) ./ (a * b + 1) .^ (jjj + 1) .* B2w(kjB));

O_ww_w = -O_w_ww;
O_ww_w(1) = 0;

O_w_0 = zeros(N, 1);
O_w_0(1) = a * b + 1;
O_w_0(2 : end) = b * w(2 : (end - 1)) + accumarray(nnn, Binom .* (-1) .^ jjj .* gamma(jjj + 1) .* b .^ (jjj + 1) .* (w(nnn - kkk + 1) + a * (kkk == nnn)) .* Bw(kjB));

O_ww_0 = zeros(N, 1);
O_ww_0(2 : end) = -ww(2 : (end - 1)) / a - accumarray(nnn, Binom .* gamma(jjj + 1) .* (ww(nnn - kkk + 1) + a * (kkk == nnn)) ./ a .^ (jjj + 1) .* Bww(kjB));

O_w_L = zeros(N, 1);
O_w_L(1) = exp(b * t) * (a * b + 1);

O_ww_L = zeros(N, 1);

x = zeros(N + 1, 1);
x(1) = b * t;

xx = zeros(N + 1, 1);
xx(1) = log(a* b / (a * b + 1 - exp(-b * t)));

y = zeros(N, 1);
y(1) = -b * t;

yy = zeros(N + 1, 1);
yy(1) = xx(1);

Psi = zeros(N, 1);
Psi(1) = -w(2) * (O_w_0(1) + O_ww_0(1) - O_w_L(1) - O_ww_L(1) + y(1) + yy(1));

Phi = zeros(N, 1);
Phi(1) = O_w_L(1) - O_ww_L(1);

L = zeros(N, 1);
L(1) = Psi(1) / Phi(1); 

for in = 1 : (N - 1)
    [kkk, jjj] = ndgrid(1 : in);
    kkk = reshape(kkk, [], 1);
    jjj = reshape(jjj, [], 1);
    ind = find(jjj > kkk);
    kkk(ind) = [];
    jjj(ind) = [];
    Binom = gamma(in + 1) ./ gamma(kkk + 1) ./ gamma(in - kkk + 1);
    kjB = in * (jjj - 1) + kkk;
    
    BwL = bellPoly(w(2 : (in + 1)) - L(1 : in));
    BwwL = bellPoly(ww(2 : (in + 1)) - L(1 : in));
    
    O_w_L(in + 1) = b * exp(b * t) * w(in + 1) + sum(Binom .* (-1) .^ jjj .* gamma(jjj + 1) .* (b * exp(b * t)) .^ (jjj + 1) .* (w(in - kkk + 1) + a * (kkk == in)) .* BwL(kjB), 1);
    O_ww_L(in + 1) = -b * ww(in + 1) / (a * b + 1 - exp(-b * t)) - sum(Binom .* gamma(jjj + 1) .* b .^ (jjj + 1) .* (ww(in - kkk + 1) + a * (kkk == in)) ./ (a * b + 1 - exp(-b * t)) .^ (jjj + 1) .* BwwL(kjB), 1);
    
    k = 1 : in;
    
    x(in + 1) = sum((-1) .^ (k - 1) .* gamma(k) .* b .^ k .* (Bw(in, k) - exp(k * b * t) .* BwL(in, k)), 2);
    xx(in + 1) = -sum(gamma(k) .* (Bww(in, k) ./ a .^ k - b .^ k .* BwwL(in, k) ./ (a * b + 1 - exp(-b * t)) .^ k), 2);
    Binomy = gamma(in + 1) ./ gamma(k + 1) ./ gamma(in - k + 1);
    y(in + 1) = -x(in + 1) - 2 * sum(Binomy .* O_w_ww(k + 1)' .* x(in - k + 1)', 2);
    yy(in + 1) = xx(in + 1) - 2 * sum(Binomy .* O_ww_w(k + 1)' .* xx(in - k + 1)', 2);

    k = (0 : in)';
    Psi(in + 1) = sum(gamma(in + 1) ./ gamma(k + 1) ./ gamma(in - k + 1) .* -w(in - k + 2) .* (O_w_0(k + 1) + O_ww_0(k + 1) - O_w_L(k + 1) - O_ww_L(k + 1) + y(k + 1) + yy(k + 1)), 1);
    Phi(in + 1) = O_w_L(in + 1) - O_ww_L(in + 1);
    
    BPhi = bellPoly(Phi(2 : (in + 1)));
    L(in + 1) = Psi(in + 1) / Phi(1) + sum(Binom .* (-1) .^ jjj .* gamma(jjj + 1) .* Psi(in - kkk + 1) ./ Phi(1) .^ (jjj + 1) .* BPhi(kjB), 1);
end

k = 1 : N;
BwL = bellPoly(w(2 : end) - L);
BwwL = bellPoly(ww(2 : end) - L);
x(end) =  sum((-1) .^ (k - 1) .* gamma(k) .* b .^ k .* (Bw(N, k) - exp(k * b * t) .* BwL(N, k)), 2);
xx(end) = -sum(gamma(k) .* (Bww(N, k) ./ a .^ k - b .^ k .* BwwL(N, k) ./ (a * b + 1 - exp(-b * t)) .^ k), 2);

[k1, k2, k3] = ndgrid(0 : N);
k1 = reshape(k1, [], 1);
k2 = reshape(k2, [], 1);
k3 = reshape(k3, [], 1);
ind = find((k1 + k2 + k3 > N) | (k1 + k2 + k3 == 0));
k1(ind) = [];
k2(ind) = [];
k3(ind) = [];
nk = k1 + k2 + k3;
BinomH = gamma(nk + 1) ./ gamma(k1 + 1) ./ gamma(k2 + 1) ./ gamma(k3 + 1);
theta = (-1) .^ (k2 == 0) * a .* (k2 <= 1);
H = r * (L - 1 / b * accumarray(nk, BinomH .* www(k1 + 1) .* ((w(k2 + 1) - theta) .* x(k3 + 1) - (ww(k2 + 1) - theta) .* xx(k3 + 1))));

BM = bellPoly(b * L / (a * b - 1));
M = H - r * (L + 1 / b * accumarray(nn, (-1) .^ (kk - 1) .* gamma(kk) ./ ((1 - exp(-b * t)) / (a * b - 1) + 1) .^ kk .* BM(nkB)));

P = zeros(N + 1, 1);
P(1) = exp(-r * t) * ((1 - exp(-b * t)) / (a * b - 1) + 1) ^ (-r / b);
P(2 : end) = 1 ./ gamma(n + 1) * P(1) .* sum(bellPoly(M), 2);

end

function B = bellPoly(x)

n = numel(x);

B = zeros(n + 1);
B(1) = 1;

for in = 1 : n
    for ik = 1 : in
        i = 1 : (in - ik + 1);
        B(in + 1, ik + 1) = sum(gamma(in) ./ gamma(i)' ./ gamma(in - i + 1)' .* x(i) .* B(in - i + 1, ik), 1);
    end
end

B(:, 1) = [];
B(1, :) = [];

end

function P = poisPMF(lambda, t, N)

n = (0 : N)';
P = exp(n * log(lambda * t) - lambda * t - gammaln(n + 1));

end

function price = callMJD(S0, K, r, q, sigma, t, gamma, delta, jumpProb)
% Call price of JD
%
% Input:
%   K: row vector
%   jumpProb: column vector
%
% Output: price(n, K)

N = numel(jumpProb) - 1;
n = (0 : N)';

muX = n * gamma;
varX = n * delta .^ 2;

eta = log(sum(jumpProb .* exp(muX + varX / 2), 1));
mu = r - q - eta / t;
muD = (mu - sigma ^ 2 / 2) * t;
varD = sigma ^ 2 * t;

muR = muD + muX;
varR = varD + varX;

sigmahat = sqrt(varR / t);
rhat = muR / t + q + sigmahat .^ 2 / 2;

price = sum(jumpProb .* exp(-(r - rhat) * t) .* callBS(S0, K, rhat, q, sigmahat, t), 1);

end

function iv = calIV(price, S0, K, r, q, t)

S = S0 * exp(-q * t);
X = K * exp(-r * t);

a = (S + X) * t;
b = sqrt(8 * pi * t) * ((S - X) / 2 - price);
c = 2 * (S - X) .* log(S ./ X);
temp = max(b .^ 2 - 4 * a .* c, 0);
iv0 = (-b + sqrt(temp)) ./ (2 * a);

iv = iv0;
f = @(x) callBS(S0, K, r, q, x, t) - price;
fd = @(x) vegaBS(S0, K, r, q, x, t);
while max(abs(f(iv))) > 1e-12
    iv0 = iv;
	iv = iv0 - f(iv0) ./ fd(iv0);
end

end

function call = callBS(S0, K, r, q, sigma, t)

d1 = (log(S0 ./ K) + (r - q + sigma .^ 2 / 2) .* t) ./ (sigma .* sqrt(t));
d2 = d1 - sigma .* sqrt(t);
nd1 = 0.5 * erfc(-d1 / sqrt(2));
nd2 = 0.5 * erfc(-d2 / sqrt(2));
call = S0 .* exp(-q .* t) .* nd1 - K .* exp(-r .* t) .* nd2;

end

function vega = vegaBS(S0, K, r, q, sigma, t)

d1 = (log(S0 ./ K) + (r - q + sigma .^ 2 / 2) .* t) ./ (sigma .* sqrt(t));
vega = S0 * exp(-q * t - d1 .^ 2 / 2) .* sqrt(t / (2 * pi));

end

function mvsk = mvskMJD(r, q, sigma, t, gamma, delta, jumpProb)
% Return mvsk of JD
%
% Input:
%   jumpProb: column vector
%
% Output:
%   mvsk: row vector

N = numel(jumpProb) - 1;
n = (0 : N)';

muX = n * gamma;
varX = n * delta .^ 2;

eta = log(sum(jumpProb .* exp(muX + varX / 2), 1));
mu = r - q - eta / t;
muD = (mu - sigma ^ 2 / 2) * t;
varD = sigma ^ 2 * t;

muR = muD + muX;
varR = varD + varX;

ER1 = sum(jumpProb .* muR, 1);
ER2 = sum(jumpProb .* (muR .^ 2 + varR), 1);
ER3 = sum(jumpProb .* (muR .^ 3 + 3 * muR .* varR), 1);
ER4 = sum(jumpProb .* (muR .^ 4 + 6 * muR .^ 2 .* varR + 3 * varR .^ 2), 1);

mvsk = zeros(1, 4);
mvsk(1) = ER1;
mvsk(2) = ER2 - ER1 .^ 2;
mvsk(3) = (ER3 - 3 * ER1 .* ER2 + 2 * ER1 .^ 3) ./ mvsk(2) .^ 1.5;
mvsk(4) = (ER4 - 4 * ER1 .* ER3 + 6 * ER1 .^ 2 .* ER2 - 3 * ER1 .^ 4) ./ mvsk(2) .^ 2;

end
