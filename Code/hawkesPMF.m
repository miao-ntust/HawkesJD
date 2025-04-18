function [probUC, probC] = hawkesPMF(ga, alpha, beta, t, lambda0, N)
% Hawkes Nt PMF (Condition & Uncondition)

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


w = [1 / beta; -round(2 .^ (n - 1) .* gamma(n - 0.5) / sqrt(pi)) .* (2 * beta) .^ (n - 1) .* alpha .^ n ./ (alpha * beta + 1) .^ (2 * n - 1)];
ww = -w;
ww(1) = -alpha;

B2w = bellPoly(2 * w(2 : end));
Bw = bellPoly(w(2 : end));
Bww = bellPoly(ww(2 : end));

www = [beta / (alpha * beta + 1); accumarray(nn, (-1) .^ kk .* gamma(kk + 1) .* (beta / (alpha * beta + 1)) .^ (kk + 1) .* B2w(nkB))];

O_w_ww = zeros(N, 1);
O_w_ww(1) = 1;
O_w_ww(2 : end) = beta * w(2 : (end - 1)) / (alpha * beta + 1) + accumarray(nnn, Binom .* (-1) .^ jjj .* gamma(jjj + 1) .* beta .^ (jjj + 1) .* (w(nnn - kkk + 1) + alpha * (kkk == nnn)) ./ (alpha * beta + 1) .^ (jjj + 1) .* B2w(kjB));

O_ww_w = -O_w_ww;
O_ww_w(1) = 0;

O_w_0 = zeros(N, 1);
O_w_0(1) = alpha * beta + 1;
O_w_0(2 : end) = beta * w(2 : (end - 1)) + accumarray(nnn, Binom .* (-1) .^ jjj .* gamma(jjj + 1) .* beta .^ (jjj + 1) .* (w(nnn - kkk + 1) + alpha * (kkk == nnn)) .* Bw(kjB));

O_ww_0 = zeros(N, 1);
O_ww_0(2 : end) = -ww(2 : (end - 1)) / alpha - accumarray(nnn, Binom .* gamma(jjj + 1) .* (ww(nnn - kkk + 1) + alpha * (kkk == nnn)) ./ alpha .^ (jjj + 1) .* Bww(kjB));

O_w_L = zeros(N, 1);
O_w_L(1) = exp(beta * t) * (alpha * beta + 1);

O_ww_L = zeros(N, 1);

x = zeros(N + 1, 1);
x(1) = beta * t;

xx = zeros(N + 1, 1);
xx(1) = log(alpha* beta / (alpha * beta + 1 - exp(-beta * t)));

y = zeros(N, 1);
y(1) = -beta * t;

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
    
    O_w_L(in + 1) = beta * exp(beta * t) * w(in + 1) + sum(Binom .* (-1) .^ jjj .* gamma(jjj + 1) .* (beta * exp(beta * t)) .^ (jjj + 1) .* (w(in - kkk + 1) + alpha * (kkk == in)) .* BwL(kjB), 1);
    O_ww_L(in + 1) = -beta * ww(in + 1) / (alpha * beta + 1 - exp(-beta * t)) - sum(Binom .* gamma(jjj + 1) .* beta .^ (jjj + 1) .* (ww(in - kkk + 1) + alpha * (kkk == in)) ./ (alpha * beta + 1 - exp(-beta * t)) .^ (jjj + 1) .* BwwL(kjB), 1);
    
    k = 1 : in;
    
    x(in + 1) = sum((-1) .^ (k - 1) .* gamma(k) .* beta .^ k .* (Bw(in, k) - exp(k * beta * t) .* BwL(in, k)), 2);
    xx(in + 1) = -sum(gamma(k) .* (Bww(in, k) ./ alpha .^ k - beta .^ k .* BwwL(in, k) ./ (alpha * beta + 1 - exp(-beta * t)) .^ k), 2);
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
x(end) =  sum((-1) .^ (k - 1) .* gamma(k) .* beta .^ k .* (Bw(N, k) - exp(k * beta * t) .* BwL(N, k)), 2);
xx(end) = -sum(gamma(k) .* (Bww(N, k) ./ alpha .^ k - beta .^ k .* BwwL(N, k) ./ (alpha * beta + 1 - exp(-beta * t)) .^ k), 2);

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
theta = (-1) .^ (k2 == 0) * alpha .* (k2 <= 1);
H = ga * (L - 1 / beta * accumarray(nk, BinomH .* www(k1 + 1) .* ((w(k2 + 1) - theta) .* x(k3 + 1) - (ww(k2 + 1) - theta) .* xx(k3 + 1))));

Mc = H - L * lambda0;

BM = bellPoly(beta * L / (alpha * beta - 1));
M = H - ga * (L + 1 / beta * accumarray(nn, (-1) .^ (kk - 1) .* gamma(kk) ./ ((1 - exp(-beta * t)) / (alpha * beta - 1) + 1) .^ kk .* BM(nkB)));

probC = zeros(N + 1, 1);
probC(1) = exp((ga - lambda0) / beta * (1 - exp(-beta * t)) - ga * t);
probC(2 : end) = 1 ./ gamma(n + 1) * probC(1) .* sum(bellPoly(Mc), 2);

probUC = zeros(N + 1, 1);
probUC(1) = exp(-ga * t) * ((1 - exp(-beta * t)) / (alpha * beta - 1) + 1) ^ (-ga / beta);
probUC(2 : end) = 1 ./ gamma(n + 1) * probUC(1) .* sum(bellPoly(M), 2);

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
