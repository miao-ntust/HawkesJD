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

