function mvsk = mvskNt(PMF)
% mvsk of Nt
%
% Input:
%   PMF: column vector
%
% Output:
%   mvsk: row vector

N = numel(PMF);
ER = sum(PMF .* ((0 : (N - 1))' .^ (1 : 4)), 1);

mvsk = zeros(1, 4);
mvsk(1) = ER(1);
mvsk(2) = ER(2) - ER(1) .^ 2;
mvsk(3) = (ER(3) - 3 * ER(1) .* ER(2) + 2 * ER(1) .^ 3) ./ mvsk(2) .^ 1.5;
mvsk(4) = (ER(4) - 4 * ER(1) .* ER(3) + 6 * ER(1) .^ 2 .* ER(2) - 3 * ER(1) .^ 4) ./ mvsk(2) .^ 2;

end

