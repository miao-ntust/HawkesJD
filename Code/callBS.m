function call = callBS(S0, K, r, q, sigma, t)

d1 = (log(S0 ./ K) + (r - q + sigma .^ 2 / 2) .* t) ./ (sigma .* sqrt(t));
d2 = d1 - sigma .* sqrt(t);
nd1 = 0.5 * erfc(-d1 / sqrt(2));
nd2 = 0.5 * erfc(-d2 / sqrt(2));
call = S0 .* exp(-q .* t) .* nd1 - K .* exp(-r .* t) .* nd2;

end

