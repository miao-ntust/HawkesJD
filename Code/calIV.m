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
