function P = poisPMF(lambda, t, N)

n = (0 : N)';
P = exp(n * log(lambda * t) - lambda * t - gammaln(n + 1));

end