% Hawkes PMF params
r       = 2.5;   
a       = 1;     
b       = 2;     
t       = 1;    
lambda0 = 8;     % initial intensity
N       = 20;    % print probabilities for counts 0..N

% get PMFs from Hawkes
[P, Pc] = hawkesPMF(r, a, b, t, lambda0, N);

% show results
format long g
n = (0:N)';

fprintf('\nPc (N_t = n):\n');
for i = 1:numel(n)
    fprintf('n=%2d\tPc=%.12g\n', n(i), Pc(i));
end

fprintf('\nP (N_t = n):\n');
for i = 1:numel(n)
    fprintf('n=%2d\tP=%.12g\n', n(i), P(i));
end

