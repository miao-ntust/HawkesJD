r       = 2.5;
a       = 1;
b       = 2;   
t       = 1;
lambda0 = 8;
N       = 20;

[P, Pc] = hawkesPMF(r, a, b, t, lambda0, N);

format long g
n = (0:N)';

fprintf('\nPc (N_t = n):\n');
for i = 1:numel(n)
    fprintf('n=%2d\tPc=%.12g\n', n(i), Pc(i));
end

fprintf('P (N_t = n):\n');
for i = 1:numel(n)
    fprintf('n=%2d\tP=%.12g\n', n(i), P(i));
end

