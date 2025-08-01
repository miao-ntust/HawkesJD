clear

lambda = 5;

a = [0.3, 0.15, 0.4, 0.2, 0.5, 0.25, 0.6, 0.3];
b = [5, 10, 5, 10, 5, 10, 5, 10];
ab = a .* b;
r = lambda * (1 - 1 ./ ab);
t = [0.25, 0.5, 1];

lambda0 = [2, 5, 8];

paths = 1000;

for i = 2 : 2
    P = hawkesSim(r(i), a(i), b(i), t, paths);
    str = ['simData_U_', num2str(i), '.mat'];
    save(str, 'P');
    for j = 1 : 3
        P = hawkesSimC(r(i), lambda0(j), a(i), b(i), t, paths);
        str = ['simData_C', num2str(lambda0(j)), '_', num2str(i), '.mat'];
        save(str, 'P');
    end
end
