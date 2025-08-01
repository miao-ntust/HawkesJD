clear

load('simData_C8_1.mat')
PSimC = P(1 : 21, 3);
load('simData_U_1.mat')
PSimU = P(1 : 21, 3);

lambda = 5;
a = 0.3;
b = 5;
ab  = a * b;
v = lambda * (1 - 1 / ab);
t = 1;
lambda0 = 8;

filename = 'Table.xlsx';
sheet = 'Simulation';
xlRange = 'D4';

N = 20;
[PAnU, PAnC] = hawkesPMF(v, a, b, t, lambda0, N);

SEC = sqrt(PSimC .* (1 - PSimC) / 10000000);
SEU = sqrt(PSimU .* (1 - PSimU) / 10000000);
CIC = [PSimC - 1.96 * SEC, PSimC + 1.96 * SEC];
CIU = [PSimU - 1.96 * SEU, PSimU + 1.96 * SEU];

tableX = cell(21, 9);
for i = 1 : 21
    tableX{i, 1} = ['$', num2str(PAnC(i), '%.6f'), '$'];
    tableX{i, 2} = ['$', num2str(PSimC(i), '%.6f'), '$'];
    tableX{i, 3} = ['$', num2str(SEC(i), '%.6f'), '$'];
    tableX{i, 4} = ['[$', num2str(CIC(i, 1), '%.6f'), ', ~', num2str(CIC(i, 2), '%.6f'), '$]'];
    tableX{i, 6} = ['$', num2str(PAnU(i), '%.6f'), '$'];
    tableX{i, 7} = ['$', num2str(PSimU(i), '%.6f'), '$'];
    tableX{i, 8} = ['$', num2str(SEU(i), '%.6f'), '$'];
    tableX{i, 9} = ['[$', num2str(CIU(i, 1), '%.6f'), ', ~', num2str(CIU(i, 2), '%.6f'), '$]'];
end

xlswrite(filename, tableX, sheet, xlRange);
