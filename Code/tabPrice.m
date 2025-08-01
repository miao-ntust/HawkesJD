clear

S0 = 100;
K = [90, 100, 110];
r = 0.05;
q = 0;
sigma = 0.2;
t = 1;
lambda = 5;
gamma = -0.05;
delta = 0.1;

lambda0 = [2, 5, 8];

a = [0.3, 0.15, 0.4, 0.2, 0.5, 0.25, 0.6, 0.3];
b = [5, 10, 5, 10, 5, 10, 5, 10];

ab  = a .* b;
v = lambda * (1 - 1 ./ ab);

filename = 'Table.xlsx';
sheet = 'Price';
xlRange = {'E6', 'E21', 'E36', 'E51'};

N = 100;

PP = poisPMF(lambda, t, N);

for i = 1 : 3
    tableX = cell(10, 17);
    call = zeros(10, 3);
    iv = zeros(10, 3);
    call(1, :) = callMJD(S0, K, r, q, sigma, t, gamma, delta, PP);
    iv(1, :) = calIV(call(1, :), S0, K, r, q, t);
    PP0 = poisPMF(lambda0(i), t, N);
    call(2, :) = callMJD(S0, K, r, q, sigma, t, gamma, delta, PP0);
    iv(2, :) = calIV(call(2, :), S0, K, r, q, t);
    for ia = 1 : 8
        str = ['simData_C', num2str(lambda0(i)), '_', num2str(ia), '.mat'];
        load(str);
        PH = P(1 : 101, 3);
        call(ia + 2, :)= callMJD(S0, K, r, q, sigma, t, gamma, delta, PH);
        iv(ia + 2, :) = calIV(call(ia + 2, :), S0, K, r, q, t);
    end
    for iK = 1 : 3
        tableX{1, 6 * iK - 5} = ['$', num2str(call(1, iK), '%.3f'), '$'];
        tableX{1, 6 * iK - 2} = ['$', num2str(iv(1, iK) * 100, '%.3f'), '%$'];
        tableX{2, 6 * iK - 5} = ['$(', num2str(call(2, iK), '%.3f'), ')$'];
        tableX{2, 6 * iK - 2} = ['($', num2str(iv(2, iK) * 100, '%.3f'), '%$)'];
        for ia = 1 : 8
            tableX{ia + 2, 6 * iK - 5} =['$', num2str(call(ia + 2, iK), '%.3f'), '$'];
            tableX{ia + 2, 6 * iK - 4} =['$', num2str((call(ia + 2, iK) - call(1, iK)) ./ call(1, iK) * 100, '%.3f'), '%$'];
            tableX{ia + 2, 6 * iK - 2} =['$', num2str(iv(ia + 2, iK) * 100, '%.3f'), '%$'];
            tableX{ia + 2, 6 * iK - 1} =['$', num2str((iv(ia + 2, iK) - iv(1, iK)) ./ iv(1, iK) * 100, '%.3f'), '%$'];
        end
    end
    xlswrite(filename, tableX, sheet, xlRange{i});
end

tableX = cell(10, 17);
call = zeros(10, 3);
iv = zeros(10, 3);
call(1, :) = callMJD(S0, K, r, q, sigma, t, gamma, delta, PP);
iv(1, :) = calIV(call(1, :), S0, K, r, q, t);
sigBS = sqrt((gamma ^ 2 + delta ^ 2) * lambda + sigma ^ 2);
call(2, :) = callBS(S0, K, r, q, sigBS, t);
iv(2, :) = sigBS;
for ia = 1 : 8
    str = ['simData_U_', num2str(ia), '.mat'];
    load(str);
    PH = P(1 : 101, 3);
    call(ia + 2, :)= callMJD(S0, K, r, q, sigma, t, gamma, delta, PH);
    iv(ia + 2, :) = calIV(call(ia + 2, :), S0, K, r, q, t);
end
for iK = 1 : 3
    tableX{1, 6 * iK - 5} = ['$', num2str(call(1, iK), '%.3f'), '$'];
    tableX{1, 6 * iK - 2} = ['$', num2str(iv(1, iK) * 100, '%.3f'), '%$'];
    tableX{2, 6 * iK - 5} = ['$(', num2str(call(2, iK), '%.3f'), ')$'];
    tableX{2, 6 * iK - 2} = ['($', num2str(iv(2, iK) * 100, '%.3f'), '%$)'];
    for ia = 1 : 8
        tableX{ia + 2, 6 * iK - 5} =['$', num2str(call(ia + 2, iK), '%.3f'), '$'];
        tableX{ia + 2, 6 * iK - 4} =['$', num2str((call(ia + 2, iK) - call(1, iK)) ./ call(1, iK) * 100, '%.3f'), '%$'];
        tableX{ia + 2, 6 * iK - 2} =['$', num2str(iv(ia + 2, iK) * 100, '%.3f'), '%$'];
        tableX{ia + 2, 6 * iK - 1} =['$', num2str((iv(ia + 2, iK) - iv(1, iK)) ./ iv(1, iK) * 100, '%.3f'), '%$'];
    end
end
xlswrite(filename, tableX, sheet, xlRange{4});
