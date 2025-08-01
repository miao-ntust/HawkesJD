clear

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
sheet = 'MVSK';
xlRange = {'E5', 'E19', 'E33', 'E47'};

N = 100;

PP = poisPMF(lambda, t, N);
PPN2 = 1 - exp(-lambda * t) * sum((lambda * t) .^ (0 : 10) ./ exp(gammaln((0 : 10) + 1)), 2);
PPN3 = 1 - exp(-lambda * t) * sum((lambda * t) .^ (0 : 20) ./ exp(gammaln((0 : 20) + 1)), 2);
[~, PPR2] = returnMJD(-0.3, r, q, sigma, t, gamma, delta, PP);
[~, PPR3] = returnMJD(-0.5, r, q, sigma, t, gamma, delta, PP);

for i = 1 : 3
    tableX = cell(10, 15);
    mvskN = zeros(10, 4);
    PN = zeros(10, 2);
    mvskR = zeros(10, 4);
    PR = zeros(10, 2);
    mvskN(1, :) = mvskNt(PP);
    PN(1, 1) = PPN2;
    PN(1, 2) = PPN3;
    mvskR(1, :) = mvskMJD(r, q, sigma, t, gamma, delta, PP);
    PR(1, 1) = PPR2;
    PR(1, 2) = PPR3;
    PP0 = poisPMF(lambda0(i), t, N);
    mvskN(2, :) = mvskNt(PP0);
    PN(2, 1) = 1 - exp(-lambda0(i) * t) * sum((lambda0(i) * t) .^ (0 : 10) ./ exp(gammaln((0 : 10) + 1)), 2);
    PN(2, 2) = 1 - exp(-lambda0(i) * t) * sum((lambda0(i) * t) .^ (0 : 20) ./ exp(gammaln((0 : 20) + 1)), 2);
    mvskR(2, :) = mvskMJD(r, q, sigma, t, gamma, delta, PP0);
    [~, PR(2, 1)] = returnMJD(-0.3, r, q, sigma, t, gamma, delta, PP0);
    [~, PR(2, 2)] = returnMJD(-0.5, r, q, sigma, t, gamma, delta, PP0);
    for ia = 1 : 8
        str = ['simData_C', num2str(lambda0(i)), '_', num2str(ia), '.mat'];
        load(str);
        PH = P(:, 3);
        mvskN(ia + 2, :) = mvskNt(PH);
        PN(ia + 2, 1) = sum(PH(12 : end), 1);
        PN(ia + 2, 2) = sum(PH(22 : end), 1);
        mvskR(ia + 2, :) = mvskMJD(r, q, sigma, t, gamma, delta, PH);
        [~, PR(ia + 2, 1)] = returnMJD(-0.3, r, q, sigma, t, gamma, delta, PH);
        [~, PR(ia + 2, 2)] = returnMJD(-0.5, r, q, sigma, t, gamma, delta, PH);
    end
    
    for j = 1 : 4
        tableX{1, j} = ['$', num2str(mvskN(1, j), '%.3f'), '\phantom{)}$'];
        tableX{2, j} = ['$(', num2str(mvskN(2, j), '%.3f'), ')$'];
        tableX{1, 8 + j} = ['$', num2str(mvskR(1, j), '%.3f'), '\phantom{)}$'];
        tableX{2, 8 + j} = ['$(', num2str(mvskR(2, j), '%.3f'), ')$'];
        for ia = 1 : 8
            tableX{ia + 2, j} =['$', num2str(mvskN(ia + 2, j), '%.3f'), '\phantom{)}$'];
            tableX{ia + 2, 8 + j} =['$', num2str(mvskR(ia + 2, j), '%.3f'), '\phantom{)}$'];
        end
    end

    for j = 1 : 2
        tableX{1, 5 + j} = ['$', num2str(PN(1, j) * 100, '%.3f'), '\phantom{)}%$'];
        tableX{2, 5 + j} = ['$(', num2str(PN(2, j) * 100, '%.3f'), ')%$'];
        tableX{1, 13 + j} = ['$', num2str(PR(1, j), '%.3f'), '\phantom{)}$'];
        tableX{2, 13 + j} = ['$(', num2str(PR(2, j), '%.3f'), ')$'];
        for ia = 1 : 8
            tableX{ia + 2, 5 + j} =['$', num2str(PN(ia + 2, j) * 100, '%.3f'), '\phantom{)}%$'];
            tableX{ia + 2, 13 + j} =['$', num2str(PR(ia + 2, j), '%.3f'), '\phantom{)}$'];
        end
    end
    xlswrite(filename, tableX, sheet, xlRange{i});
end

tableX = cell(10, 15);
mvskN = zeros(10, 4);
PN = zeros(10, 2);
mvskR = zeros(10, 4);
PR = zeros(10, 2);
mvskN(1, :) = mvskNt(PP);
PN(1, 1) = PPN2;
PN(1, 2) = PPN3;
mvskR(1, :) = mvskMJD(r, q, sigma, t, gamma, delta, PP);
PR(1, 1) = PPR2;
PR(1, 2) = PPR3;
sigBS = sqrt((gamma ^ 2 + delta ^ 2) * lambda + sigma ^ 2);
muBS = r - 0.5 * sigBS ^ 2;
mvskR(2, :) = [muBS * t, sigBS ^ 2 * t, 0, 3];
PR(2, 1) = 0.5 * erfc(-(-0.3 - muBS * t) / (sigBS * sqrt(2 * t)));
PR(2, 2) = 0.5 * erfc(-(-0.5 - muBS * t) / (sigBS * sqrt(2 * t)));
for ia = 1 : 8
    str = ['simData_U_', num2str(ia), '.mat'];
    load(str);
    PH = P(:, 3);
    mvskN(ia + 2, :) = mvskNt(PH);
    PN(ia + 2, 1) = sum(PH(12 : end), 1);
    PN(ia + 2, 2) = sum(PH(22 : end), 1);
    mvskR(ia + 2, :) = mvskMJD(r, q, sigma, t, gamma, delta, PH);
    [~, PR(ia + 2, 1)] = returnMJD(-0.3, r, q, sigma, t, gamma, delta, PH);
    [~, PR(ia + 2, 2)] = returnMJD(-0.5, r, q, sigma, t, gamma, delta, PH);
end

for j = 1 : 4
    tableX{1, j} = ['$', num2str(mvskN(1, j), '%.3f'), '\phantom{)}$'];
    tableX{1, 8 + j} = ['$', num2str(mvskR(1, j), '%.3f'), '\phantom{)}$'];
    tableX{2, 8 + j} = ['$(', num2str(mvskR(2, j), '%.3f'), ')$'];
    for ia = 1 : 8
        tableX{ia + 2, j} =['$', num2str(mvskN(ia + 2, j), '%.3f'), '\phantom{)}$'];
        tableX{ia + 2, 8 + j} =['$', num2str(mvskR(ia + 2, j), '%.3f'), '\phantom{)}$'];
    end
end

for j = 1 : 2
    tableX{1, 5 + j} = ['$', num2str(PN(1, j) * 100, '%.3f'), '\phantom{)}%$'];
    tableX{1, 13 + j} = ['$', num2str(PR(1, j), '%.3f'), '\phantom{)}$'];
    tableX{2, 13 + j} = ['$(', num2str(PR(2, j), '%.3f'), ')$'];
    for ia = 1 : 8
        tableX{ia + 2, 5 + j} =['$', num2str(PN(ia + 2, j) * 100, '%.3f'), '\phantom{)}%$'];
        tableX{ia + 2, 13 + j} =['$', num2str(PR(ia + 2, j), '%.3f'), '\phantom{)}$'];
    end
end
xlswrite(filename, tableX, sheet, xlRange{4});

