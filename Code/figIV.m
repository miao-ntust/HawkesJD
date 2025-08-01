clear

color = get(groot, 'DefaultAxesColorOrder');

S0 = 100;
nK = 201;
K = linspace(50, 150, nK);
r = 0.05;
q = 0;
sigma = 0.2;
t = 1;
lambda = 5;
gamma = -0.05;
delta = 0.1;

a = [0.3, 0.4, 0.3];
b = [5, 5, 10];
ab  = a .* b;
v = lambda * (1 - 1 ./ ab);

ind = [1, 3, 8];

lambda0 = [2, 5, 8];

N = 199;

for i = 1 : 3
    prob = zeros(N + 1, 5);
    prob(:, 1) = poisPMF(lambda, t, N);
    prob(:, 2) = poisPMF(lambda0(i), t, N);
    for ia = 1 : 3
        str = ['simData_C', num2str(lambda0(i)), '_', num2str(ind(ia)), '.mat'];
        load(str);
        prob(:, ia + 2) = P(:, 3);
    end
    call = zeros(nK, 5);
    iv = zeros(nK, 5);
    for j = 1 : 5
        call(:, j) = callMJD(S0, K, r, q, sigma, t, gamma, delta, prob(:, j));
        iv(:, j) = calIV(call(:, j), S0, K', r, q, t);
    end
    
    figure('PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 10, 8], 'PaperSize', [10, 8]);
    axes('Box', 'on', 'NextPlot', 'add', 'XGrid', 'on', 'YGrid', 'on', 'YLim', [0.2, 0.5]);
    plotH = plot(K' / S0, iv);
    set(plotH(1), 'DisplayName', ['Poisson $\lambda=', num2str(lambda), '$'], 'Color', [0.6, 0.6, 0.6]);
    set(plotH(2), 'DisplayName', ['Poisson $\lambda=', num2str(lambda0(i)), '$'], 'Color', [0, 0, 0], 'LineStyle', '--');
    for ia = 1 : 3
        str = ['$(\alpha,\beta)=(', num2str(a(ia)), ',', num2str(b(ia)), ')$'];
        set(plotH(ia + 2), 'DisplayName', str, 'Color', color(ia, :));
    end
    neworder = [2, 1, 3, 4, 5];
    labels = get(legend(), 'String');
    leg = legend(plotH(neworder), labels(neworder), 'Location', 'NorthEast');
    set(leg, 'Interpreter', 'latex', 'FontSize', 8);
    xlabel('$K/S_0$', 'Interpreter', 'latex');
    tstr = ['(', char(96 + i), ') Conditional ($\lambda_0=', num2str(lambda0(i)), '$, $\lambda=', num2str(lambda), ')$'];
    title(tstr, 'Interpreter', 'latex');
    fstr = ['Fig_IV_C', num2str(lambda0(i))];
    saveas(gcf, fstr, 'pdf');
    
end

prob = zeros(N + 1, 4);
prob(:, 1) = poisPMF(lambda, t, N);
for ia = 1 : 3
    str = ['simData_U_', num2str(ind(ia)), '.mat'];
    load(str);
    prob(:, ia + 1) = P(:, 3);
end
call = zeros(nK, 4);
iv = zeros(nK, 4);
for j = 1 : 4
    call(:, j) = callMJD(S0, K, r, q, sigma, t, gamma, delta, prob(:, j));
    iv(:, j) = calIV(call(:, j), S0, K', r, q, t);
end

figure('PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 10, 8], 'PaperSize', [10, 8]);
axes('Box', 'on', 'NextPlot', 'add', 'XGrid', 'on', 'YGrid', 'on', 'YLim', [0.2, 0.5]);
plotH = plot(K' / S0, iv);
set(plotH(1), 'DisplayName', ['Poisson $\lambda=', num2str(lambda), '$'], 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1);
for ia = 1 : 3
    str = ['$(\alpha,\beta)=(', num2str(a(ia)), ',', num2str(b(ia)), ')$'];
    set(plotH(ia + 1), 'DisplayName', str, 'Color', color(ia, :));
end
leg = legend(gca, 'show', 'Location', 'NorthEast');
set(leg, 'Interpreter', 'latex', 'FontSize', 8);
xlabel('$K/S_0$', 'Interpreter', 'latex');
tstr = ['(d) Unconditional $(\lambda=', num2str(lambda), ')$'];
title(tstr, 'Interpreter', 'latex');
fstr = 'Fig_IV_U';
saveas(gcf, fstr, 'pdf');

