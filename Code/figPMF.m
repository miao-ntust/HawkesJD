clear

color = get(groot, 'DefaultAxesColorOrder');

lambda = 5;

a = [0.3, 0.4, 0.3];
b = [5, 5, 10];
ab  = a .* b;
r = lambda * (1 - 1 ./ ab);

t = 1;

lambda0 = [2, 5, 8];

N = 20;
n = (0 : N)';

for i = 1 : 3
    P = zeros(N + 1, 5);
    P(:, 1) = poisPMF(lambda, t, N);
    P(:, 2) = poisPMF(lambda0(i), t, N);
    for ia = 1 : 3
        [~, P(:, ia + 2)] = hawkesPMF(r(ia), a(ia), b(ia), t, lambda0(i), N);
    end
    
    figure('PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 10, 8], 'PaperSize', [10, 8]);
    axes('Box', 'on', 'NextPlot', 'add', 'XGrid', 'on', 'YGrid', 'on');
    plotH = plot(n, P);
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
    xlabel('$n$', 'Interpreter', 'latex');
    tstr = ['(', char(96 + i), ') Conditional ($\lambda_0=', num2str(lambda0(i)), '$, $\lambda=', num2str(lambda), ')$'];
    title(tstr, 'Interpreter', 'latex');
    fstr = ['Fig_PMF_C', num2str(lambda0(i))];
    saveas(gcf, fstr, 'pdf');
    
end

P = zeros(N + 1, 4);
P(:, 1) = poisPMF(lambda, t, N);
for ia = 1 : 3
    [P(:, ia + 1), ~] = hawkesPMF(r(ia), a(ia), b(ia), t, 0, N);
end

figure('PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 10, 8], 'PaperSize', [10, 8]);
axes('Box', 'on', 'NextPlot', 'add', 'XGrid', 'on', 'YGrid', 'on');
plotH = plot(n, P);
set(plotH(1), 'DisplayName', ['Poisson $\lambda=', num2str(lambda), '$'], 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1);
for ia = 1 : 3
    str = ['$(\alpha,\beta)=(', num2str(a(ia)), ',', num2str(b(ia)), ')$'];
    set(plotH(ia + 1), 'DisplayName', str, 'Color', color(ia, :));
end
leg = legend(gca, 'show', 'Location', 'NorthEast');
set(leg, 'Interpreter', 'latex', 'FontSize', 8);
xlabel('$n$', 'Interpreter', 'latex');
tstr = ['(d) Unconditional $(\lambda=', num2str(lambda), ')$'];
title(tstr, 'Interpreter', 'latex');
fstr = 'Fig_PMF_U';
saveas(gcf, fstr, 'pdf');
