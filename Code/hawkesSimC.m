function Pi = hawkesSimC(r, lambda0, a, b, t, paths)
% Hawkes Processes Simulation (Conditional)
%   Input:
%       r - intensity baseline
%       lambda0 - initial lambda
%       a - intensity jump rate
%       b - decay rate
%       t - time (row vector)
%       paths - simulation paths (x10000)
%
%   Output:
%       Pi - (N, t)

nPath = 10000;
rng('shuffle');
nt = numel(t);
% maxNumJump = round(gammaincinv(1e-10, r * a * b / (a * b - 1) * t(end), 'upper'));
maxNumJump = 200;
Pi = zeros(maxNumJump, nt, paths);

for ipath = 1 : paths

    jumpY = zeros(nPath, maxNumJump + 1);
    jumpY(:, 1) = lambda0 - r;
    jumpY(:, 2 : end) = -log(rand(nPath, maxNumJump)) / a;
    jumpT = zeros(nPath, maxNumJump + 1);
    jumpLambda = r + jumpY(:, 1);
    for i = 2 : (maxNumJump + 1)
        jumpCheck = true(nPath, 1);
        tt = jumpT(:, i - 1);
        indNG = find(tt < max(t));
        nindNG = numel(indNG);
        jumpCheck(indNG) = false;
        tt(jumpCheck) = [];
        jumpT(jumpCheck, i) = jumpT(jumpCheck, i - 1);
        while any(~jumpCheck)
            tt = -log(rand(nindNG, 1)) ./ jumpLambda(indNG) + tt;
            lambda = r + sum(jumpY(indNG, 1 : (i - 1)) .* exp(b * jumpT(indNG, 1 : (i - 1))), 2) .* exp(-b * tt);
            ind = find((rand(nindNG, 1) .* jumpLambda(indNG)) < lambda);
            jumpLambda(indNG) = lambda;
            indOK = indNG(ind);
            indNG(ind) = [];
            nindNG = numel(indNG);
            jumpCheck(indOK) = true;
            jumpT(indOK, i) = tt(ind);
            tt(ind) = [];
            jumpLambda(indOK) = lambda(ind) + jumpY(indOK, i);
        end
        if ~any(jumpT(:, i) < 1)
            break
        end
    end
    
    for it = 1 : nt
        jumpFlag = (jumpT(:, 2 : end) < t(it)) & (jumpT(:, 2 : end) ~= 0);
        N = sort(sum(jumpFlag, 2), 1);
        i = [find(N(1 : (end - 1)) ~= N(2 : end)); numel(N)];
        Pi(N(i) + 1, it, ipath) = diff([0; i]);
    end
    
    ipath

end
Pi = sum(Pi, 3);
Pi = Pi ./ sum(Pi, 1);

end
