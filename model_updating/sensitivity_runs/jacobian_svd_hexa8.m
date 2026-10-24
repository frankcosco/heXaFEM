%% jacobian_svd_hexa8.m
% R3.3 / R4.2 / R5.3 evidence: SVD of the residual Jacobian at the Hexa8-Int8
% (26x234x2) optimum p*, to characterize how well-conditioned the inverse
% problem is in (rho, E, nu) and WHICH combination of parameters is poorly
% constrained by the modal data (the flat/compensating direction).
%
% Requires aligned_experimental_data.mat in the current folder --
% run compute_and_cache_experimental_data.m once first (shared with
% multi_start_hexa8.m).
%
% Re-runs the SAME calibration as the 'script_nominal' trial in
% multi_start_hexa8.m, but requests lsqnonlin's 7th output (jacobian) at the
% converged solution -- the Jacobian MATLAB's own trust-region-reflective
% solver already computed, in the SCALED variable space (each column already
% normalized by x_scale), which is exactly what SVD needs for a fair
% comparison across rho/E/nu (otherwise raw units would dominate).

clear
addpath(genpath('..\..\'))

if ~isfile('aligned_experimental_data.mat')
    error('Run compute_and_cache_experimental_data.m first (in this same folder).');
end
load('aligned_experimental_data.mat', 'ExpMV', 'ExpLambda', 'ExpMNodes');

%% Fixed problem setup (mesh = Hexa8-Int8, 26x234x2 -- identical to multi_start_hexa8.m)
a = 40; b = 360; c = 4;
EXA_type = 8; EXA_order = 8;
na = 26; nb = 234; nc = 2;

[ConstrFE,ElemFE,NodesFE] = brick_hexa8_mesher_advanced(1, a,b,-c, 13,117,1, 0,-b,0);
logIndConstFE = ~ConstrFE(:,1:3)';
[~,indexMorphingNodes] = OutputRequests(NodesFE,logIndConstFE,ExpMNodes);
NodesFE(indexMorphingNodes,2:4) = ExpMNodes(:,2:4);
NodesRef = NodesFE(indexMorphingNodes,:);

[Constraints, Elements, Nodes] = brick_hexa_mesher(EXA_type, 1, a, b, -c, na, nb, nc);
Nodes(:,2:4) = Nodes(:,2:4) + repmat([-0 -b -0*c], size(Nodes,1), 1);
logicalIndConstrained = ~Constraints(:,1:3)';
[~,indexMorphingNodes2] = OutputRequests(Nodes, logicalIndConstrained, NodesRef);
indexOutputDofs = reshape([3*(indexMorphingNodes2-1)+1; 3*(indexMorphingNodes2-1)+2; 3*(indexMorphingNodes2-1)+3],[],1);

x_scale  = [2.70e-06, 69580000, 0.33, 4.0, 40.0, 360];
lower_b  = [2.60e-06, 64e6,     0.30, 3.5,   36,  360];
upper_b  = [2.85e-06, 74e6,     0.38, 4.5,   42,  360];
lower_b_sc = lower_b ./ x_scale;
upper_b_sc = upper_b ./ x_scale;
var2optimize = [true true true false false false];
res2use = [true true true true; true true true true];

%% Re-run the calibration from the script_nominal starting point, capturing the Jacobian
x0 = [2.70e-06, 69580000, 0.33, 4.0, 40.0, 360];
x0_scaled = x0 ./ x_scale;

options = optimoptions('lsqnonlin', 'Display', 'iter-detailed', ...
    'StepTolerance',1e-4, 'Algorithm','trust-region-reflective', ...
    'UseParallel',true, 'FunctionTolerance',1e-5);

fprintf('Running calibration (same as script_nominal trial) to capture the Jacobian at p*...\n');
[x_sc_opt_red, resnorm, residual, exitflag, output, lambda, jacobian] = lsqnonlin( ...
    @(x) objectiveFunction(x, EXA_type, EXA_order, Constraints, Elements, Nodes, ...
         ExpMV, ExpLambda, indexOutputDofs, x_scale, var2optimize, res2use, ExpMNodes), ...
    x0_scaled(var2optimize), lower_b_sc(var2optimize), upper_b_sc(var2optimize), options);

x_sc_opt = ones(1,6);
x_sc_opt(var2optimize) = x_sc_opt_red;
x_opt = x_scale .* x_sc_opt;
fprintf('\np* : rho=%.4e  E=%.4e kPa  nu=%.4f  (resnorm=%.4f, exitflag=%d)\n', ...
    x_opt(1), x_opt(2), x_opt(3), resnorm, exitflag);

%% SVD of the Jacobian (scaled variable space: columns = [rho, E, nu])
J = full(jacobian);
fprintf('\nJacobian size: %d residuals x %d parameters\n', size(J,1), size(J,2));

[U,S,V] = svd(J, 'econ');
sv = diag(S);
condNumber = sv(1) / sv(end);

paramNames = {'rho','E','nu'};
fprintf('\nSingular values: %s\n', mat2str(sv,6));
fprintf('Condition number (sv_max/sv_min): %.2f\n', condNumber);
fprintf('\nRight singular vectors (columns), scaled-parameter space [rho,E,nu]:\n');
disp(array2table(V, 'RowNames', paramNames, ...
    'VariableNames', {'V1_stiffest','V2','V3_flattest'}));

% Parameter correlation matrix from inv(J'J) -- added for consistency with the
% other sensitivity_runs scripts and to make the corr(rho,E) figure discussed
% in the paper reproducible directly from this script's own output, instead
% of only from an ad-hoc follow-up calculation.
JTJ = J.' * J;
covP = inv(JTJ);
dP = sqrt(diag(covP));
corrP = covP ./ (dP * dP.');
fprintf('\nParameter correlation matrix [rho, E, nu]:\n');
disp(array2table(corrP, 'RowNames', paramNames, 'VariableNames', paramNames));

fprintf(['\nInterpretation: V3_flattest (associated with the SMALLEST singular value,\n' ...
    'sv=%.4g) is the parameter combination the modal data constrains LEAST --\n' ...
    'i.e. the direction along which rho/E/nu can trade off against each other\n' ...
    'while barely changing the residual. Compare its dominant components\n' ...
    '(rho vs E) against the expected omega^2 ~ E/rho degeneracy.\n'], sv(end));

logsDir = fullfile('..','logs');
if ~exist(logsDir,'dir'); mkdir(logsDir); end
save(fullfile(logsDir,'jacobian_svd_hexa8_results.mat'), 'J','U','S','V','sv','condNumber','corrP','x_opt','resnorm');
fprintf('\nSaved jacobian_svd_hexa8_results.mat in %s\n', logsDir);
