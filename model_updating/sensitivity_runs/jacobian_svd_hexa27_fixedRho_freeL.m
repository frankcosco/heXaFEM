%% jacobian_svd_hexa27_fixedRho_freeL.m
% Follow-up identifiability test (Francesco's question): if density rho is
% FIXED (not fitted) and length L is freed instead, does freeing a geometry
% parameter reproduce the same kind of flat/degenerate direction we found
% for rho-E, or a different one? And does fixing rho (breaking the rho-E
% degeneracy) make E, nu well identified on their own, with L coupling in
% separately (or not at all)?
%
% Mesh: Hexa27-Int27, 13x117x1 (chosen over Hexa64 for this test purely for
% speed: ~5-6s/iteration vs ~53s/iteration for Hexa64 on the same density,
% per the historical 21-04-2026 logs -- same nominal x0, same beam dims).
%
% IMPORTANT ASSUMPTION (flag before trusting the numbers): rho is fixed here
% at the same nominal placeholder value used everywhere else as x0
% (2.70e-06 kg/mm^3), NOT at an independently measured density -- that
% question (do we have a direct mass/volume measurement of the specimen?)
% is still open. For THIS test (probing the identifiability STRUCTURE of
% E, nu, L once rho's degeneracy is removed) the exact rho value shouldn't
% matter much -- but if/when a measured rho becomes available, rerun with
% that value plugged in here instead.
%
% Length bounds (lower_b/upper_b, last entry) are set to a generous +-10%
% numerical probing range (324-396mm) around the nominal 360mm, purely to
% give the optimizer room to move -- NOT a claim about the beam's true
% physical length uncertainty (which should be small/well measured). Adjust
% if a different range is wanted.
%
% Requires aligned_experimental_data.mat in the current folder (shared
% cache, mesh-independent) -- run compute_and_cache_experimental_data.m
% once first if not already present.

clear
addpath(genpath('..\..\'))

if ~isfile('aligned_experimental_data.mat')
    error('Run compute_and_cache_experimental_data.m first (in this same folder).');
end
load('aligned_experimental_data.mat', 'ExpMV', 'ExpLambda', 'ExpMNodes');

%% Fixed problem setup (mesh = Hexa27-Int27, 13x117x1)
a = 40; b = 360; c = 4;
EXA_type = 27; EXA_order = 27;
na = 13; nb = 117; nc = 1;

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

%% Scaling reference / bounds -- rho FIXED, E/nu/length FREE (thickness/width fixed)
x_scale  = [2.70e-06, 69580000, 0.33, 4.0, 40.0, 360];
lower_b  = [2.60e-06, 64e6,     0.30, 3.5,   36,  324];   % length lower: 360*0.9
upper_b  = [2.85e-06, 74e6,     0.38, 4.5,   42,  396];   % length upper: 360*1.1
lower_b_sc = lower_b ./ x_scale;
upper_b_sc = upper_b ./ x_scale;

var2optimize = [false true true false false true];   % rho FIXED; E, nu, length FREE
res2use = [true true true true; true true true true];

%% Starting point: nominal rho/E/nu + nominal geometry
x0 = [2.70e-06, 69580000, 0.33, 4.0, 40.0, 360];
x0_scaled = x0 ./ x_scale;

options = optimoptions('lsqnonlin', 'Display', 'iter-detailed', ...
    'StepTolerance',1e-4, 'Algorithm','trust-region-reflective', ...
    'UseParallel',true, 'FunctionTolerance',1e-5);

fprintf('Running calibration (Hexa27-Int27, 13x117x1, rho fixed, E/nu/length free) to capture the Jacobian at p*...\n');
[x_sc_opt_red, resnorm, residual, exitflag, output, lambda, jacobian] = lsqnonlin( ...
    @(x) objectiveFunction(x, EXA_type, EXA_order, Constraints, Elements, Nodes, ...
         ExpMV, ExpLambda, indexOutputDofs, x_scale, var2optimize, res2use, ExpMNodes), ...
    x0_scaled(var2optimize), lower_b_sc(var2optimize), upper_b_sc(var2optimize), options);

x_sc_opt = ones(1,6);
x_sc_opt(var2optimize) = x_sc_opt_red;
x_opt = x_scale .* x_sc_opt;
fprintf('\np* : rho(fixed)=%.4e  E=%.4e kPa  nu=%.4f  length=%.4f mm (resnorm=%.4f, exitflag=%d)\n', ...
    x_opt(1), x_opt(2), x_opt(3), x_opt(6), resnorm, exitflag);
fprintf('(Reference, historical materialOnly 21-04-2026 log, rho/E/nu free, length FIXED at 360: rho=2.6966e-06 E=6.6738e+07 nu=0.3797)\n');

%% SVD of the Jacobian (scaled variable space: columns = [E, nu, length] -- rho not a free column here)
J = full(jacobian);
fprintf('\nJacobian size: %d residuals x %d parameters\n', size(J,1), size(J,2));

[U,S,V] = svd(J, 'econ');
sv = diag(S);
condNumber = sv(1) / sv(end);

paramNames = {'E','nu','length'};
fprintf('\nSingular values: %s\n', mat2str(sv,6));
fprintf('Condition number (sv_max/sv_min): %.2f\n', condNumber);
fprintf('\nRight singular vectors (columns), scaled-parameter space [E,nu,length]:\n');
disp(array2table(V, 'RowNames', paramNames, ...
    'VariableNames', {'V1_stiffest','V2','V3_flattest'}));

% Correlation matrix (E, nu, length) from inv(J'J) -- same diagnostic used for rho-E-nu
JTJ = J.' * J;
covP = inv(JTJ);
dP = sqrt(diag(covP));
corrP = covP ./ (dP * dP.');
fprintf('\nParameter correlation matrix [E, nu, length]:\n');
disp(array2table(corrP, 'RowNames', paramNames, 'VariableNames', paramNames));

fprintf(['\nInterpretation: with rho fixed, this tests (a) whether E and nu become\n' ...
    'individually well identified (compare condition number/correlations here\n' ...
    'against the rho-E-nu case in jacobian_svd_hexa8_results.mat /\n' ...
    'jacobian_svd_hexa64_small_results.mat), and (b) whether length introduces\n' ...
    'a NEW flat/degenerate direction (e.g. length trading off against E, since\n' ...
    'bending frequency ~ sqrt(E*I)/(rho*A*L^4)-type scaling means length and E\n' ...
    'do NOT enter the eigenproblem with the same exact exponent as rho and E do --\n' ...
    'so this is not expected to be an EXACT degeneracy like rho-E, but could still\n' ...
    'be a strong NUMERICAL one if the modal data do not constrain L tightly enough\n' ...
    'to disambiguate it from E).\n']);

logsDir = fullfile('..','logs');
if ~exist(logsDir,'dir'); mkdir(logsDir); end
save(fullfile(logsDir,'jacobian_svd_hexa27_fixedRho_freeL_results.mat'), 'J','U','S','V','sv','condNumber','corrP','x_opt','resnorm');
fprintf('\nSaved jacobian_svd_hexa27_fixedRho_freeL_results.mat in %s\n', logsDir);
