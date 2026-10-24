%% jacobian_svd_hexa64_small.m
% Cross-mesh validation of jacobian_svd_hexa8.m: repeats the SAME SVD
% identifiability analysis (R3.3 / R4.2 / R5.3 evidence) on a DIFFERENT
% element type/order and a DIFFERENT (coarser) density, to check that the
% E-rho flat direction found for Hexa8-Int8 (26x234x2) is a property of the
% modal-data parametrization itself, not an artifact of one specific mesh.
%
% Mesh: Hexa64-Int64, 13x117x1 (the SMALLER/coarser of the two Hexa64
% densities in the DoE, Table 2 -- matches
% dev2/.../revisedFrequencies/materialOnly/MU-Singlerun-TRRBeamHexa64-Int64_13x117x1.txt,
% same x0 = [2.70e-06, 69580000, 0.33] as the Hexa8 'script_nominal' case).
%
% Requires aligned_experimental_data.mat in the current folder (SHARED with
% jacobian_svd_hexa8.m / multi_start_hexa8.m -- experimental data does not
% depend on the FE mesh, so the same cache applies) -- run
% compute_and_cache_experimental_data.m once first if not already present.

clear
addpath(genpath('..\..\'))

if ~isfile('aligned_experimental_data.mat')
    error('Run compute_and_cache_experimental_data.m first (in this same folder).');
end
load('aligned_experimental_data.mat', 'ExpMV', 'ExpLambda', 'ExpMNodes');

%% Fixed problem setup (mesh = Hexa64-Int64, 13x117x1 -- the smaller Hexa64 density)
a = 40; b = 360; c = 4;
EXA_type = 64; EXA_order = 64;
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

% Same scaling reference / bounds / var2optimize as jacobian_svd_hexa8.m --
% keeps the two analyses directly comparable (only mesh/element changes).
x_scale  = [2.70e-06, 69580000, 0.33, 4.0, 40.0, 360];
lower_b  = [2.60e-06, 64e6,     0.30, 3.5,   36,  360];
upper_b  = [2.85e-06, 74e6,     0.38, 4.5,   42,  360];
lower_b_sc = lower_b ./ x_scale;
upper_b_sc = upper_b ./ x_scale;
var2optimize = [true true true false false false];   % rho, E, nu free (geometry fixed) -- FIXED 24/10: was accidentally left as [false true true false false true] (rho fixed, length free) from reuse as a template for the later fixedRho scripts; the saved jacobian_svd_hexa64_small_results.mat (rho*=2.7254e-06 != x0's 2.70e-06, condNumber=422.0) confirms rho WAS free when that result was produced, so this restores the script to what actually generated it.
res2use = [true true true true; true true true true];

%% Re-run the calibration from the same nominal starting point, capturing the Jacobian
x0 = [2.70e-06, 69580000, 0.33, 4.0, 40.0, 360];
x0_scaled = x0 ./ x_scale;

options = optimoptions('lsqnonlin', 'Display', 'iter-detailed', ...
    'StepTolerance',1e-4, 'Algorithm','trust-region-reflective', ...
    'UseParallel',true, 'FunctionTolerance',1e-5);

fprintf('Running calibration (Hexa64-Int64, 13x117x1) to capture the Jacobian at p*...\n');
[x_sc_opt_red, resnorm, residual, exitflag, output, lambda, jacobian] = lsqnonlin( ...
    @(x) objectiveFunction(x, EXA_type, EXA_order, Constraints, Elements, Nodes, ...
         ExpMV, ExpLambda, indexOutputDofs, x_scale, var2optimize, res2use, ExpMNodes), ...
    x0_scaled(var2optimize), lower_b_sc(var2optimize), upper_b_sc(var2optimize), options);

x_sc_opt = ones(1,6);
x_sc_opt(var2optimize) = x_sc_opt_red;
x_opt = x_scale .* x_sc_opt;
fprintf('\np* : rho=%.4e  E=%.4e kPa  nu=%.4f  (resnorm=%.4f, exitflag=%d)\n', ...
    x_opt(1), x_opt(2), x_opt(3), resnorm, exitflag);
fprintf('(Reference, historical 21-04-2026 log: rho=2.7807e-06 E=6.9348e+07 nu=0.3334)\n');

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

% Parameter correlation matrix from inv(J'J) -- added 24/10, same reasoning as
% jacobian_svd_hexa8.m: makes the corr(rho,E)=1.0000 figure cited in the reviewer
% response reproducible directly from this script instead of an ad-hoc calculation.
JTJ = J.' * J;
covP = inv(JTJ);
dP = sqrt(diag(covP));
corrP = covP ./ (dP * dP.');
fprintf('\nParameter correlation matrix [rho, E, nu]:\n');
disp(array2table(corrP, 'RowNames', paramNames, 'VariableNames', paramNames));

fprintf(['\nInterpretation: compare V3_flattest here against the Hexa8-Int8 result\n' ...
    '(jacobian_svd_hexa8_results.mat). If both meshes/element types show the same\n' ...
    'flat direction (rho,E with near-zero nu component) and a comparable condition\n' ...
    'number, the E-rho non-identifiability is a property of the modal-data\n' ...
    'parametrization itself, not an artifact of one specific mesh/element order.\n']);

logsDir = fullfile('..','logs');
if ~exist(logsDir,'dir'); mkdir(logsDir); end
save(fullfile(logsDir,'jacobian_svd_hexa64_small_results.mat'), 'J','U','S','V','sv','condNumber','corrP','x_opt','resnorm');
fprintf('\nSaved jacobian_svd_hexa64_small_results.mat in %s\n', logsDir);
