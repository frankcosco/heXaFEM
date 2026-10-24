%% jacobian_svd_hexa27_fixedRho_fixedL_retry.m
% Standalone retry of the Hexa27-Int27 (13x117x1) case from
% jacobian_svd_three_meshes_fixedRho_fixedL.m: that run hit nu=0.4000
% EXACTLY (the upper bound), with resnorm=0.1852 -- worse than the
% historical free-fit resnorm for this same mesh (0.1675) -- so the true
% unconstrained nu is above 0.40. Isolated here so Hexa8/Hexa64 (which
% converged to interior points) don't need to be rerun.
%
% Same rho reference (Hexa64-Int64 26x234x2's own historical optimum,
% 2.7747e-06) and same fixed length (360mm) as the three-mesh script;
% nu upper bound raised to 0.45.

clear
addpath(genpath('..\..\'))

if ~isfile('aligned_experimental_data.mat')
    error('Run compute_and_cache_experimental_data.m first (in this same folder).');
end
load('aligned_experimental_data.mat', 'ExpMV', 'ExpLambda', 'ExpMNodes');

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

x_scale  = [2.7747e-06, 69580000, 0.33, 4.0, 40.0, 360];   % rho ref = Hexa64-Int64 26x234x2 historical optimum
lower_b  = [2.60e-06,   60e6,     0.26, 3.5,   36,  324];
upper_b  = [2.85e-06,   78e6,     0.45, 4.5,   42,  396];   % nu upper RAISED 0.40 -> 0.45 (0.40 clipped last time)
lower_b_sc = lower_b ./ x_scale;
upper_b_sc = upper_b ./ x_scale;

var2optimize = [false true true false false false];   % rho FIXED, length FIXED; only E, nu free
res2use = [true true true true; true true true true];

x0 = [2.7747e-06, 69580000, 0.33, 4.0, 40.0, 360];
x0_scaled = x0 ./ x_scale;

options = optimoptions('lsqnonlin', 'Display', 'iter-detailed', ...
    'StepTolerance',1e-4, 'Algorithm','trust-region-reflective', ...
    'UseParallel',true, 'FunctionTolerance',1e-5);

fprintf('Running Hexa27-Int27 (13x117x1), rho & length fixed, E/nu free, nu upper=0.45...\n');
[x_sc_opt_red, resnorm, residual, exitflag, output, lambda, jacobian] = lsqnonlin( ...
    @(x) objectiveFunction(x, EXA_type, EXA_order, Constraints, Elements, Nodes, ...
         ExpMV, ExpLambda, indexOutputDofs, x_scale, var2optimize, res2use, ExpMNodes), ...
    x0_scaled(var2optimize), lower_b_sc(var2optimize), upper_b_sc(var2optimize), options);

x_sc_opt = ones(1,6);
x_sc_opt(var2optimize) = x_sc_opt_red;
x_opt = x_scale .* x_sc_opt;

J = full(jacobian);
[U,S,V] = svd(J, 'econ');
sv = diag(S);
condNumber = sv(1) / sv(end);
JTJ = J.' * J;
covP = inv(JTJ);
dP = sqrt(diag(covP));
corrP = covP ./ (dP * dP.');

fprintf('\np*  : E=%.4e kPa  nu=%.4f  (rho fixed=%.4e, length fixed=%.1f, resnorm=%.4f, exitflag=%d)\n', ...
    x_opt(2), x_opt(3), x_opt(1), x_opt(6), resnorm, exitflag);
fprintf('(Reference: historical free-fit resnorm for this mesh = 0.1675; nu bound now 0.45)\n');
fprintf('singular values: %s   condition number: %.3f\n', mat2str(sv,6), condNumber);
fprintf('corr(E,nu) = %.4f\n', corrP(1,2));

logsDir = fullfile('..','logs');
if ~exist(logsDir,'dir'); mkdir(logsDir); end
save(fullfile(logsDir,'jacobian_svd_hexa27_fixedRho_fixedL_retry_results.mat'), 'J','U','S','V','sv','condNumber','corrP','x_opt','resnorm');
fprintf('\nSaved jacobian_svd_hexa27_fixedRho_fixedL_retry_results.mat in %s\n', logsDir);
