%% jacobian_svd_along_trajectory_hexa8.m
% Francesco's question: is the rho-E ill-conditioning a property of the
% WHOLE parameter space (present throughout the optimization trajectory,
% as the exact invariance argument predicts), or does it only show up near
% the final converged point (which jacobian_svd_hexa8.m alone cannot
% distinguish, since it only evaluates the Jacobian AT the optimum)?
%
% Captures the optimizer's iterate trajectory (OutputFcn) for a single
% Hexa8-Int8 (26x234x2) run, then independently recomputes a Jacobian via a
% FIXED central-difference scheme (NOT lsqnonlin's own internal Jacobian
% estimate, to keep the comparison consistent across points) at three
% points along that trajectory: the initial guess x0, an intermediate
% iterate (~half the iterations), and the final converged point x*.
%
% Requires aligned_experimental_data.mat (shared cache).

clear
clear captureTrajectory   % reset the persistent trajectory store from any previous run
addpath(genpath('..\..\'))

if ~isfile('aligned_experimental_data.mat')
    error('Run compute_and_cache_experimental_data.m first (in this same folder).');
end
load('aligned_experimental_data.mat', 'ExpMV', 'ExpLambda', 'ExpMNodes');

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

x0 = [2.70e-06, 69580000, 0.33, 4.0, 40.0, 360];
x0_scaled = x0 ./ x_scale;

objfun = @(x) objectiveFunction(x, EXA_type, EXA_order, Constraints, Elements, Nodes, ...
    ExpMV, ExpLambda, indexOutputDofs, x_scale, var2optimize, res2use, ExpMNodes);

%% Capture the iterate trajectory via OutputFcn
options = optimoptions('lsqnonlin', 'Display', 'iter-detailed', ...
    'StepTolerance',1e-4, 'Algorithm','trust-region-reflective', ...
    'UseParallel',true, 'FunctionTolerance',1e-5, ...
    'OutputFcn', @captureTrajectory);

fprintf('Running Hexa8-Int8 (26x234x2) and capturing the iterate trajectory...\n');
[x_sc_opt_red, resnorm, residual, exitflag] = lsqnonlin( ...
    objfun, x0_scaled(var2optimize), lower_b_sc(var2optimize), upper_b_sc(var2optimize), options);

trajectory = captureTrajectory([], [], '__retrieve__');
nIter = numel(trajectory);
fprintf('Captured %d iterates.\n', nIter);

idxInit = 1;
idxMid  = max(1, round(nIter/2));
idxFinal= nIter;

points = struct('label',{},'x_red',{});
points(1) = struct('label','x0 (initial)',               'x_red', trajectory{idxInit});
points(2) = struct('label', sprintf('iterate #%d (mid)',idxMid), 'x_red', trajectory{idxMid});
points(3) = struct('label','x* (final/converged)',       'x_red', trajectory{idxFinal});

%% Fixed central-difference Jacobian (independent of lsqnonlin's own estimate)
h = 1e-6;
paramNames = {'rho','E','nu'};
summaryRows = table();

for ip = 1:numel(points)
    pt = points(ip);
    fprintf('\n=== %s ===\n', pt.label);
    x_full_sc = ones(1,6); x_full_sc(var2optimize) = pt.x_red;
    x_full = x_scale .* x_full_sc;
    fprintf('  point (rho,E,nu) = %.4e, %.4e, %.4f\n', x_full(1), x_full(2), x_full(3));

    J = fdJacobian(objfun, pt.x_red, h);
    [~,S,V] = svd(J, 'econ');
    sv = diag(S);
    condNumber = sv(1) / sv(end);
    fprintf('  singular values: %s\n', mat2str(sv,6));
    fprintf('  condition number: %.3f\n', condNumber);
    fprintf('  V3 (flattest) = [%.4f %.4f %.4f]  (rho,E,nu)\n', V(1,end), V(2,end), V(3,end));

    summaryRows = [summaryRows; table({pt.label}, x_full(1), x_full(2), x_full(3), ...
        sv(1), sv(end), condNumber, V(1,end), V(2,end), V(3,end), ...
        'VariableNames', {'label','rho','E','nu','sv_max','sv_min','condNumber','V3_rho','V3_E','V3_nu'})]; %#ok<AGROW>
end

disp(summaryRows);
logsDir = fullfile('..','logs');
if ~exist(logsDir,'dir'); mkdir(logsDir); end
writetable(summaryRows, fullfile(logsDir,'jacobian_svd_along_trajectory_hexa8_summary.csv'));
fprintf('\nSaved jacobian_svd_along_trajectory_hexa8_summary.csv in %s\n', logsDir);
fprintf(['\nIf condNumber stays consistently high (same order of magnitude) and V3 stays\n' ...
    'close to the same (rho,E)-dominated direction across all three points, this confirms\n' ...
    'the ill-conditioning is a structural property present throughout parameter space,\n' ...
    'not something specific to the converged optimum.\n']);

function out = captureTrajectory(x, ~, state)
    persistent traj
    if isempty(traj); traj = {}; end
    if ischar(state) && strcmp(state,'__retrieve__')
        out = traj;
        return
    end
    if strcmp(state,'iter') || strcmp(state,'init')
        traj{end+1} = x; %#ok<AGROW>
    end
    out = false;   % 'stop' flag expected by lsqnonlin's OutputFcn
end

function J = fdJacobian(objfun, x, h)
    n = numel(x);
    f0 = objfun(x); %#ok<NASGU>
    m = numel(objfun(x));
    J = zeros(m,n);
    for i = 1:n
        xp = x; xp(i) = xp(i) + h;
        xm = x; xm(i) = xm(i) - h;
        fp = objfun(xp);
        fm = objfun(xm);
        J(:,i) = (fp - fm) / (2*h);
    end
end
