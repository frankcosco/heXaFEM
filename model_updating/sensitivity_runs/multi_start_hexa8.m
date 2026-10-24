%% multi_start_hexa8.m
% R3.1b evidence: does Hexa8-Int8 (mesh 26x234x2) re-converge to the same
% "compensated" (rho,E,nu) region regardless of starting point, or does it
% depend on initialization (i.e. is the documented compensated solution an
% artifact of always starting from the same nominal x0)?
%
% Requires aligned_experimental_data.mat in the current folder --
% run compute_and_cache_experimental_data.m once first.
%
% Fixed throughout: mesh = Hexa8-Int8, 26x234x2 (k=2 in the original script's
% DoE table); thickness/width/length FIXED at [4, 40, 360] (material-only,
% matching the paper's Fig. 8 configuration -- see README.md / the header
% comment in execution/Modal_updating_clean_LM.m for why this was corrected);
% only rho, E, nu are free (var2optimize = [true true true false false false]).
%
% Varies: the starting point p0 for rho, E, nu only. Same bounds and same
% scaling reference (x_scale) for every trial, so results are directly comparable.

clear
addpath(genpath('..\..\'))

if ~isfile('aligned_experimental_data.mat')
    error('Run compute_and_cache_experimental_data.m first (in this same folder).');
end
load('aligned_experimental_data.mat', 'ExpMV', 'ExpLambda', 'ExpMNodes');

%% Fixed problem setup (mesh = Hexa8-Int8, 26x234x2 -- identical to k=2 in the DoE table)
a = 40; b = 360; c = 4;      % beam dims, mm (Length/Width/Thickness convention as in the original scripts)
EXA_type = 8; EXA_order = 8;
na = 26; nb = 234; nc = 2;

% "design mesh prototype" used only to snap a reference node set to the DIC grid
% (identical logic to Modal_updating_clean_LM.m lines ~72-80)
[ConstrFE,ElemFE,NodesFE] = brick_hexa8_mesher_advanced(1, a,b,-c, 13,117,1, 0,-b,0);
logIndConstFE = ~ConstrFE(:,1:3)';
[~,indexMorphingNodes] = OutputRequests(NodesFE,logIndConstFE,ExpMNodes);
NodesFE(indexMorphingNodes,2:4) = ExpMNodes(:,2:4);
NodesRef = NodesFE(indexMorphingNodes,:);

% actual analysis mesh for this case
[Constraints, Elements, Nodes] = brick_hexa_mesher(EXA_type, 1, a, b, -c, na, nb, nc);
Nodes(:,2:4) = Nodes(:,2:4) + repmat([-0 -b -0*c], size(Nodes,1), 1);
logicalIndConstrained = ~Constraints(:,1:3)';
[~,indexMorphingNodes2] = OutputRequests(Nodes, logicalIndConstrained, NodesRef);
indexOutputDofs = reshape([3*(indexMorphingNodes2-1)+1; 3*(indexMorphingNodes2-1)+2; 3*(indexMorphingNodes2-1)+3],[],1);

%% Fixed scaling reference, bounds, var2optimize, residual selection
x_scale  = [2.70e-06, 69580000, 0.33, 4.0, 40.0, 360];   % same reference for every trial
lower_b  = [2.60e-06, 64e6,     0.30, 3.5,   36,  360];
upper_b  = [2.85e-06, 74e6,     0.38, 4.5,   42,  360];
lower_b_sc = lower_b ./ x_scale;
upper_b_sc = upper_b ./ x_scale;

var2optimize = [true true true false false false];   % rho, E, nu free; thickness/width/length fixed

res2use = [true true true true; true true true true];

%% Starting points to test (rho [kg/mm3], E [kPa], nu) -- thickness/width/length always [4,40,360]
trials = struct( ...
    'label', {}, 'rho0', {}, 'E0', {}, 'nu0', {});
trials(1) = struct('label','nominal_paper',      'rho0',2.711e-06, 'E0',68980000, 'nu0',0.330);
trials(2) = struct('label','script_nominal',     'rho0',2.70e-06,  'E0',69580000, 'nu0',0.330);
trials(3) = struct('label','hexa64_asymptotic',  'rho0',2.772e-06, 'E0',69405000, 'nu0',0.3314);
trials(4) = struct('label','bounds_low_corner',  'rho0',2.61e-06,  'E0',64.5e6,   'nu0',0.305);
trials(5) = struct('label','bounds_high_corner', 'rho0',2.84e-06,  'E0',73.5e6,   'nu0',0.375);

results = table();

options = optimoptions('lsqnonlin', 'Display', 'iter-detailed', ...
    'StepTolerance',1e-4, 'Algorithm','trust-region-reflective', ...
    'FunctionTolerance',1e-5, 'UseParallel',true);

fprintf('Starting multi-start sensitivity study: %d trials.\n', numel(trials));

for it = 1:numel(trials)
    tr = trials(it);
    x0 = [tr.rho0, tr.E0, tr.nu0, 4.0, 40.0, 360];
    x0_scaled = x0 ./ x_scale;

    fprintf('\n>>> Trial %d/%d: %s -- starting (rho0=%.4e, E0=%.4e, nu0=%.4f)...\n', ...
        it, numel(trials), tr.label, tr.rho0, tr.E0, tr.nu0);

    [Error_ini,Err_pi_reali_ini] = objectiveFunction(x0_scaled(var2optimize), ...
        EXA_type, EXA_order, Constraints, Elements, Nodes, ...
        ExpMV, ExpLambda, indexOutputDofs, x_scale, var2optimize, res2use, ExpMNodes);

    fprintf('    initial E_pi: %s (resnorm=%.4f)\n', mat2str(Err_pi_reali_ini,4), sum(Err_pi_reali_ini));

    started = tic;
    [x_sc_opt_red, resnorm, residual, exitflag] = lsqnonlin( ...
        @(x) objectiveFunction(x, EXA_type, EXA_order, Constraints, Elements, Nodes, ...
             ExpMV, ExpLambda, indexOutputDofs, x_scale, var2optimize, res2use, ExpMNodes), ...
        x0_scaled(var2optimize), lower_b_sc(var2optimize), upper_b_sc(var2optimize), options);
    elapsed = toc(started);

    x_sc_opt = ones(1,6);
    x_sc_opt(var2optimize) = x_sc_opt_red;
    x_opt = x_scale .* x_sc_opt;

    fprintf('\n<<< Trial %d/%d DONE: %s ===\n', it, numel(trials), tr.label);
    fprintf('  p0  : rho=%.4e  E=%.4e kPa  nu=%.4f\n', tr.rho0, tr.E0, tr.nu0);
    fprintf('  p*  : rho=%.4e  E=%.4e kPa  nu=%.4f  (resnorm=%.4f, exitflag=%d, %.1fs)\n', ...
        x_opt(1), x_opt(2), x_opt(3), resnorm, exitflag, elapsed);

    results = [results; table({tr.label}, tr.rho0, tr.E0, tr.nu0, ...
        x_opt(1), x_opt(2), x_opt(3), resnorm, exitflag, elapsed, ...
        'VariableNames', {'label','rho0','E0','nu0','rho_opt','E_opt','nu_opt','resnorm','exitflag','elapsed_s'})]; %#ok<AGROW>
end

disp(results);
logsDir = fullfile('..','logs');
if ~exist(logsDir,'dir'); mkdir(logsDir); end
writetable(results, fullfile(logsDir,'multi_start_hexa8_results.csv'));
save(fullfile(logsDir,'multi_start_hexa8_results.mat'), 'results');
fprintf('\nSaved multi_start_hexa8_results.csv / .mat in %s\n', logsDir);
