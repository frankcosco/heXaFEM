%% verify_rho_E_gauge_invariance.m
% Direct empirical check of the exact rho-E gauge invariance (R3.3/R4.2/R5.3
% evidence): if the modal data truly only constrain the COMBINATION E/rho
% (and nu, and geometry), then for a given optimum (rho*, E*, nu*) we should
% be able to rescale rho* back to ANY reference value rho_ref -- in
% particular the nominal rho_ref = 2.70e-06 kg/mm^3 used as x0 everywhere --
% by setting E_new = E* * (rho_ref/rho*), nu_new = nu* (unchanged), and get
% EXACTLY the same residual/resnorm as at the original optimum. No
% re-optimization needed: this is a single extra function evaluation per
% case, directly testing the theoretical invariance (not an approximation)
% against the actual FE model + experimental data pipeline.
%
% Also reports E_new for each mesh/order case: this is the "E normalized to
% a common rho reference" quantity that can be compared directly ACROSS
% mesh/order cases without the rho-E gauge ambiguity contaminating the
% comparison (useful for presenting compensation results more cleanly).
%
% Requires aligned_experimental_data.mat (shared cache) and the existing
% jacobian_svd_hexa8_results.mat / jacobian_svd_hexa64_small_results.mat
% (today's optima) in this folder / ../logs.

clear
addpath(genpath('..\..\'))

if ~isfile('aligned_experimental_data.mat')
    error('Run compute_and_cache_experimental_data.m first (in this same folder).');
end
load('aligned_experimental_data.mat', 'ExpMV', 'ExpLambda', 'ExpMNodes');

rho_ref = 2.70e-06;   % nominal reference rho -- NOT a measured value, just the common
                      % gauge-fixing point everyone already uses as x0

cases = struct('label',{},'EXA_type',{},'EXA_order',{},'na',{},'nb',{},'nc',{},'resultsFile',{});
cases(1) = struct('label','Hexa8-Int8_26x234x2',   'EXA_type',8,  'EXA_order',8,  'na',26,'nb',234,'nc',2, 'resultsFile','../logs/jacobian_svd_hexa8_results.mat');
cases(2) = struct('label','Hexa64-Int64_13x117x1', 'EXA_type',64, 'EXA_order',64, 'na',13,'nb',117,'nc',1, 'resultsFile','../logs/jacobian_svd_hexa64_small_results.mat');

a = 40; b = 360; c = 4;
results = table();

for ic = 1:numel(cases)
    cs = cases(ic);
    fprintf('\n=== %s ===\n', cs.label);

    [ConstrFE,ElemFE,NodesFE] = brick_hexa8_mesher_advanced(1, a,b,-c, 13,117,1, 0,-b,0);
    logIndConstFE = ~ConstrFE(:,1:3)';
    [~,indexMorphingNodes] = OutputRequests(NodesFE,logIndConstFE,ExpMNodes);
    NodesFE(indexMorphingNodes,2:4) = ExpMNodes(:,2:4);
    NodesRef = NodesFE(indexMorphingNodes,:);

    [Constraints, Elements, Nodes] = brick_hexa_mesher(cs.EXA_type, 1, a, b, -c, cs.na, cs.nb, cs.nc);
    Nodes(:,2:4) = Nodes(:,2:4) + repmat([-0 -b -0*c], size(Nodes,1), 1);
    logicalIndConstrained = ~Constraints(:,1:3)';
    [~,indexMorphingNodes2] = OutputRequests(Nodes, logicalIndConstrained, NodesRef);
    indexOutputDofs = reshape([3*(indexMorphingNodes2-1)+1; 3*(indexMorphingNodes2-1)+2; 3*(indexMorphingNodes2-1)+3],[],1);

    S = load(cs.resultsFile, 'x_opt', 'resnorm');
    x_opt = S.x_opt;            % [rho*, E*, nu*, thickness, width, length]
    rho_star = x_opt(1); E_star = x_opt(2); nu_star = x_opt(3);

    x_scale = [2.70e-06, 69580000, 0.33, 4.0, 40.0, 360];
    var2optimize_all = [true true true false false false];
    res2use = [true true true true; true true true true];

    % Evaluate at the ORIGINAL optimum
    x_orig_scaled = x_opt ./ x_scale;
    Err_orig = objectiveFunction(x_orig_scaled(var2optimize_all), cs.EXA_type, cs.EXA_order, ...
        Constraints, Elements, Nodes, ExpMV, ExpLambda, indexOutputDofs, x_scale, var2optimize_all, res2use, ExpMNodes);
    resnorm_orig = sum(Err_orig.^2);

    % Rescale: rho -> rho_ref, E -> E*(rho_ref/rho*), nu unchanged
    alpha = rho_ref / rho_star;
    E_rescaled = E_star * alpha;
    x_rescaled = [rho_ref, E_rescaled, nu_star, x_opt(4), x_opt(5), x_opt(6)];
    x_resc_scaled = x_rescaled ./ x_scale;
    Err_resc = objectiveFunction(x_resc_scaled(var2optimize_all), cs.EXA_type, cs.EXA_order, ...
        Constraints, Elements, Nodes, ExpMV, ExpLambda, indexOutputDofs, x_scale, var2optimize_all, res2use, ExpMNodes);
    resnorm_resc = sum(Err_resc.^2);

    fprintf('  original : rho=%.4e  E=%.4e kPa  nu=%.4f  resnorm=%.6f  (saved resnorm=%.6f)\n', ...
        rho_star, E_star, nu_star, resnorm_orig, S.resnorm);
    fprintf('  rescaled : rho=%.4e  E=%.4e kPa  nu=%.4f  resnorm=%.6f  (alpha=%.4f)\n', ...
        rho_ref, E_rescaled, nu_star, resnorm_resc, alpha);
    fprintf('  |resnorm_orig - resnorm_resc| = %.3e  (should be ~0 if exact gauge invariance holds)\n', ...
        abs(resnorm_orig - resnorm_resc));

    results = [results; table({cs.label}, rho_star, E_star, nu_star, resnorm_orig, ...
        rho_ref, E_rescaled, resnorm_resc, abs(resnorm_orig-resnorm_resc), ...
        'VariableNames', {'label','rho_star','E_star','nu_star','resnorm_orig', ...
                           'rho_ref','E_at_rho_ref','resnorm_rescaled','abs_diff'})]; %#ok<AGROW>
end

disp(results);
logsDir = fullfile('..','logs');
if ~exist(logsDir,'dir'); mkdir(logsDir); end
writetable(results, fullfile(logsDir,'gauge_invariance_check_results.csv'));
fprintf('\nSaved gauge_invariance_check_results.csv in %s\n', logsDir);
fprintf(['\nIf abs_diff is numerically negligible (~1e-6 or smaller, i.e. finite-precision\n' ...
    'noise) for every case, this directly confirms -- on the real FE+experimental\n' ...
    'pipeline, not just in theory -- that rho and E can be re-gauged to any common\n' ...
    'reference without changing the fit quality, and that E_at_rho_ref is the\n' ...
    'correct quantity to compare across mesh/order cases instead of the raw,\n' ...
    'gauge-dependent E_star.\n']);
