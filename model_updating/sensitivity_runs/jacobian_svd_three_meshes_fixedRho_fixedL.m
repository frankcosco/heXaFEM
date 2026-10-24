%% jacobian_svd_three_meshes_fixedRho_fixedL.m
% Closes the loop on the identifiability tests: with rho AND all geometry
% (thickness, width, length) FIXED at nominal, only E and nu are free.
% Repeats this 2-parameter SVD/correlation check on the SAME three
% mesh/order configurations already used in the other tests, so all four
% configurations (rho-E-nu free / rho fixed+E-nu-length free / this one)
% are directly comparable:
%   - Hexa8-Int8,  26x234x2   (same mesh as jacobian_svd_hexa8.m)
%   - Hexa27-Int27,13x117x1   (same mesh as jacobian_svd_hexa27_fixedRho_freeL.m)
%   - Hexa64-Int64,13x117x1   (same mesh as jacobian_svd_hexa64_small.m)
%
% rho fixed at the nominal placeholder (2.70e-06 kg/mm^3, same as every
% other script here); length/thickness/width fixed at the nominal beam
% dimensions (4,40,360). Expectation: with BOTH known degenerate directions
% (rho-E and E-length) removed by fixing rho and length, the remaining
% 2-parameter (E,nu) problem should be noticeably better conditioned than
% either 3-parameter case -- this run quantifies exactly how much.
%
% Requires aligned_experimental_data.mat in the current folder (shared
% cache) -- run compute_and_cache_experimental_data.m once first if not
% already present.

clear
addpath(genpath('..\..\'))

if ~isfile('aligned_experimental_data.mat')
    error('Run compute_and_cache_experimental_data.m first (in this same folder).');
end
load('aligned_experimental_data.mat', 'ExpMV', 'ExpLambda', 'ExpMNodes');

a = 40; b = 360; c = 4;

cases = struct('label',{},'EXA_type',{},'EXA_order',{},'na',{},'nb',{},'nc',{});
cases(1) = struct('label','Hexa8-Int8_26x234x2',    'EXA_type',8,  'EXA_order',8,  'na',26,'nb',234,'nc',2);
cases(2) = struct('label','Hexa27-Int27_13x117x1',  'EXA_type',27, 'EXA_order',27, 'na',13,'nb',117,'nc',1);
cases(3) = struct('label','Hexa64-Int64_13x117x1',  'EXA_type',64, 'EXA_order',64, 'na',13,'nb',117,'nc',1);

x_scale  = [2.7747e-06, 69580000, 0.33, 4.0, 40.0, 360];   % rho ref = Hexa64-Int64 26x234x2's own
                                                            % historical optimum (rho*=2.7747e-06),
                                                            % matching the manuscript's asymptotic
                                                            % reference p*inf (page 13) -- NOT an
                                                            % arbitrary nominal placeholder.
lower_b  = [2.60e-06, 60e6,     0.26, 3.5,   36,  324];   % small safety margin (not the original 64e6/0.30)
upper_b  = [2.85e-06, 78e6,     0.40, 4.5,   42,  396];   % small safety margin (not the original 74e6/0.38)
lower_b_sc = lower_b ./ x_scale;
upper_b_sc = upper_b ./ x_scale;

var2optimize = [false true true false false false];   % rho FIXED, length FIXED; only E, nu free
res2use = [true true true true; true true true true];

x0 = [2.7747e-06, 69580000, 0.33, 4.0, 40.0, 360];   % rho component unused for optimization (fixed), kept consistent with x_scale
x0_scaled = x0 ./ x_scale;

options = optimoptions('lsqnonlin', 'Display', 'iter-detailed', ...
    'StepTolerance',1e-4, 'Algorithm','trust-region-reflective', ...
    'UseParallel',true, 'FunctionTolerance',1e-5);

summaryRows = table();
allResults = struct();

for ic = 1:numel(cases)
    cs = cases(ic);
    fprintf('\n================ %s (rho & length fixed; E,nu free) ================\n', cs.label);

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

    started = tic;
    [x_sc_opt_red, resnorm, residual, exitflag, output, lambda, jacobian] = lsqnonlin( ...
        @(x) objectiveFunction(x, cs.EXA_type, cs.EXA_order, Constraints, Elements, Nodes, ...
             ExpMV, ExpLambda, indexOutputDofs, x_scale, var2optimize, res2use, ExpMNodes), ...
        x0_scaled(var2optimize), lower_b_sc(var2optimize), upper_b_sc(var2optimize), options);
    elapsed = toc(started);

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

    fprintf('  p*  : E=%.4e kPa  nu=%.4f  (rho fixed=%.4e, length fixed=%.1f, resnorm=%.4f, exitflag=%d, %.1fs)\n', ...
        x_opt(2), x_opt(3), x_opt(1), x_opt(6), resnorm, exitflag, elapsed);
    fprintf('  singular values: %s   condition number: %.3f\n', mat2str(sv,6), condNumber);
    fprintf('  corr(E,nu) = %.4f\n', corrP(1,2));

    summaryRows = [summaryRows; table({cs.label}, x_opt(2), x_opt(3), resnorm, exitflag, elapsed, ...
        sv(1), sv(end), condNumber, corrP(1,2), ...
        'VariableNames', {'label','E_opt','nu_opt','resnorm','exitflag','elapsed_s', ...
                           'sv_max','sv_min','condNumber','corr_E_nu'})]; %#ok<AGROW>

    allResults.(matlab.lang.makeValidName(cs.label)) = struct('x_opt',x_opt,'resnorm',resnorm, ...
        'J',J,'U',U,'S',S,'V',V,'sv',sv,'condNumber',condNumber,'corrP',corrP);
end

fprintf('\n================ SUMMARY (rho & length fixed; E, nu free) ================\n');
disp(summaryRows);

fprintf(['\nCompare condNumber/corr_E_nu here against:\n' ...
    '  - rho,E,nu free (geometry fixed):      Hexa8 condNumber=172.6, Hexa64 condNumber=422.0\n' ...
    '  - rho fixed, E,nu,length free (Hexa27): condNumber=56.6, corr(E,length)=0.993\n' ...
    'to see how much fixing BOTH rho and length improves (or does not improve)\n' ...
    'the identifiability of E and nu specifically.\n']);

logsDir = fullfile('..','logs');
if ~exist(logsDir,'dir'); mkdir(logsDir); end
writetable(summaryRows, fullfile(logsDir,'jacobian_svd_three_meshes_fixedRho_fixedL_summary.csv'));
save(fullfile(logsDir,'jacobian_svd_three_meshes_fixedRho_fixedL_results.mat'), 'allResults','summaryRows');
fprintf('\nSaved jacobian_svd_three_meshes_fixedRho_fixedL_summary.csv / _results.mat in %s\n', logsDir);
