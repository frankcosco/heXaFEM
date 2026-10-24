function S = parseModelUpdatingReport(filename)
%PARSEMODELUPDATINGREPORT Parse MATLAB text report for model updating runs
%
% S = parseModelUpdatingReport(filename)
%
% Output struct S contains:
%   S.file
%   S.beamDims
%   S.mesh
%   S.formulation
%   S.dofs.free / fixed / total
%   S.iterTable
%   S.errHistory
%   S.optimized
%   S.summary
%   S.elapsedTime
%
% Example:
%   S = parseModelUpdatingReport('ModelUpdatingLMBeamHexa27-Int27_13x117x1.txt');

    txt = fileread(filename);

    S = struct();
    S.file = filename;

    %% Beam dimensions
    tok = regexp(txt, ...
        'Beam Dimensions:\s*([0-9Ee+.\-]+)\s*x\s*([0-9Ee+.\-]+)\s*x\s*([0-9Ee+.\-]+)', ...
        'tokens', 'once');
    if ~isempty(tok)
        S.beamDims = str2double(tok);
    else
        S.beamDims = [NaN NaN NaN];
    end

    %% Mesh / formulation
    tok = regexp(txt, ...
        'Discretized with:\s*([0-9]+)x([0-9]+)x([0-9]+)\s+([A-Za-z0-9\-]+)\s+elements', ...
        'tokens', 'once');
    if ~isempty(tok)
        S.mesh.nx = str2double(tok{1});
        S.mesh.ny = str2double(tok{2});
        S.mesh.nz = str2double(tok{3});
        S.formulation = strtrim(tok{4});
    else
        S.mesh = struct('nx',NaN,'ny',NaN,'nz',NaN);
        S.formulation = "";
    end

    %% DOFs
    tok = regexp(txt, ...
        'free DOFs:\s*([0-9]+)\s*-\s*fixed DOFs:\s*([0-9]+)\s*-\s*total DOFs:\s*([0-9]+)', ...
        'tokens', 'once');
    if ~isempty(tok)
        S.dofs.free  = str2double(tok{1});
        S.dofs.fixed = str2double(tok{2});
        S.dofs.total = str2double(tok{3});
    else
        S.dofs = struct('free',NaN,'fixed',NaN,'total',NaN);
    end

    %% Initial objective function line
    tok = regexp(txt, ...
        ['First Objective Function Evaluation:\s*initial E_pi:\s*' ...
         '([0-9Ee+.\-]+)\s+([0-9Ee+.\-]+)\s+([0-9Ee+.\-]+)\s+([0-9Ee+.\-]+)\s*' ...
         '>\s*resNorm\(sum E_pi\):\s*([0-9Ee+.\-]+)'], ...
        'tokens', 'once');
    if ~isempty(tok)
        S.summary.initialErrPI_fromHeader = str2double(tok(1:4));
        S.summary.initialResNorm_fromHeader = str2double(tok{5});
    else
        S.summary.initialErrPI_fromHeader = nan(1,4);
        S.summary.initialResNorm_fromHeader = NaN;
    end

    %% Iteration table
    % Rows like:
    %  0   7   0.716789    0.227528   0.309
    patIter = ['^\s*' ...
        '([0-9]+)\s+' ...              % iteration
        '([0-9]+)\s+' ...              % func-count
        '([0-9Ee+.\-]+)\s*' ...        % resnorm
        '([0-9Ee+.\-]*)\s*' ...        % step (may be empty at iter 0)
        '([0-9Ee+.\-]*)\s*$'];         % optimality (may be empty)
    toks = regexp(txt, patIter, 'tokens', 'lineanchors');

    n = numel(toks);
    iter = nan(n,1);
    funcCount = nan(n,1);
    resnorm = nan(n,1);
    step = nan(n,1);
    optimality = nan(n,1);

    for i = 1:n
        ti = toks{i};
        iter(i) = str2double(ti{1});
        funcCount(i) = str2double(ti{2});
        resnorm(i) = str2double(ti{3});

        if ~isempty(ti{4})
            step(i) = str2double(ti{4});
        end
        if ~isempty(ti{5})
            optimality(i) = str2double(ti{5});
        end
    end

    S.iterTable = table(iter, funcCount, resnorm, step, optimality);

    %% Err_PI history
    % Example:
    % Err_PI: [ 0.1050  0.0940  0.3730  0.1449 ]. ResNorm (sum(ErrP)=)0.7167894. RMS_pi=0.2120990
    patErr = ['Err_PI:\s*\[\s*' ...
        '([0-9Ee+.\-]+)\s+' ...
        '([0-9Ee+.\-]+)\s+' ...
        '([0-9Ee+.\-]+)\s+' ...
        '([0-9Ee+.\-]+)\s*\]\.\s*' ...
        'ResNorm\s*\(sum\(ErrP\)=\)\s*([0-9Ee+.\-]+)\.\s*' ...
        'RMS_pi\s*=\s*([0-9Ee+.\-]+)'];

    toks = regexp(txt, patErr, 'tokens');
    n = numel(toks);

    err1 = nan(n,1);
    err2 = nan(n,1);
    err3 = nan(n,1);
    err4 = nan(n,1);
    resNormErr = nan(n,1);
    rmsPi = nan(n,1);

    for i = 1:n
        ti = toks{i};
        err1(i) = str2double(ti{1});
        err2(i) = str2double(ti{2});
        err3(i) = str2double(ti{3});
        err4(i) = str2double(ti{4});
        resNormErr(i) = str2double(ti{5});
        rmsPi(i) = str2double(ti{6});
    end

    S.errHistory = table(err1, err2, err3, err4, resNormErr, rmsPi, ...
        'VariableNames', {'err_f1','err_f2','err_f3','err_f4','resNorm','RMS'});

    %% Elapsed times
    toks = regexp(txt, 'elapsed Time\s*([0-9Ee+.\-]+)\s*\[s\]', 'tokens');
    if ~isempty(toks)
        S.elapsedTime.samples = cellfun(@(c) str2double(c{1}), toks(:));
    else
        S.elapsedTime.samples = [];
    end



    tok = regexp(txt, 'Model updating\s*([0-9Ee+.\-]+)\s*seconds', 'tokens', 'once');
    if ~isempty(tok)
        S.elapsedTime.total = str2double(tok{1});
    else
        S.elapsedTime.total = NaN;
    end

    % switch S.formulation
    %     case {"Hexa8-Int8_78x702x6","Hexa20-Int27_39x351x3","Hexa27-Int27_39x351x3","Hexa64-Int64_26x234x2"}            
    %         S.elapsedTime.samples=S.elapsedTime.samples/3;
    % end
    %% Combined iteration history
    nHist = min([height(S.iterTable), height(S.errHistory), numel(S.elapsedTime.samples)]);

    S.history = S.iterTable(1:nHist,:);
    S.history.err_f1 = S.errHistory.err_f1(1:nHist);
    S.history.err_f2 = S.errHistory.err_f2(1:nHist);
    S.history.err_f3 = S.errHistory.err_f3(1:nHist);
    S.history.err_f4 = S.errHistory.err_f4(1:nHist);
    S.history.resNormErr = S.errHistory.resNorm(1:nHist);
    S.history.RMS = S.errHistory.RMS(1:nHist);
    S.history.time_s = S.elapsedTime.samples(1:nHist);

    S.history.resnorm_norm = S.history.resnorm ./ S.history.resnorm(1);
    S.history.RMS_norm = S.history.RMS ./ S.history.RMS(1);
    S.history.resNormErr_norm = S.history.resNormErr ./ S.history.resNormErr(1);

    %% Convergence-gap normalization
    J0 = S.history.resnorm(1);
    Jstar = min(S.history.resnorm);

    den =  Jstar;

    if den > eps
        S.history.gap_norm = (S.history.resnorm )/den;
    else
        S.history.gap_norm = zeros(height(S.history),1);
    end

    S.history.gap_norm = max(S.history.gap_norm,0);

    %% Final optimized parameters
    S.optimized = struct();
    S.optimized.density = extractScalar(txt, ...
        'Optimized density:\s*([0-9Ee+.\-]+)');
    S.optimized.young_kPa = extractScalar(txt, ...
        'Optimized young''s Modulus:\s*([0-9Ee+.\-]+)');
    S.optimized.poisson = extractScalar(txt, ...
        'Optimized poisson''s Ratio:\s*([0-9Ee+.\-]+)');
    S.optimized.thickness_mm = extractScalar(txt, ...
        'Optimized thikness\s*:\s*([0-9Ee+.\-]+)');
    S.optimized.width_mm = extractScalar(txt, ...
        'Optimized width\s*:\s*([0-9Ee+.\-]+)');
    S.optimized.length_mm = extractScalar(txt, ...
        'Optimized lenght\s*:\s*([0-9Ee+.\-]+)');

    %% Final summary lines
    tok = regexp(txt, ...
        ['initial E_pi:\s*' ...
         '([0-9Ee+.\-]+)\s+([0-9Ee+.\-]+)\s+([0-9Ee+.\-]+)\s+([0-9Ee+.\-]+)\s*' ...
         '->\s*resNorm=sum\(E_pi\)=:\s*([0-9Ee+.\-]+)\s*-\s*rms\(E_pi\):\s*([0-9Ee+.\-]+)'], ...
        'tokens', 'once');
    if ~isempty(tok)
        S.summary.initialErrPI = str2double(tok(1:4));
        S.summary.initialResNorm = str2double(tok{5});
        S.summary.initialRMS = str2double(tok{6});
    else
        S.summary.initialErrPI = nan(1,4);
        S.summary.initialResNorm = NaN;
        S.summary.initialRMS = NaN;
    end

    tok = regexp(txt, ...
        ['optimal E_pi:\s*' ...
         '([0-9Ee+.\-]+)\s+([0-9Ee+.\-]+)\s+([0-9Ee+.\-]+)\s+([0-9Ee+.\-]+)\s*' ...
         '->\s*resNorm=sum\(E_pi\)=:\s*([0-9Ee+.\-]+)\s*-\s*rms\(E_pi\):\s*([0-9Ee+.\-]+)'], ...
        'tokens', 'once');
    if ~isempty(tok)
        S.summary.optimalErrPI = str2double(tok(1:4));
        S.summary.optimalResNorm = str2double(tok{5});
        S.summary.optimalRMS = str2double(tok{6});
    else
        S.summary.optimalErrPI = nan(1,4);
        S.summary.optimalResNorm = NaN;
        S.summary.optimalRMS = NaN;
    end
end

function val = extractScalar(txt, pattern)
    tok = regexp(txt, pattern, 'tokens', 'once');
    if isempty(tok)
        val = NaN;
    else
        val = str2double(tok{1});
    end
end