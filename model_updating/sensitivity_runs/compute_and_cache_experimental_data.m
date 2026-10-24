%% compute_and_cache_experimental_data.m
% Runs the experimental-data pipeline (UFF read + DIC realignment) EXACTLY ONCE
% and caches ExpMV, ExpLambda, ExpMNodes, ExpModal_freqs to a .mat file.
%
% WHY THIS EXISTS: riallineaDICtests() internally runs its own lsqnonlin fit
% (Levenberg-Marquardt, in fineTuneHomography/residui_vincolati) to refine the
% pixel->physical homography. That inner fit is a separate nonlinear least-squares
% problem from the main model-updating optimization, and re-running it from
% scratch on every multi-start trial would add its own run-to-run numerical noise
% on top of the thing we actually want to study (sensitivity of the OUTER
% model-updating optimizer to different starting points p0). Computing it once
% and reusing the cached result for every trial isolates that effect.
%
% Run this ONCE from model_updating/sensitivity_runs/, then run multi_start_hexa8.m.

clear
addpath(genpath('..\..\'))  % fe_core/, experimental_data/, model_updating/*

%% Experimental data (identical logic to Modal_updating_clean_LM.m lines ~13-58)
[ExpV,ExpModal_freqs,ExpNodes] = readuff_modalshapes('real_modeset.unv'); %#ok<ASGLU>
ExpModal_freqs = [25.3118 157.6535 410.1009 440.2394];

ExpMNodes = riallineaDICtests();

ExpV2 = reshape(ExpV(:,2),3,[]);
ExpV(:,2) = -ExpV2(:);

ExpV3 = reshape(ExpV(:,3),3,[]);
ExpV(:,3) = -ExpV3(:);

ExpMV = ExpV;

ExpLambda = (2*pi * ExpModal_freqs').^2;

fprintf('Cached experimental data: %d nodes, %d mode-shape DOFs, freqs = %s\n', ...
    size(ExpMNodes,1), size(ExpMV,1), mat2str(ExpModal_freqs,6));

save('aligned_experimental_data.mat', 'ExpMV', 'ExpLambda', 'ExpMNodes', 'ExpModal_freqs', 'ExpNodes');
fprintf('Saved aligned_experimental_data.mat in %s\n', pwd);
