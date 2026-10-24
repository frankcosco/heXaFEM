%% CLEAN REFERENCE COPY -- derived from __final__ModelUpdatingBeam__AnnaPHD/Modal_updating_clean_LM.m
% Correction applied vs. the live script in dev2 (as of 23/10/2026):
%   var2optimize was [true true true true false false] (thickness free) in the
%   live script, but ALL saved result files in revisedFrequencies/materialOnly/
%   report the optimized thickness identical to the initial value (4.0000 mm)
%   across all 12 mesh/formulation cases -- proof that the historical "material
%   only" runs (the ones consistent with the paper's Fig. 8, which reports only
%   E, rho, nu) were run with thickness FIXED. Corrected here (line below, in the
%   k-loop) to [true true true false false false] (rho, E, nu free; thickness/
%   width/length fixed) to match the paper. VALIDATE by re-running k=2
%   (Hexa8-Int8, 26x234x2) and comparing against
%   revisedFrequencies/materialOnly/MU-Singlerun-TRRBeamHexa8-Int8_26x234x2.txt
%   (rho=2.8317e-06, E=6.4781e+07 kPa, nu=0.3159) before trusting this for new runs.
%% ---------------------------------------------------------------------------


%% Clear previous data
clear
addpath(genpath('..\..\'))  % pulls in fe_core/, experimental_data/, model_updating/post_processing/ (now 2 levels up to dev_SR_github root)

%%
% Beam dimensions in mm
a = 40; % Length (L)
b = 360;  % Width (W)
c = 4;   % Thickness (T)

%% Experimental data
[ExpV,ExpModal_freqs,ExpNodes] = readuff_modalshapes('real_modeset.unv'); %#ok<ASGLU>
ExpModal_freqs = [25.3118 157.6535 410.1009 440.2394];

% ExpNodes=ExpNodes';
% ExpNodes(:,2:4)=ExpNodes(:,2:4)*1000;
% ExpNodes(:,2) = ExpNodes(:,2)+a/2;
% ExpNodes(:,3) = ExpNodes(:,3)-b/2;
% ExpNodes(:,4)= -c;
ExpMNodes=riallineaDICtests();

% 
% 

m3_torsional = reshape(ExpV(:,3),3,[])';

% figure(13)
% subplot(2,2,1),plot(m3_torsional(:,1),m3_torsional(:,2),'*'), axis equal
% subplot(2,2,2),plot(m3_torsional(:,2),m3_torsional(:,3),'*'), 
% subplot(2,2,3),plot(m3_torsional(:,1),m3_torsional(:,3),'*'),

% 
DT = delaunayTriangulation(ExpMNodes(:,2:3));
figure(20),cla,
patch('Faces',DT.ConnectivityList,'Vertices',ExpMNodes(:,2:4)+m3_torsional*10000,...
    'FaceColor', 'interp','FaceVertexCData',m3_torsional(:,3),'LineStyle','none')
axis equal
hold on
patch('Faces',DT.ConnectivityList,'Vertices',ExpMNodes(:,2:4),...
    'FaceColor', 'g','LineStyle','none')

% ExpNodes(:,4) = -c;
%%
% % 
% % 
% % oneSideNodes = ExpNodes(:,2)>30 | ExpNodes(:,2)<10;
% % 
ExpV2=reshape(ExpV(:,2),3,[]);
% ExpV2(3,:)=-ExpV2(3,:);
ExpV(:,2)=-ExpV2(:);

ExpV3=reshape(ExpV(:,3),3,[]);
% ExpV3(3,:)=-ExpV3(3,:);
ExpV(:,3)=-ExpV3(:);

% ExpMV=ExpV(reshape(repmat(oneSideNodes,1,3)',[],1),:);
ExpMV=ExpV;


% ExpMNodes = ExpNodes(oneSideNodes,:);
ExpMNodes=riallineaDICtests();
% PgridCorrected(:,3)=-c;

% 
% ExpMNodes(:,2) = ExpMNodes(:,2)-ExpMNodes(:,3)./360*2.4 ;
% ExpMNodes(:,2) = ExpMNodes(:,2)-ExpMNodes(:,2)./40*1.6 ;
% ExpMNodes(:,3) = ExpMNodes(:,3)+1.15 +(ExpMNodes(:,2)-a/2)/(a/2).*(ExpMNodes(:,3)+b/2)/(b/2)*1.0;
% % % 


[ConstrFE,ElemFE,NodesFE] = brick_hexa8_mesher_advanced( ...
    1,...
    a,b,-c, 13,117,1,...
    0,-b,0);
logIndConstFE  = ~ConstrFE(:,1:3)';
 [~,indexMorphingNodes] = OutputRequests(NodesFE,logIndConstFE,ExpMNodes);

NodesFE(indexMorphingNodes,2:4)=ExpMNodes(:,2:4);
NodesRef = NodesFE(indexMorphingNodes,:);
% ExpMNodes(:,2:4)=NodesRef(:,2:4);



%  plot experiments
figure(14),cla
title('design mesh prototipe')
  FEplot(ElemFE,NodesFE,ConstrFE); %model

 hold on
scatter3(ExpNodes(:,2),ExpNodes(:,3),ExpNodes(:,4),'*b');
scatter3(ExpMNodes(:,2),ExpMNodes(:,3),ExpMNodes(:,4),'ob');
axis equal


% 
% xlabel('x')
% ylabel('y')
% zlabel('z')





 axis equal


 h_design=gca;
 h_design.View=[-89.7 -89];
 h_design.Projection="orthographic";
 h_design.CameraPosition=[-1444,-800,-870];
 h_design.CameraTarget = [20 -180 -4];
 h_design.CameraUpVector = [0 0 -1];
 h_design.CameraViewAngle = 3.5;
 h_deesign.CameraUpVector=[0 0 -1];

%%
ExpLambda  = (2*pi * ExpModal_freqs').^2;
% ExpV = ExpV
%% Material and FE properties (initial guess for aluminum)
RHO_vec = 2.70e-06; % Density (kg/mm³)
Ey_vec = 69580000;   % Young's Modulus (kPa)
Ni_vec = 0.33;       % Poisson's Ratio

%% DoE
nW = 13;
nL = 117;
nH = 1;

% W_values  = [10  20 10  30 10  40  20   50  60  30  20   ];
% L_values  = [90 180 90 270 90 360 180  450 540 270 180   ];
% H_values  = [1    2 1    3  1   4   2    5   6   3   2   ];

%             1   2  3   4  5   6  7   8  9  10  11  12  13
S_values  = [ 1   2  1  1   3  1   4  2   2  5   6   3   3   2   ];
EXA_types = [ 8   8  20 27  8 64   8  20 27  8   8  27   20  64   ];
EXA_orders= [ 8   8  27 27  8 64   8  27 27  8   8  27   27  64   ];
% W_values  = [            10               20  30  50  60 ];
% L_values  = [            90              180 270 450 540 ];
% H_values  = [             1                2   3   5   6 ];
% EXA_types = [            64               64  27   8   8 ];
% EXA_orders= [            64               64  27   8   8 ];



numberOfModes = 5;  % Number of modes for modal updating
numberofH = numel(S_values);  % Mesh levels
% lamda = zeros(numberOfModes, numberofH);  % Eigenvalues
% f_estim = zeros(numberOfModes, numberofH);  % Estimated frequencies
% time_h = zeros(numberofH, 4);  % Time statistics
% Ndim = zeros(1, numberofH);  % Total DOFs
% Edim = zeros(1, numberofH);  % Number of elements



x0 =      [RHO_vec, Ey_vec, Ni_vec,4.0,  40.0, 360];
% x0= [2.76635128057767e-06,6.79008602719778e+07,0.32931001949276,3.873,40,360 ];
% x0 =      [2.79e-06, 68.8e6, 0.3358,3.85, 40.0, 360 ];%opt trials
% lower_b = [2.68e-06, 66.5e6, 0.29,  3.6,   36, 360];  % Lower bound
% upper_b = [2.82e-06, 70e6,   0.37,  4.1,   42, 360];  % Upper bound
lower_b = [2.60e-06, 64e6,  0.30,  3.5,   36, 360];  % Lower bound
upper_b = [2.85e-06, 74e6,  0.38,  4.5,   42, 360];  % Upper bound


% x0 = [RHO_vec, Ey_vec, Ni_vec];
% lower_b = [2.640e-06, 67e6, 0.31];  % Lower bound
% upper_b = [2.810e-06, 74e6, 0.34];  % Upper bound
    % 
    % lower_b = [2e-06, 30e6, 0.25, 3.5];  % Lower bound
    % upper_b = [5e-06, 100e6, 0.4, 4.5];  % Upper bound


% x_optimal=x0;


x0_scaled = x0./x0;
lower_b_sc = lower_b./x0;
upper_b_sc = upper_b./x0;


%%%            %1st bending   %2nd bending  %1st torsional  %3rd bending 
res2use      = [true             true         true                true ;
                true             true         true                true ];

for k =[2] % VALIDATION: Hexa8-Int8, 26x234x2 -- compare vs materialOnly reference. Restore [11:14] or others to reproduce other figures.
    % Mesh settings
    na = S_values(k) *nW;%W_values(k);
    nb = S_values(k) *nL;%L_values(k);
    nc = S_values(k) *nH;%H_values(k);

    % Element type and order
    EXA_type  = EXA_types(k);
    EXA_order = EXA_orders(k);


    namefile_str = sprintf("MU-Singlerun-TRRBeamHexa%i-Int%i_%ix%ix%i",EXA_type,EXA_order,na,nb,nc);

    logsDir = fullfile("..","logs");
    if ~exist(logsDir,"dir"); mkdir(logsDir); end
    namefile_log = fullfile(logsDir,sprintf("%s.txt",namefile_str));



    %% Mesh domain with hexa elements (Initial discretization)
    [Constraints, Elements, Nodes] = brick_hexa_mesher(EXA_type, 1, a, b, -c, na, nb, nc);  % Mesh generator
    Nodes(:,2:4) = Nodes(:,2:4)    + repmat([-0 -b  -0*c],size(Nodes,1),1);

    %Constrains 
    % Constraints = false(NN,3);
    %Choose where constrain
    % indClamp = Nodes(:,3)==0;
    % Constraints(Nodes(:,3)==0,:)=true;
    logicalIndConstrained = ~Constraints(:,1:3)';  % Define constraints


    %%
    [~,indexMorphingNodes] = OutputRequests(Nodes,logicalIndConstrained,NodesRef);
    indexOutputDofs = reshape([3*(indexMorphingNodes-1)+1;
        3*(indexMorphingNodes-1)+2;
        3*(indexMorphingNodes-1)+3],[],1);


    %% PLOT



 % % % % %    figure(111),cla,
 % % % % %    P1=plot3(Nodes(:,2),Nodes(:,3),Nodes(:,4),'*b');
 % % % % %    FEplot(Elements,Nodes,Constraints); %model
 % % % % %    axis equal
 % % % % % 
 % % % % %    hold on
 % % % % %    scatter3(ExpMNodes(:,2),ExpMNodes(:,3),ExpMNodes(:,4),'*m');
 % % % % % 
 % % % % %     h_design=gca;
 % % % % % h_design.View=[-89.7 -89];
 % % % % % h_design.Projection="orthographic";
 % % % % % h_design.CameraPosition=[-1444,-800,-870];
 % % % % % h_design.CameraTarget = [20 -180 -4];
 % % % % % h_design.CameraUpVector = [0 0 -1];
 % % % % % h_design.CameraViewAngle = 3.5;
 % % % % % 
 % % % % % 
 % % % % % 
 % % % % % 
 % % % % %    drawnow;
 % % % % %        h_cube=gca;
%%
    if exist( namefile_log,"file") 
        [path, file,ext ]=fileparts(namefile_log);
        delete(namefile_log);
    end    
    diary(namefile_log);
    diary("on");
    fprintf("\n\n######################################## trust-region-reflective\n"); ...
    
    fprintf("Model Updating Problem.\n"); ...
    fprintf("Beam Dimensions:    %i x %i x %i\n",a,b,c);
    fprintf("Discretized with: %ix%ix%i Hexa%i-Int%i elements\n",na,nb,nc,EXA_type,EXA_order);
   
    
    %% Initial parameters (x0: [density, Young's modulus, Poisson's ratio])
    % tic
    x0_scaled = x0./x0;
%%%                   rho   Ey    Ni    c    a      b
var2optimize =       [true true true false false false];

    [Error_ini,Err_pi_reali_ini,Err_lam_ini,Err_V_ini ,f_ini,n_DOFS] = objectiveFunction(x0_scaled(var2optimize),...
        EXA_type, EXA_order, Constraints, Elements, Nodes,...
        ExpMV,  ExpLambda,indexOutputDofs,x0,var2optimize,res2use,ExpMNodes);
    % toc

    fprintf("free DOFs: %i - fixed DOFs: %i - total DOFs: %i \n",n_DOFS, sum(Constraints,"all"),n_DOFS + sum(Constraints,"all"));
    fprintf('First Objective Function Evaluation: ');
   

    fprintf(' initial E_pi: %1.4f %1.4f %1.4f %1.4f -> resNorm=sum(E_pi)=:%1.4f - rms(E_pi):%1.4f\n', Err_pi_reali_ini,sum(Err_pi_reali_ini),rms(Err_pi_reali_ini));
    

    started_optimization = tic;
%% lsqnonlin (least squares non-linear optimization) - step 1
    % NOTE: 'ScaleProblem' was removed here -- it is a Levenberg-Marquardt-only
    % option (MathWorks docs group it under the LM algorithm table, not under
    % Trust-Region-Reflective), so it had no effect with Algorithm set to
    % 'trust-region-reflective' below; it was dead/inert in the original script.
    options = optimoptions('lsqnonlin', 'Display', 'iter-detailed', ...
        'StepTolerance',1e-4, 'Algorithm','trust-region-reflective', ...
        'UseParallel',true, 'FunctionTolerance',1e-5, ...
        'OutputFcn',@(x,optimValues,state) outFun(x,optimValues,state,started_optimization,res2use,namefile_str,x0,var2optimize) ...
        );

    if k>10
        options.UseParallel=false;
    end
        % 
        % 'OptimalityTolerance',1e-9, ...,'StepTolerance',1e-6,'FunctionTolerance',1e-9

% x_sc_optimal_red=x0_scaled(var2optimize);
    [x_sc_optimal_red,resnorm, residual, exitflag, output] = ...
        lsqnonlin(@(x) objectiveFunction(x ,...
        EXA_type, EXA_order, Constraints, Elements, Nodes,...
        ExpMV,  ExpLambda,indexOutputDofs,x0,var2optimize,res2use,ExpMNodes), ...
        x0_scaled(var2optimize) , lower_b_sc(var2optimize), upper_b_sc(var2optimize), options);

% % %  x0_scaled(var2optimize) = x_sc_optimal_red;
% % % 
% % % %%%                   rho   Ey    Ni   c    a      b
% % % var2optimize =       [true true true true false false];
% % % 
% % % 
% % % options.OutputFcn = @(x,optimValues,state) outFun(x,optimValues,state,started_optimization,res2use,namefile_str,x0,var2optimize);
% % % % options.StepTolerance = 1e-3;
% % % % options.Algorithm = "trust-region-reflective";
% % % 
% % % % x_sc_optimal_red=x0_scaled(var2optimize);
% % %     [x_sc_optimal_red, resnorm, residual, exitflag, output] = ...
% % %         lsqnonlin(@(x) objectiveFunction(x ,...
% % %         EXA_type, EXA_order, Constraints, Elements, Nodes,...
% % %         ExpMV,  ExpLambda,indexOutputDofs,x0,var2optimize,res2use,ExpMNodes), ...
% % %         x0_scaled(var2optimize) , lower_b_sc(var2optimize), upper_b_sc(var2optimize), options);


    elapsedTime = toc(started_optimization);





    %% fminunc
    % options_fminunc = optimoptions("fminunc","Algorithm","quasi-newton","Display","iter-detailed");
    % [x_sc_optimal, fval,exitflag,output] = ...
    %     fminunc(@(x) sum(abs(objectiveFunction(x,...
    %     EXA_type, EXA_order, Constraints, Elements, Nodes,...
    %     ExpV,  ExpLambda,indexOutputDofs,x0))), ...
    %     x0_scaled, optimset('Display','iter','maxIter',10));
    
    %% fminsearch 
    % 
    % [x_sc_optimal, fval,exitflag,output] = ...
    %     fminsearch(@(x) sum(abs(objectiveFunction(x,...
    %     EXA_type, EXA_order, Constraints, Elements, Nodes,...
    %     ExpV,  ExpLambda,indexOutputDofs,x0))), ...
    %     x0_scaled, optimset('Display','iter','maxIter',10));

    % 
    % 
    % fprintf('Optimized density: %1.14e [kg/mm³]\n', x_optimal(1));
    % fprintf('Optimized young''s Modulus: %1.14e [kPa]\n', x_optimal(2));
    % fprintf('Optimized poisson''s Ratio: %.14f \n', x_optimal(3));
  
    % x_sc_optimal_=[x_sc_optimal x0_scaled(3:4)];
    x_sc_optimal = 0*var2optimize;
    x_sc_optimal(var2optimize) = x_sc_optimal_red;
    x_sc_optimal(~var2optimize)=1;
    x_optimal = x0.*x_sc_optimal;

    fprintf('Optimized density: %1.4e [kg/mm³]\n', x_optimal(1));
    fprintf('Optimized young''s Modulus: %1.4e [kPa]\n', x_optimal(2));
    fprintf('Optimized poisson''s Ratio: %.4f \n', x_optimal(3));
    fprintf('Optimized thikness        : %.4f [mm]\n', x_optimal(4));
    fprintf('Optimized width           : %.4f [mm]\n', x_optimal(5));
    fprintf('Optimized lenght          : %.4f [mm]\n', x_optimal(6));



    % 
    % 


    % Optimized material properties
    rho_optimal = x_optimal(1);
    Ey_optimal = x_optimal(2);
    Ni_optimal = x_optimal(3);
    c_optimal  = x_optimal(4);

    fprintf('Model updating  %.2f seconds.\n', elapsedTime);




    %%

    [Error_opt,Err_pi_reali_opt,Err_lam_opt,Err_V_opt ,f_opt] = objectiveFunction(x_sc_optimal_red,...
        EXA_type, EXA_order, Constraints, Elements, Nodes,...
        ExpMV,  ExpLambda,indexOutputDofs,x0,var2optimize,res2use,ExpMNodes);

    fprintf(' initial E_pi: %1.4f %1.4f %1.4f %1.4f -> resNorm=sum(E_pi)=:%1.4f - rms(E_pi):%1.4f\n', Err_pi_reali_ini,sum(Err_pi_reali_ini),rms(Err_pi_reali_ini));
    fprintf(' optimal E_pi: %1.4f %1.4f %1.4f %1.4f -> resNorm=sum(E_pi)=:%1.4f - rms(E_pi):%1.4f\n', Err_pi_reali_opt,sum(Err_pi_reali_opt),rms(Err_pi_reali_opt));

    diary("off");

    % namefile = sprintf("%s.mat",namefile_str);
    % save(namefile,...
    %     "Err_pi_reali_ini",  "Err_V_ini", "Err_lam_ini","f_ini", ...
    %     "Err_pi_reali_opt",  "Err_V_opt", "Err_lam_opt","f_opt", ...
    %     "resnorm", "output", "elapsedTime", ...
    %     "rho_optimal", "Ey_optimal","Ni_optimal" ...
    %      ,"c_optimal"...
    %     );
 end
%%
% % ExpModal_freqs
% for i=1:numberOfModes
% fprintf('mode %d: Initial = %f [Hz]; Optimal = %f [Hz]\n', i,f_ini(i),f_opt(i));
% end
%
% for i=1:4
% fprintf('Error mode %d: Initial = %f ; Optimal = %f \n', i,Error_ini(i),Error_opt(i));
% end
%%%%%
%save Error_ini,Err_V_ini, Err_lam_ini,f_ini, Error_opt,Err_V_opt,
%Err_lam_opt,f_opt, resnorm, output, elapsedTime, rho_optimal, Ey_optimal,
%Ni_optimal