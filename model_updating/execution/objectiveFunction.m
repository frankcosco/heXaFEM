
function [Error,Err_pi_reali,Err_lam,Err_V,f,n_DOFS] = objectiveFunction(x_sc_red,...
    EXA_type, EXA_order, Constraints, Elements, Nodes_, ...
     V_ref,lamda_ref,ind1cc,x_scale,logic2opt,res2use,ExpNodes)
    % x(1) = rho (density), x(2) = Ey (Young's modulus), x(3) = Ni (Poisson's ratio)

    x_sc(logic2opt)=x_sc_red;
    x_sc(~logic2opt)=1;

    x = x_sc.*x_scale;

    % Extract the material properties from optimization vector
    rho = (x(1));  % Density
    Ey = (x(2));   % Young's Modulus
    Ni = (x(3));   % Poisson's Ratio


    c  = x(4); %thickness


    a = x(5);
    b = x(6);

    Nodes = Nodes_;
    Nodes(:,4) = Nodes_(:,4)*(c/4);
    Nodes(:,2) = Nodes_(:,2)*(a/40);
    Nodes(:,3) = Nodes_(:,3)*(b/360);

    logicalIndConstrained = ~Constraints(:,1:3)';  % Define constraints
    [~,indexMorphingNodes] = OutputRequests(Nodes,logicalIndConstrained,ExpNodes);

    Nodes(indexMorphingNodes,2:3)=ExpNodes(:,2:3);

    % figure(1337),cla,
    % FEplot(Elements,Nodes,Constraints);
    % axis equal

    % Mesh domain with hexa elements
    % logicalIndConstrained = ~Constraints(:,1:3)'; % Constraints

    % Assembling Matrices
%     fprintf('Assembling matrices: ');
%     tic;
    isConsistent = true;  % Use lumped mass matrix
    % time_objf = tic;

      [M, K] = assmK_MFem(Constraints, Elements, Nodes, rho, Ey, Ni, EXA_type, EXA_order, isConsistent);

    % time_h = toc;
%     toc;
    % assembly_objf = toc(time_objf),

    % Number of modes to compute
    numberOfModes = 5;
    
    % Solve Eigenvalue Problem
    % fprintf('Solve Eigenvalue Problem: ');
%     tic;
    [V, D] = eigs(K, M, numberOfModes, 'smallestabs');%, 'Tolerance', 1e-15
%     time_h(2) = toc;

% eigenSolved_objf = toc(time_objf),
    
    % Extract eigenvalues and convert to frequencies
    Ds = diag(D);
%     lamda_estim = Ds;
    omega = sqrt(Ds);
    f = omega / (2 * pi);  % Natural frequencies in Hz
    
    % Remove 3rd mode (as per given requirement)
%     f_estim = f;
    V_all = zeros(numel(logicalIndConstrained),numberOfModes);
    V_all(logicalIndConstrained,:) = V;


    V_estim = V_all(ind1cc,[1 2 4 5]);
    lamda_estim = Ds([1 2 4 5],1);

 % % % %    figure(113),
 % % % % 
 % % % %    subplot(221),cla,
 % % % %    Nodes_def = Nodes;
 % % % %    Nodes_def(:,2:4)=Nodes_def(:,2:4)+1000*reshape(V_all(:,1)/norm(V_all(ind1cc,1)),3,[])';
 % % % %    FEplot(Elements,Nodes_def,Constraints);
 % % % % 
 % % % %    Nodes_ref = reshape(Nodes(:,2:4)',[],1);
 % % % %    Nodes_ref_def=reshape(Nodes_ref(ind1cc)+1000*V_ref(:,1)/norm(V_ref(:,1)),3,[])';
 % % % %    hold on,
 % % % %    plot3(Nodes_ref_def(:,1),Nodes_ref_def(:,2),Nodes_ref_def(:,3),'*r');
 % % % %    axis equal
 % % % % h_design1=gca;
 % % % % h_design1.View=[-89.7 -89];
 % % % % h_design1.Projection="orthographic";
 % % % % h_design1.CameraPosition=[-1444,-800,-870];
 % % % % h_design1.CameraTarget = [20 -180 -4];
 % % % % h_design1.CameraUpVector = [0 0 -1];
 % % % % h_design1.CameraViewAngle = 3.5;
 % % % % 
 % % % %        subplot(222),cla,
 % % % %         Nodes_def = Nodes;
 % % % %    Nodes_def(:,2:4)=Nodes_def(:,2:4)-1000*reshape(V_all(:,2)/norm(V_all(ind1cc,2)),3,[])';
 % % % %    FEplot(Elements,Nodes_def,Constraints);
 % % % % 
 % % % %    Nodes_ref = reshape(Nodes(:,2:4)',[],1);
 % % % %    Nodes_ref_def=reshape(Nodes_ref(ind1cc)+1000*V_ref(:,2)/norm(V_ref(:,2)),3,[])';
 % % % %    hold on,
 % % % %    plot3(Nodes_ref_def(:,1),Nodes_ref_def(:,2),Nodes_ref_def(:,3),'*r');
 % % % %    axis equal
 % % % % h_design1=gca;
 % % % % h_design1.View=[-89.7 -89];
 % % % % h_design1.Projection="orthographic";
 % % % % h_design1.CameraPosition=[-1444,-800,-870];
 % % % % h_design1.CameraTarget = [20 -180 -4];
 % % % % h_design1.CameraUpVector = [0 0 -1];
 % % % % h_design1.CameraViewAngle = 3.5;
 % % % % 
 % % % %        subplot(223),cla,
 % % % %         Nodes_def = Nodes;
 % % % %    Nodes_def(:,2:4)=Nodes_def(:,2:4)-100*reshape(V_all(:,4)/norm(V_all(ind1cc,4)),3,[])';
 % % % %    FEplot(Elements,Nodes_def,Constraints);
 % % % % 
 % % % %    Nodes_ref = reshape(Nodes(:,2:4)',[],1);
 % % % %    Nodes_ref_def=reshape(Nodes_ref(ind1cc)+100*V_ref(:,3)/norm(V_ref(:,3)),3,[])';
 % % % %    hold on,
 % % % %    plot3(Nodes_ref_def(:,1),Nodes_ref_def(:,2),Nodes_ref_def(:,3),'*r');
 % % % %    axis equal
 % % % % h_design1=gca;
 % % % % h_design1.View=[-89.7 -89];
 % % % % h_design1.Projection="orthographic";
 % % % % h_design1.CameraPosition=[-1444,-800,-870];
 % % % % h_design1.CameraTarget = [20 -180 -4];
 % % % % h_design1.CameraUpVector = [0 0 -1];
 % % % % h_design1.CameraViewAngle = 3.5;
 % % % % 
 % % % %        subplot(224),cla,
 % % % %  Nodes_def = Nodes;
 % % % %    Nodes_def(:,2:4)=Nodes_def(:,2:4)+100*reshape(V_all(:,5)/norm(V_all(ind1cc,5)),3,[])';
 % % % %    FEplot(Elements,Nodes_def,Constraints);
 % % % % 
 % % % %    Nodes_ref = reshape(Nodes(:,2:4)',[],1);
 % % % %    Nodes_ref_def=reshape(Nodes_ref(ind1cc)+100*V_ref(:,4)/norm(V_ref(:,4)),3,[])';
 % % % %    hold on,
 % % % %    plot3(Nodes_ref_def(:,1),Nodes_ref_def(:,2),Nodes_ref_def(:,3),'*r');
 % % % %    axis equal
 % % % % h_design1=gca;
 % % % % h_design1.View=[-89.7 -89];
 % % % % h_design1.Projection="orthographic";
 % % % % h_design1.CameraPosition=[-1444,-800,-870];
 % % % % h_design1.CameraTarget = [20 -180 -4];
 % % % % h_design1.CameraUpVector = [0 0 -1];
 % % % % h_design1.CameraViewAngle = 3.5;

    % 
    % V_estim = V_all(ind1cc,[1 2 5]);
    % lamda_estim = Ds([1 2 5],1);
    % 
    % V_ref = V_ref(:,[1 2 4]);
    % lamda_ref = lamda_ref([1 2 4],1);
%     [~, ind1nodes]=OutputRequests(Nodes,logicalIndConstrained,Nodes_ref);

%     ind1cc = reshape([3*(ind1nodes-1)+1;
%           3*(ind1nodes-1)+2;
%           3*(ind1nodes-1)+3],[],1);

    
    % Error estimation using MAC (Modal Assurance Criterion)
    [Err_pi_reali,Err_lam,Err_V,ErrorsV] = MAC_Reali_Diff_Modal_updating(...
            V_estim,lamda_estim,...
            V_ref,lamda_ref);
    
    % Output the error for optimization
    % Error = sqrt(MAC_reali); % Prova 1

    % Error = Err_pi_reali; % Prova 1

    % Error = max(MAC_reali);% Prova 2
    % Error = [Err_lam;Err_V([1 2 4])];


if size(res2use,1)>1
    reslam2use = res2use(1,:);
    resV2use = res2use(2,:);
else
    reslam2use = res2use;
    resV2use = res2use;
end


    Error = [Err_lam(reslam2use);reshape(ErrorsV(:,resV2use),[],1)];%blocked Exa27
    % Error = sqrt(Err_pi_reali(res2use));
    % Error = [Err_lam(res2use)];
    % Error = [Err_lam(selected);Err_V(selected) ];
    if nargout>5
        n_DOFS = size(M,1) ;
    end
end

