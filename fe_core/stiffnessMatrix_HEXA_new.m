function [K,M] = stiffnessMatrix_HEXA_new (...
    Nodes_el,ey,nu,rho, ...
    X,W,N_,dN_)
% STIFFNESSMATRIX_HEXA con MS 24x24, input(E,nu,nodes)for hexahedral finite
% element;
% ey - modulo elastico;
% nu - coefficiente di Poisson;
% Nodes_el [x1 ..x8; y1..y8;z1..z8]
%          [x1 ..x20; y1..y20;z1..z20]
%          [x1 ..x27; y1..y27;z1..z27]
%          [x1 ..x32; y1..y32;z1..z32]

% Matrice delle costanti elastiche
D= ey./((1+nu)*(1-2*nu))*...
    [1-nu   nu   nu     0          0          0;
       nu 1-nu   nu     0          0          0;...
       nu   nu 1-nu     0          0          0;
      0    0    0   (1-2*nu)/2     0          0;
      0    0    0       0      (1-2*nu)/2     0;...
      0    0    0       0          0      (1-2*nu)/2];
n_nodes = size(Nodes_el,1);
n_DOF=n_nodes*3;
% 
% Dxx = ey./((1+nu)*(1-2*nu))*...
%     [1-nu   nu   nu   
%        nu 1-nu   nu   
%        nu   nu 1-nu ];
% Dxy = ey./((1+nu)*(1-2*nu))*...
%     [(1-2*nu)/2     0          0;
%           0      (1-2*nu)/2     0;...
%           0          0      (1-2*nu)/2];

%% Integrazione numerica  con punti Gauss
%1x1x1 1 punto
%2x2x2 8 punti
% 14 punti
%3x3x3x 27 punti
%4x4x4 64 puni

% Gauss points coordinates on each direction
% switch (hexa_order)
%      case 1
%         [X,W] = HEXA_integration_1();
%      case 8
%         [X,W]= HEXA_integration_8();
%      case 27
%         [X,W]= HEXA_integration_27();
%      case 14
%        [X,W]= HEXA_integration_14();
%      case 64
%     [X,W]= HEXA_integration_64();
% end
%% ADD

% Number of Gauss Points
Ng=size(X,1);

coordinates = Nodes_el;

% 
% % Preallocate memory for stiffness matrix
% K = zeros (n_DOF,n_DOF,Ng );
% 
% % Preallocate memory for mass matrix
% M = zeros (n_DOF,n_DOF,Ng);

% Preallocate memory for stiffness matrix
K = zeros (n_DOF,n_DOF);

% Preallocate memory for mass matrix
M = zeros (n_DOF,n_DOF);


% Jh_test=zeros(3,3,Ng);
for ig=1:Ng
    % xi1=X(ig,1);
    % xi2=X(ig,2);
    % xi3=X(ig,3);
    wg=W(ig);

    %% Funzioni di forma
    % File Shape_funct - Ruhr Universität Bochum_Shape Functions
    % generation,requirements, etc. PP-38-42)
    % switch (hexa_type)
    % case 8
    %     [N,dShape] = HEXA_shapes_8(xi1,xi2,xi3);
    % case 20
    %     [N,dShape] = HEXA_shapes_20(xi1,xi2,xi3);
    % case 27
    %     [N,dShape] = HEXA_shapes_27(xi1,xi2,xi3);
    % case 32
    %     [N,dShape] = HEXA_shapes_32(xi1,xi2,xi3);
    % case 64
    %     [N,dShape] = HEXA_shapes_64(xi1,xi2,xi3);
    %
    % end
    N_ig = N_(:,:,ig);
    dShape = dN_(:,:,ig);

    %% Сalcolo delle matrici di massa e rigidezza
    %Compute Jacobian matrix
    %             JacobianMatrix = dShape *coordinates;
    %             dN = inv(dShape *coordinates ) * dShape ;
    J = dShape *coordinates;

    %    Jh_test(:,:,ig)=J;
    invJ = inv(J);
    dN = invJ*dShape ;   %#ok<MINV>

    wg_detJ = wg*det(J);

    % % Preallocate memory for B-Operator
    B   = zeros(6,n_DOF);
    % Bxx = zeros(3,n_DOF);
    % Bxy = zeros(3,n_DOF);
    % 
    % B(1,1:3:end)=dN(1,:);
    for j=0:n_nodes-1

         % Construct first three rows
        B(1,3*j+1) = dN(1,j+1);
        B(2,3*j+2) = dN(2,j+1);
        B(3,3*j+3) = dN(3,j+1);

        % Construct fourth row
        B(4,3*j+1) = dN(2,j+1);
        B(4,3*j+2) = dN(1,j+1);
        % Construct fifth row
        B(5,3*j+3) = dN(2,j+1);
        B(5,3*j+2) = dN(3,j+1);
        % Construct sixth row
        B(6,3*j+1) = dN(3,j+1);
        B(6,3*j+3) = dN(1,j+1);
        % 
        % % Construct first three rows
        % Bxx(1,3*j+1) = dN(1,j+1);
        % Bxx(2,3*j+2) = dN(2,j+1);
        % Bxx(3,3*j+3) = dN(3,j+1);
        % 
        % % Construct fourth row
        % Bxy(1,3*j+1) = dN(2,j+1);
        % Bxy(1,3*j+2) = dN(1,j+1);
        % % Construct fifth row
        % Bxy(2,3*j+3) = dN(2,j+1);
        % Bxy(2,3*j+2) = dN(3,j+1);
        % % Construct sixth row
        % Bxy(3,3*j+1) = dN(3,j+1);
        % Bxy(3,3*j+3) = dN(1,j+1);
    end
    %             B = sparse(B);
    
    % Add to stiffness matrix
    % Kg = wg_detJ*B'*D*B;
    K=K+wg_detJ*B'*D*B;  



    % 
    % 
    % for i=0:n_nodes-1
    %     % dNxi = dN(1,i+1);
    %     % dNyi = dN(2,i+1);
    %     % dNzi = dN(3,i+1);
    %     for j=0:n_nodes-1
    %         % dNxj = dN(1,j+1);
    %         % dNyj = dN(1,j+1);
    %         % dNzj = dN(1,j+1);
    % 
    %         K(3*i+1,3*j+1,ig)= wg_detJ*dNxi*dNxj*(1-nu) + wg_detJ*(dNyi*dNyj+dNzi*dNzj)*(1-2*nu)/2;
    %         K(3*i+2,3*j+1,ig)= wg_detJ*dNyi*dNxj*nu     + wg_detJ* dNxi*dNyj*(1-2*nu)/2;
    %         K(3*i+3,3*j+1,ig)= wg_detJ*dNzi*dNxj*nu     + wg_detJ* dNxi*dNzj*(1-2*nu)/2;
    %         K(3*i+1,3*j+2,ig)= wg_detJ*dNxi*dNyj*nu     + wg_detJ* dNyi*dNxj*(1-2*nu)/2;
    %         K(3*i+2,3*j+2,ig)= wg_detJ*dNyi*dNyj*(1-nu) + wg_detJ*(dNxi*dNxj+dNzi*dNzj)*(1-2*nu)/2;
    %         K(3*i+3,3*j+2,ig)= wg_detJ*dNzi*dNyj*nu     + wg_detJ* dNyi*dNzj*(1-2*nu)/2;
    %         K(3*i+1,3*j+3,ig)= wg_detJ*dNxi*dNzj*nu     + wg_detJ* dNzi*dNxj*(1-2*nu)/2;
    %         K(3*i+2,3*j+3,ig)= wg_detJ*dNyi*dNzj*nu     + wg_detJ* dNzi*dNyj*(1-2*nu)/2;
    %         K(3*i+3,3*j+3,ig)= wg_detJ*dNzi*dNzj*(1-nu) + wg_detJ*(dNxi*dNxj+dNyi*dNyj)*(1-2*nu)/2;
    % 
    %     end
    % end

     % K = K + (Kg);
    % K = K + (Kg+Kg')/2 ;
    % M(:,:,ig) = wg_detJ*(N_ig'*N_ig)*rho;
    % Mg = 
    
    M=  M  +wg_detJ*(N_ig'*N_ig)*rho;

    %         end
    %     end
end
% K=sum(K,3);
% M=sum(M,3);

%std(Jh_test,0,3)