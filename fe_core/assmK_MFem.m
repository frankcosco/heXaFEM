function [m_mat,k_mat]=assmK_MFem(...
    Constraints,...
    Elements,Nodes,...
    RHO_vec,Ey_vec,Ni_vec,...
    hexa_type, hexa_order , ...
    isConsistent)
%% 
%  Elements - Ne x (3 + 8) matrix
%    [ EID PID TypeID G1 G2 G3 G4 [g5 g6 g7 g8 .... gN ]]
%  EID - element ID
%  PID - material property ID 
%  TypeID = {0- TETRA4,  1- TETRA10, 2-TETRA20, 
%             3-HEXA8, 
%             4-HEXA20(QUAD), 5-HEXA27(QUAD),
%             6-HEXA32(CUBIC). 7-HEXA64(CUBIC) .... }
%  Rho_vec - [rho(1) ... rho(np)],  np - number of material properties
%  Ey_vec
%  Ni_vec

% se il parametro "isConsistent" non 'e passato come argomento, usiamo
% isConsistent=false per avere lumped
if nargin<9
    isConsistent=false;
end
Ndof_node=3;
N_node = size(Nodes,1);
N_element=size(Elements,1);

%% Hexa order
switch(hexa_type)
    case 8
    N_dofs_x_element=24;
    case 20
    N_dofs_x_element=60;  
    case 27
    N_dofs_x_element=81;
    case 32
    N_dofs_x_element=96;
    case 64
    N_dofs_x_element=192;      
end

% dfg
i_mat  =zeros(N_element,N_dofs_x_element*N_dofs_x_element);
j_mat  =zeros(N_element,N_dofs_x_element*N_dofs_x_element);
v_matK =zeros(N_element,N_dofs_x_element*N_dofs_x_element);
v_matM =zeros(N_element,N_dofs_x_element*N_dofs_x_element);


%% Initialize Mesh Invariants
% Gauss points coordinates and weights
switch (hexa_order)
     case 1
        [X,W] = HEXA_integration_1();        
     case 8
        [X,W]= HEXA_integration_8();
     case 27
        [X,W]= HEXA_integration_27();    
     case 14
       [X,W]= HEXA_integration_14(); 
     case 64
       [X,W]= HEXA_integration_64(); 
end

% Loop over each Gauss point
Ng=size(X,1);
% Jh_test=zeros(3,3,Ng);
N = zeros(3,N_dofs_x_element,Ng);
dShape = zeros(3,N_dofs_x_element/Ndof_node,Ng);

shape_times=tic;
for ig=1:Ng
    xi1=X(ig,1);
    xi2=X(ig,2);
    xi3=X(ig,3);
    % wg=W(ig);

    switch (hexa_type)
        case 8
            [N(:,:,ig),dShape(:,:,ig)] = HEXA_shapes_8(xi1,xi2,xi3);
        case 20
            [N(:,:,ig),dShape(:,:,ig)] = HEXA_shapes_20(xi1,xi2,xi3);
        case 27
            [N(:,:,ig),dShape(:,:,ig)] = HEXA_shapes_27(xi1,xi2,xi3);
        case 32
            [N(:,:,ig),dShape(:,:,ig)] = HEXA_shapes_32(xi1,xi2,xi3);
        case 64
            [N(:,:,ig),dShape(:,:,ig)] = HEXA_shapes_64(xi1,xi2,xi3);

    end


end

 % switch (hexa_type)
 %     case 64
 %         [N,~,dShape]=HEXA_shapes_64_newww(X(:,1),X(:,2),X(:,3));
 % end

%%
% elapsed_time_shapes = toc(shape_times);

for index_el=1:N_element
    
    Elements_idx = Elements(index_el,:);

    PIDu = Elements_idx(1,2);
    ey = Ey_vec(PIDu);
    rho = RHO_vec(PIDu);
    nu = Ni_vec(PIDu);
    
    TypeID = Elements_idx(1,3);
    
    switch TypeID
        case 0 % u element is a tetra4
             new_j_vec =  Elements_idx(1,4:7);             
        case 3 % u element is a hexa8
            new_j_vec = Elements_idx(1,4:11);
        case 4 % u element is a hexa20
            new_j_vec = Elements_idx(1,4:23);
        case 5 % u element is a hexa27
            new_j_vec = Elements_idx(1,4:30);
        case 6 % u element is a hexa32
            new_j_vec = Elements_idx(1,4:35); 
        case 7 % u element is a hexa64
            new_j_vec = Elements_idx(1,4:67);  
    end
     
    Nodes_el = Nodes(new_j_vec,2:4);
    indexes_DOFs_element =zeros(Ndof_node,numel(new_j_vec));
    indexes_DOFs_element (1,:)=(new_j_vec-1)*3+1;
    indexes_DOFs_element (2,:)=(new_j_vec-1)*3+2;
    indexes_DOFs_element (3,:)=(new_j_vec-1)*3+3;
    indexes_DOFs_element = indexes_DOFs_element(:);
    
    [new_jj_vec,new_ii_vec]=meshgrid(indexes_DOFs_element,indexes_DOFs_element);

    
 %% K,M matrices 
     switch TypeID
%       case 0 % u element is a tetra4                  
%              [K, Vel] = stiffnessMatrix_tetra(Nodes_el',ey,nu);
%              [M] = MassMatrix_tetra(Nodes_el',rho);          
        case {3,4,5,6,7} % u element is a hexa8
              % [k_mat,m_mat] = stiffnessMatrix_HEXA (Nodes_el,ey,nu,rho, hexa_type, hexa_order);    
              [k_mat,m_mat] = stiffnessMatrix_HEXA_new (Nodes_el,ey,nu,rho, X,W,N,dShape);  



              i_mat( index_el,:)=new_ii_vec(:)';
              j_mat( index_el,:)=new_jj_vec(:)';
              v_matK(index_el,:)=k_mat(:)';



              v_matM(index_el,:)=m_mat(:)';
%             
     end
  
end

% 
% k_mat = sparse(i_mat,j_mat,v_matK,Ndof_node*N_node,Ndof_node*N_node);
% k_mat = (k_mat+k_mat')*0.5;
% figure,spy(k_mat)
% % elapsed_global_ass_matrices = toc(shape_times),
% largest_K=max(abs(v_matK),[],'all');
% k_filter_zeros = abs(v_matK)<eps(largest_K);
% v_matK(k_filter_zeros)=0;

k_mat = sparse(i_mat,j_mat,v_matK,Ndof_node*N_node,Ndof_node*N_node);
k_mat = (k_mat+k_mat')*0.5;
clear v_matK
% figure,spy(k_mat)
% k_mat = sparse(i_mat(k_filter_zeros),j_mat(k_filter_zeros),v_matK(k_filter_zeros),Ndof_node*N_node,Ndof_node*N_node);

m_mat = sparse(i_mat,j_mat,v_matM,Ndof_node*N_node,Ndof_node*N_node);
m_mat = (m_mat+m_mat')*0.5;
clear v_matM

% elapsed_global_sparse_matrices = toc(shape_times),

% % % k_mat = (k_mat+k_mat')*0.5;
% % % m_mat = (m_mat+m_mat')*0.5;

if(~isConsistent)
    m_mat = diag(sum(m_mat,2));
end

%% remove constrained Dofs
logicalIndConstraints = find(~Constraints(:,1:3)');
m_mat = m_mat(logicalIndConstraints,logicalIndConstraints);
k_mat = k_mat(logicalIndConstraints,logicalIndConstraints);

 % elapsed_global_constrained_matrices = toc(shape_times),
