function [Constraints,Elements,Nodes] = brick_hexa8_mesher_advanced( ...
    PID,...
    a,b,c, na,nb,nc,...
    a0,b0,c0)
%% create a mesh for a rectangular parallelepiped
% a,b,c    - dimensions
% na,nb,nc - number of elements along eac direction
% Nodes [GID x y z] by NN
% Elements - Ne x (3 + 8) matrix
%   [ EID PID TypeID G1 G2 G3 G4 [g5 g6 g7 g8 .... gN ]]
%   EID - element ID
%   PID - material property ID
%   TypeID = {0- TETRA4,  1- TETRA10, 2-TETRA20, 3-HEXA8, 4-HEXA(QUAD), 5-HEXA(CUBIC), 6-.... }

% EXA_type = 8;


% dimensionality
nna = na+1;%% n nodes along a
nnb = nb+1;%% n nodes along b
nnc = nc+1;%% n nodes along c
nnab = nna*nnb;

%Element dimensions
ha = a/na;
hb = b/nb;
hc = c/nc;

NN = nna * nnb * nnc;%number of nodes
NE =  na *  nb *  nc;%number of elements
Nodes = zeros(NN,4);
Elements = zeros(NE,3+8);
%         Constraints = false(NN,3);
% build a line of nodes along direction a,b,c
for ia = 0:nna-1
    for ib = 0:nnb-1
        for ic = 0:nnc-1
            % ID del nodo sulla mesh
            gid = ia + nna * ib + nnab * ic + 1 ;
            Nodes(gid,:) = [gid a0+ia*ha b0+ib*hb c0+ic*hc];

            if ia<na && ib<nb && ic<nc
                % ID dell'elemento sulla mesh
                eid = ia + na * (ib + nb * ic) +1;
                Elements(eid,:) = [eid PID 3 gid      gid+1      gid+1+nna      gid+nna  gid+nnab gid+1+nnab gid+1+nna+nnab gid+nna+nnab];
            end
            %                     Constraints(gid,:) = (ia==0);
            %                     Constraints(gid,:) = (ic==0);
        end
    end
end




%Constrains
Constraints = false(NN,3);
%Choose where constrain
indClamp = Nodes(:,3)==(b+b0);
Constraints(indClamp,:)=true;

end