function [Constraints,Elements,Nodes] = brick_hexa_mesher( ...
    EXA_type, PID,...
    a,b,c, na,nb,nc)
%% create a mesh for a rectangular parallelepiped
% a,b,c    - dimensions
% na,nb,nc - number of elements along eac direction
% Nodes [GID x y z] by NN
% Elements - Ne x (3 + 8) matrix
%   [ EID PID TypeID G1 G2 G3 G4 [g5 g6 g7 g8 .... gN ]]
%   EID - element ID
%   PID - material property ID
%   TypeID = {0- TETRA4,  1- TETRA10, 2-TETRA20, 3-HEXA8, 4-HEXA(QUAD), 5-HEXA(CUBIC), 6-.... }

switch(EXA_type)
     case 8  % trilinear shapes element
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
                    Nodes(gid,:) = [gid ia*ha ib*hb ic*hc];
                    
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
                
     case 20 % quadratic (element)
        nna = 2*na+1;%% n nodes along a
        nnb = 2*nb+1;%% n nodes along b
        nnc = 2*nc+1;%% n nodes along c
        nnab = nna*nnb;
        
        ha = a/na;
        hb = b/nb;
        hc = c/nc;
                
        NN = nna * nnb * nnc;%number of nodes
        NE =  na *  nb *  nc;%number of elements
        Nodes = zeros(NN,4);
        Elements = zeros(NE,3+20);
        idToremove = zeros(NE,7);
       
        for ia = 0:0.5:na
            for ib = 0:0.5:nb
                for ic = 0:0.5:nc
                    % ID del nodo sulla mesh
                    gid =2*( ia + nna * ib + nnab * ic)+ 1 ;
                    Nodes(gid,:) = [gid ia*ha ib*hb ic*hc];
                    
                    if ia<na && ib<nb && ic<nc && mod(ia*2,2)==0 && mod(ib*2,2)==0 && mod(ic*2,2)==0
                        % ID dell'elemento sulla mesh
                        eid = ia + na * (ib + nb * ic) +1;
                        Elements(eid,:) = [eid PID 4 ...
                            gid            gid+2             gid+2+2*nna            gid+2*nna  ...
                            gid+2*nnab     gid+2+2*nnab      gid+2+2*nna+2*nnab     gid+2*nna+2*nnab ...
                            gid+1          gid+2+nna         gid+1+2*nna            gid+nna ...
                            gid+1+2*nnab   gid+2+nna+2*nnab  gid+1+2*nna+2*nnab     gid+nna+2*nnab...
                            gid+nnab       gid+2+nnab        gid+2+2*nna+nnab       gid+2*nna+nnab  ...
                            ];%27
                        
       idToremove(eid,:) = [gid+1+nna      gid+1+nna+2*nnab ... 21 22
                            gid+1+nnab     gid+2+nna+nnab    gid+1+2*nna+nnab       gid+nna+nnab ... 23 24 25 26
                            gid+1+nna+nnab];
                                            
                    end
                    %                     Constraints(gid,:) = (ia==0);
                    %                     Constraints(gid,:) = (ic==0);
                end
            end
        end
%Code to cleanup not used nodes
        id2rem = unique(idToremove(:));
        N2rem =numel(id2rem);
        Nodes(id2rem,1) = NN+1;
        
        [~,indSorted, indSortedBack] = unique( Nodes(:,1));
        indSortedUtile = indSorted(1:NN-N2rem);
        Nodes = Nodes( indSortedUtile, :);
        ElementsNew=[Elements(:,1:3) indSortedBack(Elements(:,4:end))];
        Elements = ElementsNew;
        NN=size(Nodes,1);
        
     case 27 % tri quadratic element
        nna = 2*na+1;%% n nodes along a
        nnb = 2*nb+1;%% n nodes along b
        nnc = 2*nc+1;%% n nodes along c
        nnab = nna*nnb;
        
        ha = a/na;%%element dimensions
        hb = b/nb;
        hc = c/nc;
                
        NN = nna * nnb * nnc;%number of nodes
        NE =  na *  nb *  nc;%number of elements
        Nodes = zeros(NN,4);
        Elements = zeros(NE,3+27);
        %         Constraints = false(NN,3);
        for ia = 0:0.5:na
            for ib = 0:0.5:nb
                for ic = 0:0.5:nc
                    % ID del nodo sulla mesh
                    gid =2*( ia + nna * ib + nnab * ic)+ 1 ;
                    Nodes(gid,:) = [gid ia*ha ib*hb ic*hc];
                    
                    if ia<na && ib<nb && ic<nc && mod(ia*2,2)==0 && mod(ib*2,2)==0 && mod(ic*2,2)==0
                        % ID dell'elemento sulla mesh
                        eid = ia + na * (ib + nb * ic) +1;
                        Elements(eid,:) = [eid PID 5 ...
                            gid             gid+2               gid+2+2*nna             gid+2*nna  ...
                            gid+2*nnab      gid+2+2*nnab        gid+2+2*nna+2*nnab      gid+2*nna+2*nnab ...
                            gid+1           gid+2+nna           gid+1+2*nna             gid+nna ...
                            gid+1+2*nnab    gid+2+nna+2*nnab    gid+1+2*nna+2*nnab      gid+nna+2*nnab...
                            gid+nnab        gid+2+nnab          gid+2+2*nna+nnab        gid+2*nna+nnab  ...
                            gid+1+nna       gid+1+nna+2*nnab ... 21 22
                            gid+1+nnab      gid+2+nna+nnab      gid+1+2*nna+nnab        gid+nna+nnab ... 23 24 25 26
                            gid+1+nna+nnab];%27
                       
                    end
                    %                     Constraints(gid,:) = (ia==0);
                    %                     Constraints(gid,:) = (ic==0);
                end
            end
        end
        
     case 32%  tri cubic element
        nna = 3*na+1;%% n nodes along a
        nnb = 3*nb+1;%% n nodes along b
        nnc = 3*nc+1;%% n nodes along c
        nnab = nna*nnb;
        
        ha = a/na;%%element dimensions
        hb = b/nb;
        hc = c/nc;
        
        NN = nna * nnb * nnc;%number of nodes
        NE =  na *  nb *  nc;%number of elements
        Nodes = zeros(NN,4);
        Elements = zeros(NE,3+32);
        idToremove = zeros(NE,32);
        %         Constraints = false(NN,3);
       
        for ic = 0:1/3:nc
            for ib = 0:1/3:nb
                for ia = 0:1/3:na
                    
                    % ID del nodo sulla mesh
                    gid =round(3*( ia + nna * ib + nnab * ic)+ 1) ;
                    %                     disp(gid);
                    Nodes(gid,:) = [gid ia*ha ib*hb ic*hc];
                    
                    if ia<na && ib<nb && ic<nc && mod(ia*3,3)==0 && mod(ib*3,3)==0 && mod(ic*3,3)==0
                        % ID dell'elemento sulla mesh
                        %                         disp(gid)
                        eid = ia + na * (ib + nb * ic) +1;
                                       
                        Elements(eid,:) = [eid PID 6 ...
                            gid                 gid+3               gid+3+3*nna        gid+3*nna  ...         % 1 2 3 4
                            gid+3*nnab          gid+3+3*nnab        gid+3+3*nna+3*nnab gid+3*nna+3*nnab ...   % 5 6 7 8
                            gid+nnab            gid+3+nnab          gid+3+3*nna+nnab   gid+3*nna+nnab ...     % 9 10 11 12
                            gid+2*nnab          gid+3+2*nnab        gid+3+3*nna+2*nnab gid+3*nna+2*nnab ...   %13 14 15 16
                            gid+1               gid+2               gid+3+nna          gid+3+2*nna ...        % 17 18 19 20
                            gid+2+3*nna         gid+1+3*nna         gid+2*nna          gid+nna ...            % 21 22 23 24 
                            gid+1+3*nnab        gid+2+3*nnab        gid+3+nna+3*nnab   gid+3+2*nna+3*nnab ... % 25 26 27 28
                            gid+2+3*nna+3*nnab  gid+1+3*nna+3*nnab  gid+2*nna+3*nnab   gid+nna+3*nnab];       % 29 30 31 32
                         
        idToremove(eid,:) =[gid+1+nna           gid+2+nna           gid+2+2*nna        gid+1+2*nna ...        %25 26 27 28 of mesh 64
                            gid+1+nna+3*nnab    gid+2+nna+3*nnab    gid+2+2*nna+3*nnab gid+1+2*nna+3*nnab ... %37 38 39 40 of mesh 64
                            gid+1+nnab          gid+2+nnab          gid+3+nna+nnab     gid+3+2*nna+nnab      gid+2+nnab+3*nna   gid+1+nnab+3*nna ... %41 42 43 44 45 46 of mesh 64
                            gid+2*nna+nnab      gid+nna+nnab        gid+1+nna+nnab     gid+2+nna+nnab        gid+2+2*nna+nnab   gid+1+2*nna+nnab ... %47 48 49 50 51 52 of mesh 64
                            gid+1+2*nnab        gid+2+2*nnab        gid+3+nna+2*nnab   gid+3+2*nna+2*nnab    gid+2+3*nna+2*nnab gid+1+3*nna+2*nnab ... %53 54 55 56 57 58 of mesh 64
                            gid+2*nna+2*nnab    gid+nna+2*nnab      gid+1+nna+2*nnab   gid+2+nna+2*nnab      gid+2+2*nna+2*nnab gid+1+2*nna+2*nnab]; %59 60 61 62 63 64 of mesh 64
%                                     
                            
                            
                        
                    end
                    %                     Constraints(gid,:) = (ia==0);
                    %                     Constraints(gid,:) = (ic==0);
                end
            end
        end             
%Code to cleanup not used nodes
        id2rem = unique(idToremove(:));
        N2rem =numel(id2rem);
        Nodes(id2rem,1) = NN+1;
        
        [~,indSorted, indSortedBack] = unique( Nodes(:,1));
        indSortedUtile = indSorted(1:NN-N2rem);
        Nodes = Nodes( indSortedUtile, :);
        ElementsNew=[Elements(:,1:3) indSortedBack(Elements(:,4:end))];
        Elements = ElementsNew;
        NN=size(Nodes,1);      
              
     case 64%  tri cubic element
        nna = 3*na+1;%% n nodes along a
        nnb = 3*nb+1;%% n nodes along b
        nnc = 3*nc+1;%% n nodes along c
        nnab = nna*nnb;
        
        ha = a/na;%%element dimensions
        hb = b/nb;
        hc = c/nc;
        
        NN = nna * nnb * nnc;%number of nodes
        NE =  na *  nb *  nc;%number of elements
        Nodes = zeros(NN,4);
        Elements = zeros(NE,3+64);
        %         Constraints = false(NN,3);
       
        for ic = 0:1/3:nc
            for ib = 0:1/3:nb
                for ia = 0:1/3:na
                    
                    gid =round(3*( ia + nna * ib + nnab * ic)+ 1) ;
                    %                     disp(gid);
                    Nodes(gid,:) = [gid ia*ha ib*hb ic*hc];
                    
                    if ia<na && ib<nb && ic<nc && mod(ia*3,3)==0 && mod(ib*3,3)==0 && mod(ic*3,3)==0
                        % ID dell'elemento sulla mesh
                        %                         disp(gid)
                        eid = ia + na * (ib + nb * ic) +1;
                        ID_elements_17_to_28 = [...
                           gid+1           gid+2          gid+3+nna          gid+3+2*nna ...       % 17 18 19 20
                           gid+2+3*nna     gid+1+3*nna    gid+2*nna          gid+nna ...           % 24
                           gid+1+nna       gid+2+nna      gid+2+2*nna        gid+1+2*nna];         % 25 26 27 28 
                        
                        Elements(eid,:) = [eid PID 7 ...
                           gid             gid+3          gid+3+3*nna        gid+3*nna  ...       % 1 2 3 4
                           gid+3*nnab      gid+3+3*nnab   gid+3+3*nna+3*nnab gid+3*nna+3*nnab ... % 5 6 7 8
                           gid+nnab        gid+3+nnab     gid+3+3*nna+nnab   gid+3*nna+nnab ...   % 9 10 11 12 (1 2 3 4 + nnab)
                           gid+2*nnab      gid+3+2*nnab   gid+3+3*nna+2*nnab gid+3*nna+2*nnab ... %16
                            ID_elements_17_to_28 ... 
                            ID_elements_17_to_28+3*nnab ... 29 to 40
                            ID_elements_17_to_28+nnab   ... 41 to 52
                            ID_elements_17_to_28+2*nnab]; % 53 to 64
%                                                   
                    end
                    %                     Constraints(gid,:) = (ia==0);
                    %                     Constraints(gid,:) = (ic==0);
                end
            end
        end
        %         indClamp = Nodes(:,2)==0;
        %         Constraints(indClamp,:)=true;
        
end

%Constrains 
Constraints = false(NN,3);
%Choose where constrain
indClamp = (Nodes(:,3)-b).^2<1e-6;
Constraints(indClamp,:)=true;

end