function [N,dShape, Shape] = HEXA_shapes_8(xi1,xi2,xi3)
% File Shape funct.P38 (File 211517491 Three Dimensional Finite Element Optimization 
%Using the Partial p-Adaptive Method for Stress Analysis of Underground Excavations with Prismatic Cross-sections, P.24) 
% (File CMESUnifiedIsoparametric3DLagrangeFiniteElements - Unified Isoparametric 3D LagrangeFinite Elements, P.3)
% ( File fem_zienkiewicz PP.185-186)

%% Shape function 
Shape= (1/8)*[(1-xi1)*(1-xi2)*(1-xi3);
              (1+xi1)*(1-xi2)*(1-xi3);
              (1+xi1)*(1+xi2)*(1-xi3);
              (1-xi1)*(1+xi2)*(1-xi3);
              (1-xi1)*(1-xi2)*(1+xi3);
              (1+xi1)*(1-xi2)*(1+xi3);
              (1+xi1)*(1+xi2)*(1+xi3);
              (1-xi1)*(1+xi2)*(1+xi3)];

%Shape=[n1..n8] 
%             N=[n1   0    0  n2    0    0  n3   0    0  n4    0    0  n5    0    0  n6   0    0  n7    0  0  n8 0  0;
%                 0  n1   0     0  n2    0   0  n3    0    0  n4    0    0  n5    0    0 n6    0    0  n7 0  0  n8 0;
%                 0    0  n1    0    0  n2   0    0  n3    0    0  n4    0    0  n5    0   0  n6    0    0  n7 0  0  n8;];
%% Compute shape functions derivatives
 N = zeros(3,24);
 N(1,1:3:22)=Shape';
 N(2,2:3:23)=Shape';
 N(3,3:3:24)=Shape';
 
 dShape = (1/8)*[-(1-xi2)*(1-xi3), (1-xi2)*(1-xi3),(1+xi2)*(1-xi3),-(1+xi2)*(1-xi3),-(1-xi2)*(1+xi3),(1-xi2)*(1+xi3),(1+xi2)*(1+xi3),-(1+xi2)*(1+xi3);
                            -(1-xi1)*(1-xi3),-(1+xi1)*(1-xi3),(1+xi1)*(1-xi3),(1-xi1)*(1-xi3),-(1-xi1)*(1+xi3),-(1+xi1)*(1+xi3),(1+xi1)*(1+xi3),(1-xi1)*(1+xi3);
                            -(1-xi1)*(1-xi2),-(1+xi1)*(1-xi2),-(1+xi1)*(1+xi2),-(1-xi1)*(1+xi2),(1-xi1)*(1-xi2),(1+xi1)*(1-xi2),(1+xi1)*(1+xi2),(1-xi1)*(1+xi2)];