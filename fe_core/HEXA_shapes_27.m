function [N,dShape, Shape] = HEXA_shapes_27(xi1,xi2,xi3)

 % File Shape funct.P40 (File Browning-2020-Asecond-order19-nodepyramidfiniteelementsuit
% A SECOND-ORDER 19-NODE PYRAMID FINITE ELEMENT SUITABLE FOR
% LUMPED MASS EXPLICIT DYNAMIC METHODS IN NONLINEAR SOLID MECHANICS, P.277) 
     %(File 211517491 Three Dimensional Finite Element Optimization 
%Using the Partial p-Adaptive Method for Stress Analysis of Underground Excavations with Prismatic Cross-sections, P.24)       
% (File CMESUnifiedIsoparametric3DLagrangeFiniteElements - Unified Isoparametric 3D LagrangeFinite Elements, P.4)   
%% Shape function 
Shape      =[... 
     -1/8*xi1*xi2*xi3*(1-xi1)*(1-xi2)*(1-xi3),
      1/8*xi1*xi2*xi3*(1+xi1)*(1-xi2)*(1-xi3), 
     -1/8*xi1*xi2*xi3*(1+xi1)*(1+xi2)*(1-xi3),
      1/8*xi1*xi2*xi3*(1-xi1)*(1+xi2)*(1-xi3), 
      1/8*xi1*xi2*xi3*(1-xi1)*(1-xi2)*(1+xi3),
     -1/8*xi1*xi2*xi3*(1+xi1)*(1-xi2)*(1+xi3),   
      1/8*xi1*xi2*xi3*(1+xi1)*(1+xi2)*(1+xi3),
     -1/8*xi1*xi2*xi3*(1-xi1)*(1+xi2)*(1+xi3),     
      1/4*xi2*xi3*(1-xi1^2)*(1-xi2)*(1-xi3),%9
     -1/4*xi1*xi3*(1+xi1)*(1-xi2^2)*(1-xi3),
     -1/4*xi2*xi3*(1-xi1^2)*(1+xi2)*(1-xi3),
      1/4*xi1*xi3*(1-xi1)*(1-xi2^2)*(1-xi3),
     -1/4*xi2*xi3*(1-xi1^2)*(1-xi2)*(1+xi3),
      1/4*xi1*xi3*(1+xi1)*(1-xi2^2)*(1+xi3),
      1/4*xi2*xi3*(1-xi1^2)*(1+xi2)*(1+xi3),
     -1/4*xi1*xi3*(1-xi1)*(1-xi2^2)*(1+xi3),
      1/4*xi1*xi2*(1-xi1)*(1-xi2)*(1-xi3^2),%17
     -1/4*xi1*xi2*(1+xi1)*(1-xi2)*(1-xi3^2),
      1/4*xi1*xi2*(1+xi1)*(1+xi2)*(1-xi3^2),
     -1/4*xi1*xi2*(1-xi1)*(1+xi2)*(1-xi3^2),     
     -1/2*xi3*(1-xi1^2)*(1-xi2^2)*(1-xi3),%21
      1/2*xi3*(1-xi1^2)*(1-xi2^2)*(1+xi3),
     -1/2*xi2*(1-xi1^2)*(1-xi2)*(1-xi3^2),
      1/2*xi1*(1+xi1)*(1-xi2^2)*(1-xi3^2),
      1/2*xi2*(1-xi1^2)*(1+xi2)*(1-xi3^2),
     -1/2*xi1*(1-xi1)*(1-xi2^2)*(1-xi3^2),     
     (1-xi1^2)*(1-xi2^2)*(1-xi3^2)]; %#ok<*COMNL>

%% Compute shape functions derivatives
 dShape = zeros(3,27);
 N = zeros(3,81); 
 
 dShape(1,:) =[...
 (xi2*xi3*(2*xi1 - 1)*(xi2 - 1)*(xi3 - 1))/8,
 (xi2*xi3*(2*xi1 + 1)*(xi2 - 1)*(xi3 - 1))/8,
 (xi2*xi3*(2*xi1 + 1)*(xi2 + 1)*(xi3 - 1))/8,
 (xi2*xi3*(2*xi1 - 1)*(xi2 + 1)*(xi3 - 1))/8,
 (xi2*xi3*(2*xi1 - 1)*(xi2 - 1)*(xi3 + 1))/8,
 (xi2*xi3*(2*xi1 + 1)*(xi2 - 1)*(xi3 + 1))/8,
 (xi2*xi3*(2*xi1 + 1)*(xi2 + 1)*(xi3 + 1))/8,
 (xi2*xi3*(2*xi1 - 1)*(xi2 + 1)*(xi3 + 1))/8,
-(xi1*xi2*xi3*(xi2 - 1)*(xi3 - 1))/2,
-(xi3*(2*xi1 + 1)*(xi2^2 - 1)*(xi3 - 1))/4,
-(xi1*xi2*xi3*(xi2 + 1)*(xi3 - 1))/2,
-(xi3*(2*xi1 - 1)*(xi2^2 - 1)*(xi3 - 1))/4,
-(xi1*xi2*xi3*(xi2 - 1)*(xi3 + 1))/2,
-(xi3*(2*xi1 + 1)*(xi2^2 - 1)*(xi3 + 1))/4,
-(xi1*xi2*xi3*(xi2 + 1)*(xi3 + 1))/2,
-(xi3*(2*xi1 - 1)*(xi2^2 - 1)*(xi3 + 1))/4,
-(xi2*(2*xi1 - 1)*(xi3^2 - 1)*(xi2 - 1))/4,
-(xi2*(2*xi1 + 1)*(xi3^2 - 1)*(xi2 - 1))/4,
-(xi2*(2*xi1 + 1)*(xi3^2 - 1)*(xi2 + 1))/4,
-(xi2*(2*xi1 - 1)*(xi3^2 - 1)*(xi2 + 1))/4,
  xi1*xi3*(xi2^2 - 1)*(xi3 - 1),
  xi1*xi3*(xi2^2 - 1)*(xi3 + 1),
  xi1*xi2*(xi3^2 - 1)*(xi2 - 1),
((2*xi1 + 1)*(xi2^2 - 1)*(xi3^2 - 1))/2,
  xi1*xi2*(xi3^2 - 1)*(xi2 + 1),
((2*xi1 - 1)*(xi2^2 - 1)*(xi3^2 - 1))/2,
 -2*xi1*(xi2^2 - 1)*(xi3^2 - 1)];

dShape(2,:) =[...
    (xi1*xi3*(2*xi2 - 1)*(xi1 - 1)*(xi3 - 1))/8,
    (xi1*xi3*(2*xi2 - 1)*(xi1 + 1)*(xi3 - 1))/8,
    (xi1*xi3*(2*xi2 + 1)*(xi1 + 1)*(xi3 - 1))/8,
    (xi1*xi3*(2*xi2 + 1)*(xi1 - 1)*(xi3 - 1))/8,
    (xi1*xi3*(2*xi2 - 1)*(xi1 - 1)*(xi3 + 1))/8,
    (xi1*xi3*(2*xi2 - 1)*(xi1 + 1)*(xi3 + 1))/8,
    (xi1*xi3*(2*xi2 + 1)*(xi1 + 1)*(xi3 + 1))/8,
    (xi1*xi3*(2*xi2 + 1)*(xi1 - 1)*(xi3 + 1))/8,
   -(xi3*(2*xi2 - 1)*(xi1^2 - 1)*(xi3 - 1))/4,
   -(xi1*xi2*xi3*(xi1 + 1)*(xi3 - 1))/2,
   -(xi3*(2*xi2 + 1)*(xi1^2 - 1)*(xi3 - 1))/4,
   -(xi1*xi2*xi3*(xi1 - 1)*(xi3 - 1))/2,
   -(xi3*(2*xi2 - 1)*(xi1^2 - 1)*(xi3 + 1))/4,
   -(xi1*xi2*xi3*(xi1 + 1)*(xi3 + 1))/2,
   -(xi3*(2*xi2 + 1)*(xi1^2 - 1)*(xi3 + 1))/4,
   -(xi1*xi2*xi3*(xi1 - 1)*(xi3 + 1))/2,
   -(xi1*(2*xi2 - 1)*(xi3^2 - 1)*(xi1 - 1))/4,
   -(xi1*(2*xi2 - 1)*(xi3^2 - 1)*(xi1 + 1))/4,
   -(xi1*(2*xi2 + 1)*(xi3^2 - 1)*(xi1 + 1))/4,
   -(xi1*(2*xi2 + 1)*(xi3^2 - 1)*(xi1 - 1))/4,
     xi2*xi3*(xi1^2 - 1)*(xi3 - 1),
     xi2*xi3*(xi1^2 - 1)*(xi3 + 1),
   ((2*xi2 - 1)*(xi1^2 - 1)*(xi3^2 - 1))/2,
     xi1*xi2*(xi3^2 - 1)*(xi1 + 1),
   ((2*xi2 + 1)*(xi1^2 - 1)*(xi3^2 - 1))/2,
     xi1*xi2*(xi3^2 - 1)*(xi1 - 1),
    -2*xi2*(xi1^2 - 1)*(xi3^2 - 1)];

dShape(3,:) =[...
    (xi1*xi2*(2*xi3 - 1)*(xi1 - 1)*(xi2 - 1))/8,
    (xi1*xi2*(2*xi3 - 1)*(xi1 + 1)*(xi2 - 1))/8,
    (xi1*xi2*(2*xi3 - 1)*(xi1 + 1)*(xi2 + 1))/8,
    (xi1*xi2*(2*xi3 - 1)*(xi1 - 1)*(xi2 + 1))/8,
    (xi1*xi2*(2*xi3 + 1)*(xi1 - 1)*(xi2 - 1))/8,
    (xi1*xi2*(2*xi3 + 1)*(xi1 + 1)*(xi2 - 1))/8,
    (xi1*xi2*(2*xi3 + 1)*(xi1 + 1)*(xi2 + 1))/8,
    (xi1*xi2*(2*xi3 + 1)*(xi1 - 1)*(xi2 + 1))/8,
   -(xi2*(2*xi3 - 1)*(xi1^2 - 1)*(xi2 - 1))/4,
   -(xi1*(2*xi3 - 1)*(xi2^2 - 1)*(xi1 + 1))/4,
   -(xi2*(2*xi3 - 1)*(xi1^2 - 1)*(xi2 + 1))/4,
   -(xi1*(2*xi3 - 1)*(xi2^2 - 1)*(xi1 - 1))/4,
   -(xi2*(2*xi3 + 1)*(xi1^2 - 1)*(xi2 - 1))/4,
   -(xi1*(2*xi3 + 1)*(xi2^2 - 1)*(xi1 + 1))/4,
   -(xi2*(2*xi3 + 1)*(xi1^2 - 1)*(xi2 + 1))/4,
   -(xi1*(2*xi3 + 1)*(xi2^2 - 1)*(xi1 - 1))/4,
   -(xi1*xi2*xi3*(xi1 - 1)*(xi2 - 1))/2,
   -(xi1*xi2*xi3*(xi1 + 1)*(xi2 - 1))/2,
   -(xi1*xi2*xi3*(xi1 + 1)*(xi2 + 1))/2,
   -(xi1*xi2*xi3*(xi1 - 1)*(xi2 + 1))/2,
   ((2*xi3 - 1)*(xi1^2 - 1)*(xi2^2 - 1))/2,
   ((2*xi3 + 1)*(xi1^2 - 1)*(xi2^2 - 1))/2,
     xi2*xi3*(xi1^2 - 1)*(xi2 - 1),
     xi1*xi3*(xi2^2 - 1)*(xi1 + 1),
     xi2*xi3*(xi1^2 - 1)*(xi2 + 1),
     xi1*xi3*(xi2^2 - 1)*(xi1 - 1),
    -2*xi3*(xi1^2 - 1)*(xi2^2 - 1)];

for i=1:27
N(1,(i-1)*3+1)=Shape(i);
N(2,(i-1)*3+2)=Shape(i);
N(3,(i-1)*3+3)=Shape(i);
end       
   
end
