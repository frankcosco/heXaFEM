function [N,dShape,Shape] = HEXA_shapes_20(xi1,xi2,xi3)
 % File Shape funct.P39 (File 211517491 Three Dimensional Finite Element Optimization 
%Using the Partial p-Adaptive Method for Stress Analysis of Underground Excavations with Prismatic Cross-sections, P.26)
% ( File fem_zienkiewicz PP.185-186) 
% (File liu_gr_quek_ss_the_finite_element_method_a_practical_course P-219-221)
% (chari_mvk_salon_sj_numerical_methods_in_electromagnetism PP.251-253)

%% Shape function     
 Shape=(1/4)*[...
     (1/2)*(1-xi1)*(1-xi2)*(1-xi3)*(-xi1-xi2-xi3-2),
     (1/2)*(1+xi1)*(1-xi2)*(1-xi3)*(xi1-xi2-xi3-2),
     (1/2)*(1+xi1)*(1+xi2)*(1-xi3)*(xi1+xi2-xi3-2),
     (1/2)*(1-xi1)*(1+xi2)*(1-xi3)*(-xi1+xi2-xi3-2),
     (1/2)*(1-xi1)*(1-xi2)*(1+xi3)*(-xi1-xi2+xi3-2),
     (1/2)*(1+xi1)*(1-xi2)*(1+xi3)*(xi1-xi2+xi3-2),
     (1/2)*(1+xi1)*(1+xi2)*(1+xi3)*(xi1+xi2+xi3-2),  
     (1/2)*(1-xi1)*(1+xi2)*(1+xi3)*(-xi1+xi2+xi3-2),
     (1-xi1^2)*(1-xi2)*(1-xi3), 
     (1-xi2^2)*(1+xi1)*(1-xi3), 
     (1-xi1^2)*(1+xi2)*(1-xi3), 
     (1-xi2^2)*(1-xi1)*(1-xi3),
     (1-xi1^2)*(1-xi2)*(1+xi3), 
     (1-xi2^2)*(1+xi1)*(1+xi3), 
     (1-xi1^2)*(1+xi2)*(1+xi3),  
     (1-xi2^2)*(1-xi1)*(1+xi3),  
     (1-xi3^2)*(1-xi1)*(1-xi2), 
     (1-xi3^2)*(1+xi1)*(1-xi2),  
     (1-xi3^2)*(1+xi1)*(1+xi2),  
     (1-xi3^2)*(1-xi1)*(1+xi2)];

%% Compute shape functions derivatives
dShape = zeros(3,20);
     N = zeros(3,60);
dShape(1,:) =(1/4)*[(1/2)*((xi2 - 1)*(xi3 - 1)*(2*xi1 + xi2 + xi3 + 1)),  -(1/2)*((xi2 - 1)*(xi3 - 1)*(xi2 - 2*xi1 + xi3 + 1)),  -(1/2)*((xi2 + 1)*(xi3 - 1)*(2*xi1 + xi2 - xi3 - 1)),  -(1/2)*((xi2 + 1)*(xi3 - 1)*(2*xi1 - xi2 + xi3 + 1)),  -(1/2)*((xi2 - 1)*(xi3 + 1)*(2*xi1 + xi2 - xi3 + 1)),  -(1/2)*((xi2 - 1)*(xi3 + 1)*(2*xi1 - xi2 + xi3 - 1)),  (1/2)*((xi2 + 1)*(xi3 + 1)*(2*xi1 + xi2 + xi3 - 1)),  (1/2)*((xi2 + 1)*(xi3 + 1)*(2*xi1 - xi2 - xi3 + 1)),  -2*(xi1*(xi2 - 1)*(xi3 - 1)),  ((xi2^2 - 1)*(xi3 - 1)),       2*(xi1*(xi2 + 1)*(xi3 - 1)),  -((xi2^2 - 1)*(xi3 - 1)),      2*(xi1*(xi2 - 1)*(xi3 + 1)),  -((xi2^2 - 1)*(xi3 + 1)),      -2*(xi1*(xi2 + 1)*(xi3 + 1)),  ((xi2^2 - 1)*(xi3 + 1)),      -((xi3^2 - 1)*(xi2 - 1)),    ((xi3^2 - 1)*(xi2 - 1)),     -((xi3^2 - 1)*(xi2 + 1)),     ((xi3^2 - 1)*(xi2 + 1))]; 
dShape(2,:)= (1/4)*[(1/2)*((xi1 - 1)*(xi3 - 1)*(xi1 + 2*xi2 + xi3 + 1)),  -(1/2)*((xi1 + 1)*(xi3 - 1)*(2*xi2 - xi1 + xi3 + 1)),  -(1/2)*((xi1 + 1)*(xi3 - 1)*(xi1 + 2*xi2 - xi3 - 1)),  -(1/2)*((xi1 - 1)*(xi3 - 1)*(xi1 - 2*xi2 + xi3 + 1)),  -(1/2)*((xi1 - 1)*(xi3 + 1)*(xi1 + 2*xi2 - xi3 + 1)),  -(1/2)*((xi1 + 1)*(xi3 + 1)*(xi1 - 2*xi2 + xi3 - 1)),  (1/2)*((xi1 + 1)*(xi3 + 1)*(xi1 + 2*xi2 + xi3 - 1)),  (1/2)*((xi1 - 1)*(xi3 + 1)*(xi1 - 2*xi2 - xi3 + 1)),  -((xi1^2 - 1)*(xi3 - 1)),       2*(xi2*(xi1 + 1)*(xi3 - 1)),  (xi1^2 - 1)*(xi3 - 1),       -2*(xi2*(xi1 - 1)*(xi3 - 1)),  ((xi1^2 - 1)*(xi3 + 1)),      -2*(xi2*(xi1 + 1)*(xi3 + 1)),  -((xi1^2 - 1)*(xi3 + 1)),       2*(xi2*(xi1 - 1)*(xi3 + 1)), -((xi3^2 - 1)*(xi1 - 1)),    ((xi3^2 - 1)*(xi1 + 1)),     -((xi3^2 - 1)*(xi1 + 1)),     ((xi3^2 - 1)*(xi1 - 1))];            
dShape(3,:)= (1/4)*[(1/2)*((xi1 - 1)*(xi2 - 1)*(xi1 + xi2 + 2*xi3 + 1)),  -(1/2)*((xi1 + 1)*(xi2 - 1)*(xi2 - xi1 + 2*xi3 + 1)),  -(1/2)*((xi1 + 1)*(xi2 + 1)*(xi1 + xi2 - 2*xi3 - 1)),  -(1/2)*((xi1 - 1)*(xi2 + 1)*(xi1 - xi2 + 2*xi3 + 1)),  -(1/2)*((xi1 - 1)*(xi2 - 1)*(xi1 + xi2 - 2*xi3 + 1)),  -(1/2)*((xi1 + 1)*(xi2 - 1)*(xi1 - xi2 + 2*xi3 - 1)),  (1/2)*((xi1 + 1)*(xi2 + 1)*(xi1 + xi2 + 2*xi3 - 1)),  (1/2)*((xi1 - 1)*(xi2 + 1)*(xi1 - xi2 - 2*xi3 + 1)),  -((xi1^2 - 1)*(xi2 - 1)),      ((xi2^2 - 1)*(xi1 + 1)),       (xi1^2 - 1)*(xi2 + 1),       -((xi2^2 - 1)*(xi1 - 1)),      ((xi1^2 - 1)*(xi2 - 1)),      -((xi2^2 - 1)*(xi1 + 1)),      -((xi1^2 - 1)*(xi2 + 1)),      ((xi2^2 - 1)*(xi1 - 1)),     -2*(xi3*(xi1 - 1)*(xi2 - 1)), 2*(xi3*(xi1 + 1)*(xi2 - 1)), -2*(xi3*(xi1 + 1)*(xi2 + 1)), 2*(xi3*(xi1 - 1)*(xi2 + 1))]; 

for i=1:20
N(1,(i-1)*3+1)=Shape(i);
N(2,(i-1)*3+2)=Shape(i);
N(3,(i-1)*3+3)=Shape(i);
end
          
   
end
