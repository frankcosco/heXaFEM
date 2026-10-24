function [X,W] = HEXA_integration_8()
vc =[-1/sqrt(3), 1/sqrt(3)];
vw =[1,1];

[X1,X2,X3]= meshgrid(vc,vc,vc);
[W1,W2,W3]= meshgrid(vw,vw,vw);

X = [X1(:) X2(:) X3(:)];
W = W1(:).*W2(:).*W3(:);
end
