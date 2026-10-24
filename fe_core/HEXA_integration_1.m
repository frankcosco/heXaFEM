function [X,W] = HEXA_integration_1()
vc =0;
vw =8;

[X1,X2,X3]= meshgrid(vc,vc,vc);
[W1,W2,W3]= meshgrid(vw,vw,vw);

X = [X1(:) X2(:) X3(:)];
W = W1(:).*W2(:).*W3(:);
end
