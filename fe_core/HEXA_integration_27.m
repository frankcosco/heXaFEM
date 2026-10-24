function [X,W] = HEXA_integration_27()
vc =[-sqrt(3/5), 0,sqrt(3/5)];
vw =[5/9, 8/9, 5/9];
[X1,X2,X3]= meshgrid(vc,vc,vc);
[W1,W2,W3]= meshgrid(vw,vw,vw);

X = [X1(:) X2(:) X3(:)];
W = W1(:).*W2(:).*W3(:);
end
