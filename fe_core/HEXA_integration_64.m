function [X,W] = HEXA_integration_64()
   
   vw = 0.5 + [-sqrt(30) +sqrt(30) +sqrt(30) -sqrt(30)]/36;
   
   c23 = sqrt(3/7 - 2/7 *sqrt(6/5));
   c14 = sqrt(3/7 + 2/7 *sqrt(6/5));
   vc = [-c14 -c23 c23 c14];
   
   [X1,X2,X3]= meshgrid(vc,vc,vc);
   [W1,W2,W3]= meshgrid(vw,vw,vw);

    X = [X1(:) X2(:) X3(:)];
    W = W1(:).*W2(:).*W3(:);
   
%(File CMESUnifiedIsoparametric3DLagrangeFiniteElements - Unified Isoparametric 3D LagrangeFinite Elements, P.18)
end
