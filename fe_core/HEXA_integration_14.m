function [X,W] = HEXA_integration_14()
 cx1 = 0.75878;
 cx2 = 0.75868;
 cx = 0.79582;
 wx1=0.33516;
 wx=0.88645;
 
    X=[cx2  cx1  cx2;
      -cx2  cx1  cx2;
       cx2 -cx1  cx2;
       cx2  cx1 -cx2;
      -cx2 -cx1 -cx2;
      -cx2  cx1 -cx2;
       cx2 -cx1  cx2;
      -cx2 -cx1  cx2;
       cx   0    0;
      -cx   0    0;
       0    cx   0;
       0   -cx   0;
       0    0    cx;
       0    0   -cx];
    W=[wx1;wx1;wx1;wx1;wx1;wx1;wx1;wx1;
       wx;wx;wx;wx;wx;wx];
end
