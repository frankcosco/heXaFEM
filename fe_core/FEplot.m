function h = FEplot(Elements,Nodes,Constraints)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% nodePlotEnabled=true;


NodesIds = Elements(:,4:11);


     el6=NodesIds(:,[5 8 7  6 ]);% Elements(:,[8 9 10 11]);
     el1=NodesIds(:,[1 4 3  2 ]);
     el2=NodesIds(:,[1 5 6  2 ] );
     el3=NodesIds(:,[2 6 7  3 ]);
     el4=NodesIds(:,[3 7 8  4 ]);
     el5=NodesIds(:,[4 8 5  1 ]);


    % elCrown = [el1;el2;el3;el4;el5;el6;];
   elCrown = [el1;el2;el3;el4;el5;el6;];
    % elpatch = [elpatch; elpatch(:,[1 3 2 4])];
   h = patch('Faces',elCrown,'Vertices',Nodes(:,2:4), ...
              'FaceColor','g','DisplayName','Elements','FaceAlpha',0.9);

    if nargin>2
   hold on
   plot3(Nodes(Constraints(:,1),2),Nodes(Constraints(:,1),3),Nodes(Constraints(:,1),4),'r', ...
       Marker='square',Color='k',MarkerSize=6,MarkerFaceColor='r', LineStyle='none',   ...
       DisplayName='Constraint');
    end

    % if (nodePlotEnabled)
    % plot3(Nodes(:,2),Nodes(:,3),Nodes(:,4),'.b');
    % % text(Nodes(:,2),Nodes(:,3),Nodes(:,4),num2str((1:size(Nodes,1))'));
    % text(Nodes(:,2),Nodes(:,3),Nodes(:,4),num2str(Elements(4:11)'));
    % end
end