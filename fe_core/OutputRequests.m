function [indexOutputDofs,indexOutputNodes]=OutputRequests(Nodes,logicalIndConstrained,grid_c)
Nc = size(grid_c,1); % number of nodes for which we request output
Nn = size(Nodes,1);  % number of nodes  in mesh

logicalIndControls = false(size(logicalIndConstrained));
indexOutputNodes = zeros(1,Nc);
%find outputNodes into Nodes
for i=1:Nc
    err = sqrt(sum( (repmat(grid_c(i,2:4),Nn,1)-Nodes(:,2:4)).^2,2));
    [~,ind_c_i] = min(err);
    %    ind_c_i =ind_c(1);

    if logicalIndConstrained(1,ind_c_i)
        logicalIndControls(:,ind_c_i)=true;
        indexOutputNodes(i)=ind_c_i;
    end

    % figure(1),cla,
    % plot3(Nodes(:,2),Nodes(:,3),Nodes(:,4),'.k');hold on
    % plot3(grid_c(:,2),grid_c(:,3),grid_c(:,4),'.b');
    %
    % plot3(Nodes(ind_c_i,2),Nodes(ind_c_i,3),Nodes(ind_c_i,4),'or');
    % plot3(grid_c(i,2),grid_c(i,3),grid_c(i,4),'om');
    % axis equal
end
indexOutputDofs=find(logicalIndControls(logicalIndConstrained));
