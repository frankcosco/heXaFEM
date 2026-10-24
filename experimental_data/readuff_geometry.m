function Nodes = readuff_geometry(file_input)
%% readuff_geometry 
if nargin<1 %for debug purpose only
    file_input='geometry4.unv';
end
% dofXnode = 3;

[dataBlocks] = readuff_generic(file_input);
%%
typeBlocksArray=cell2mat( dataBlocks(:,2));
indUnitBlock = typeBlocksArray==164;


[lenghtGain,forceGain,tempGain,tempOffs]=readuff_164(dataBlocks{indUnitBlock,1});
% end

indGeometryBlock = typeBlocksArray==2411;

Nodes=readuff_2411(dataBlocks{indGeometryBlock,1});

figure(1),
cla,
scatter3(Nodes(2,:),Nodes(3,:),Nodes(4,:));
% axis equal
axis auto