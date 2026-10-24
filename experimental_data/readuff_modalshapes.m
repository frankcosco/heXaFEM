function [Mode_shapes,Modal_freqs,Nodes] = readuff_modalshapes(file_input)
%% readuff_geometry 
if nargin<1 %for debug purpose only
    file_input='real_modeset.unv';
end
% dofXnode = 3;

[dataBlocks] = readuff_generic(file_input);
%%
typeBlocksArray=cell2mat( dataBlocks(:,2));
indUnitBlock = typeBlocksArray==164;


[lenghtGain,forceGain,tempGain,tempOffs]=readuff_164(dataBlocks{indUnitBlock,1});
% end

Nodes=readuff_2411(dataBlocks{typeBlocksArray==2411,1});
Nn = size(Nodes,2);

indModesBlock = find(typeBlocksArray==55);
Nm = numel(indModesBlock);
Mode_shapes = zeros(Nn*3,Nm);
Modal_freqs = zeros(1,Nm);


m_counter=1;
for i=1:Nm
    im=indModesBlock(i);    
    [Mode_shapes(:,m_counter),...
       Modal_freqs(:,m_counter)] =readuff_55(dataBlocks{im,1});
    m_counter = m_counter+1;


end


Nodes  = Nodes';

