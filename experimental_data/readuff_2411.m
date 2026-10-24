function Nodes=readuff_2411(dataBlock)
% dataBlock is a cell array of char lines
nLines = size(dataBlock,1);
Nodes = zeros(4,nLines/2);

for i=1:2:nLines
    tline=dataBlock{i,1};
    %read id
    id_i = getUint(tline(1:10));
    tline=dataBlock{i+1,1};
    nx_i = getScalar(tline(1:25));
    ny_i = getScalar(tline(26:50));
    nz_i = getScalar(tline(51:75));
    Nodes(:,(i+1)/2)=[id_i nx_i ny_i nz_i]';
end

end


function s=getUint(tfield)
% tfield = strrep(tfield,'d-','e-');
s=sscanf(tfield,'%u');
end
function s=getScalar(tfield)
% tfield = strrep(tfield,'d-','e-');
s=sscanf(tfield,'%f');
end