function [V,Freq_m]=readuff_55(dataBlock)
% dataBlock is a cell array of char lines
nLines     = size(dataBlock,1);
nLinesData = nLines-8;
Shapes = zeros(4,nLinesData/2);

% UFF.d1 = dataBlock{1,1};
% UFF.d2 = dataBlock{2,1};
% UFF.date = dataBlock{3,1};
% UFF.label = dataBlock{4,1};

dataFrDamp = dataBlock{8,1};
s = sscanf(dataFrDamp(1:13),'%g',1);
w = sscanf(dataFrDamp(14:26),'%g',1);

Freq_m=w/2/pi;

for i=1:2:nLinesData
    %read id
    tline_id = dataBlock{8+i,1};
    tline    = dataBlock{8+i+1,1};
    
    id_i = getUint(tline_id(1:10));
    mx_i = getScalar(tline(14:26));
    my_i = getScalar(tline(40:52));
    mz_i = getScalar(tline(66:78));
    Shapes(:,(i+1)/2)=[id_i mx_i my_i mz_i]';
end

[~,inds]=sort(Shapes(1,:));
Shapes=Shapes(2:4,inds);
V = Shapes(:);

end


function s=getUint(tfield)
% tfield = strrep(tfield,'d-','e-');
s=sscanf(tfield,'%i');
end
function s=getScalar(tfield)
% tfield = strrep(tfield,'d-','e-');
s=sscanf(tfield,'%g');
end