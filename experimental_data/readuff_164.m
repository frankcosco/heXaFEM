function [lenghtGain,forceGain,tempGain,tempOffs]=...
    readuff_164(dataBlock)
% dataBlock is a cell array of char lines
% nLines = size(dataBlock,1);

tline=dataBlock{2,1};
lenghtGain = getScalar(tline(1:25));
forceGain  = getScalar(tline(26:50));
tempGain   = getScalar(tline(51:75));
tline=dataBlock{3,1};
tempOffs   = getScalar(tline(1:25));

end
function s=getScalar(tfield)
tfield = strrep(tfield,'d-','e-');
s=sscanf(tfield,'%f');
end