function [dataBlocks] = readuff_generic(file_input)
%% readuff_geometry 
if nargin<1 %for debug purpose only
file_input='geometry4.unv';
end
% dofXnode = 3;

fid=fopen(file_input);
dataBlocks=cell(0,2);
% typeBlocks=cell(0,1);
% dataBlocksCell
%  
dataBlockCounter=0;

validTypes=[151,164,2411,55];

while feof(fid)==0
    [dataBlock,blockType,fid] = readOneBlock(fid);
    if any(validTypes == blockType)
        dataBlockCounter=dataBlockCounter+1;     
        dataBlocks{dataBlockCounter,1}=dataBlock;
        dataBlocks{dataBlockCounter,2}=blockType;
    end
end


end


function logicValue = isaBlockDelimiterLine(tline)
       logicValue= size(tline,2)>=6 && strcmp(tline(1:6),'    -1');
end


function [dataBlock,blockType,fid] = readOneBlock(fid)
    startLine=fgetl(fid);
    %check this is the start of a new block
    if isaBlockDelimiterLine(startLine)
        typeLine = fgetl(fid);
        blockType=sscanf(typeLine(1:6),'%6u');

        dataBlock={};
        newDataLine=fgetl(fid);
        while ~isaBlockDelimiterLine(newDataLine)
            dataBlock{end+1,1}=newDataLine;
            newDataLine=fgetl(fid);
        end
    else
        % format not correct
    end

end
function read151()
end