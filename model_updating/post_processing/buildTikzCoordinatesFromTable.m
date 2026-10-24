function [tikzText] = buildTikzCoordinatesFromTable( ...
    Tall, xVar, yVar, groupVar, orderList, styleMap, precision, addLegend,figNum)
% buildTikzCoordinatesFromTable
%
% Tall       : table con tutti i dati
% xVar       : nome variabile asse x, es. 'DOFs_free'
% yVar       : nome variabile asse y, es. 'RMS_init'
% groupVar   : nome variabile gruppo, es. 'formulation'
% orderList  : ordine desiderato dei gruppi, es. ["Hexa8-Int8","Hexa20-Int27","Hexa27-Int27","Hexa64-Int64"]
% styleMap   : containers.Map gruppo -> stile tikz, es. hexa8, hexa20, ...
% precision  : cifre decimali per y
% addLegend  : boolean - (true: aggiunge legendEntry)
%
% Output:
%   tikzText : stringa completa da copiare in LaTeX

    if nargin <8 
        addLegend=false;
    end
    if nargin < 7 || isempty(precision)
        precision = 4;
    end

    if nargin < 6 || isempty(styleMap)
        styleMap = containers.Map;
        for i = 1:numel(orderList)
            styleMap(char(orderList(i))) = char(orderList(i));
        end
    end

    tikzLines = strings(0,1);

    % dataValCell={};
    if ~isempty(figNum)
        figure(figNum)
        set(gcf,"Name",sprintf("%d-%s",figNum,yVar),"NumberTitle","off","WindowStyle","docked")
    else
        figure(Name=yVar,NumberTitle="off",WindowStyle="docked"),
    end
    
    gca,cla, hold on,

    for i = 1:numel(orderList)
        g = orderList(i);

        idx = strcmp(string(Tall.(groupVar)), g);
        T = Tall(idx,:);

        if isempty(T)
            continue;
        end

        T = sortrows(T, xVar);

        coords = strings(height(T),1);
        xVal=[];
        yVal=[];
        for k = 1:height(T)
            x = T.(xVar)(k);
            y = T.(yVar)(k);

            if isnan(x) || isnan(y)
                continue;
            end

            coords(k) = sprintf('(%g,%.*f)', x, precision, y);
            xVal = [xVal x];
            yVal = [yVal y];
        
        end
        
        % yVal = [yVal ]

        coords = coords(coords ~= "");

        styleName = string(styleMap(char(g)));

        block = [
            string(sprintf('%% %s', g));
            sprintf('\\addplot[%s] coordinates {', styleName);
            strcat('    ',sprintf('%s',coords));
            '};';           
        ];
        if addLegend
            block=[block; sprintf('\\addlegendentry{%s}',g)]; %#ok<*AGROW>
        end



        tikzLines = [tikzLines; string(block)];


        plot(xVal,yVal, '-*','DisplayName',g);

            
        % dataValCell{end+1,1}=xVal;
        % dataValCell{end,  2}=yVal;


    end

    tikzText = strjoin(tikzLines, newline);
    fprintf('%s\n', tikzText);
end