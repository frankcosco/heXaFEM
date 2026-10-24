% Tall = parseModelUpdatingFolder("..\lastRun\fixed-t", "MU-Singlerun-TRRBeam*.txt");
% Tall = parseModelUpdatingFolder("..\lastRun\fixed-nu", "MU-Singlerun-TRRBeam*.txt");

% Tall = parseModelUpdatingFolder("..\revisedFrequencies", "MU-Singlerun-TRRBeam*.txt");
[Tall, TallHistory] = parseModelUpdatingFolder("revisedFrequencies\materialOnly", "MU-Singlerun-TRRBeam*.txt");
% Tall = parseModelUpdatingFolder("revisedFrequencies\allParams2", "MU-Singlerun-TRRBeam*.txt");

% close all,

%% 2) Mappa stili coerente con i tuoi plot
styleMap = containers.Map;
styleMap('Hexa8-Int8')   = 'hexa8';
styleMap('Hexa20-Int27') = 'hexa20';
styleMap('Hexa27-Int27') = 'hexa27';
styleMap('Hexa64-Int64') = 'hexa64';

orderList = ["Hexa8-Int8","Hexa20-Int27","Hexa27-Int27","Hexa64-Int64"];

%% Figure 7 - RMS Iniziale

txtRMSinit = buildTikzCoordinatesFromTable( ...
    Tall, ...
    'DOFs_total', ...
    'RMS_init', ...
    'formulation', ...
    orderList, ...
    styleMap, ...
    4,true,1);

txtRMSinit = buildTikzCoordinatesFromTable( ...
    Tall, ...
    'DOFs_total', ...
    'ResNorm_init', ...
    'formulation', ...
    orderList, ...
    styleMap, ...
    4,true,2);

txtRMSinit = buildTikzCoordinatesFromTable( ...
    Tall, ...
    'DOFs_total', ...
    'ResNorm_opt', ...
    'formulation', ...
    orderList, ...
    styleMap, ...
    4,true,4);




%% Figure 9 - RMS finale
txtRMSopt = buildTikzCoordinatesFromTable( ...
    Tall, ...
    'DOFs_total', ...
    'RMS_opt', ...
    'formulation', ...
    orderList, ...
    styleMap, ...
    4, true,3);

%% figure 8-a
xtRho = buildTikzCoordinatesFromTable( ...
    Tall, 'DOFs_total',   'rho_opt', 'formulation', orderList, styleMap, 9,false,5);

txtY = buildTikzCoordinatesFromTable( ...
    Tall, 'DOFs_total',   'Y_opt', 'formulation', orderList, styleMap, 0,false,6);

txtNU = buildTikzCoordinatesFromTable( ...
    Tall, 'DOFs_total',   'nu_opt', 'formulation', orderList, styleMap, 4,true,7);

% txtT = buildTikzCoordinatesFromTable( ...
%     Tall, 'DOFs_total',   't_opt', 'formulation', orderList, styleMap, 4,false,8);

% txtT = buildTikzCoordinatesFromTable( ...
%     Tall, 'DOFs_total',   'b_opt', 'formulation', orderList, styleMap, 4,false);

% buildTikzCoordinatesFromTable( ...
%     Tall, 'time_s',   'ResNorm_opt', 'formulation', orderList, styleMap, 4,false,10);

%% figure 10
selected = ismember(TallHistory.seriesLabel, [
     "Hexa8-Int8-52x468x4"
     "Hexa8-Int8-65x585x5"
    % "Hexa8-Int8_78x702x6"
    "Hexa20-Int27-26x234x2"
    "Hexa27-Int27-26x234x2"
    "Hexa64-Int64-13x117x1"
]);
%
figure(10);clf, hold on; grid on;

labels = unique(TallHistory.seriesLabel(selected),'stable');

for i = 1:numel(labels)
    idx = selected & TallHistory.seriesLabel == labels(i);

    plot(TallHistory.iter(idx), TallHistory.gap_norm(idx), '-o', ...
        'DisplayName', labels(i));
end

xlabel('Time [s]');
ylabel('$(J_k-J^*)/(J_0-J^*)$', 'Interpreter','latex');
legend('Location','best');
set(gca,"XScale","linear","YScale","log")

%%
%%
selected = ismember(TallHistory.seriesLabel, [
    "Hexa8-Int8-52x468x4"
    "Hexa20-Int27-26x234x2"
    "Hexa27-Int27-26x234x2"  
    "Hexa8-Int8-39x351x3"
    "Hexa64-Int64-13x117x1"

    % 
    % "Hexa8-Int8-78x702x6"
    % "Hexa20-Int27-39x351x3"
    % "Hexa27-Int27-39x351x3"
    % "Hexa64-Int64-26x234x2"
]);
figure(10);clf, hold on; grid on;

% labels = unique(TallHistory.seriesLabel(selected),'stable');

labels = [
    "Hexa8-Int8-52x468x4"
    "Hexa20-Int27-26x234x2"
    "Hexa27-Int27-26x234x2"
    "Hexa8-Int8-39x351x3"
    "Hexa64-Int64-13x117x1"
    ];
for i = 1:numel(labels)
    idx = selected & TallHistory.seriesLabel == labels(i);
    time_s_i = TallHistory.time_s(idx);
    gap_norm_i  = TallHistory.gap_norm(idx);
    gapStart = gap_norm_i(1);

    plot_i = plot(time_s_i, gap_norm_i , '-o', ...
        'DisplayName', labels(i));

    colorPlot(i,:)=plot_i.Color;              %#ok<SAGROW>
    lastPoint(i,:)=[time_s_i(end),gap_norm_i(end)];%#ok<SAGROW>

end

xlabel('Time [s]');
ylabel('$J_k/J^*$', 'Interpreter','latex');
legend(NumColumns=2,Location="northeast");
set(gca,"XScale","linear","YScale","lin")
ylim([0.9 3.3])
box on

for i = 1:numel(labels)
p=plot([lastPoint(i,1) lastPoint(i,1)],[1.00 lastPoint(i,2)],Color=colorPlot(i,:), ...
    Marker='pentagram', MarkerSize=8,LineWidth=2,DisplayName='Best',LineStyle='none');
set( get(get(p,'Annotation'), 'LegendInformation'),'IconDisplayStyle','off')

end

%%
selected = ismember(TallHistory.seriesLabel, [
    %  "Hexa8-Int8-52x468x4"
    % "Hexa20-Int27-26x234x2"
    % "Hexa27-Int27-26x234x2"  
    % 
    % "Hexa64-Int64-13x117x1"

    % 
    "Hexa8-Int8-78x702x6"
    "Hexa20-Int27-39x351x3"
    "Hexa27-Int27-39x351x3"
    "Hexa64-Int64-26x234x2"
]);
figure(11);clf, hold on; grid on;

labels = unique(TallHistory.seriesLabel(selected),'stable');
labels = [
    "Hexa8-Int8-78x702x6"
    "Hexa20-Int27-39x351x3"
    "Hexa27-Int27-39x351x3"
    "Hexa64-Int64-26x234x2"
    ];


for i = 1:numel(labels)
    idx = selected & TallHistory.seriesLabel == labels(i);
    time_s_i = TallHistory.time_s(idx);
    gap_norm_i  = TallHistory.gap_norm(idx);

    plot_i = plot(time_s_i,gap_norm_i , '-*', ...
        'DisplayName', labels(i));

    colorPlot(i,:)=plot_i.Color;              %#ok<SAGROW>
    lastPoint(i,:)=[time_s_i(end),gap_norm_i(end)];%#ok<SAGROW>

end

xlabel('Time [s]');
ylabel('$J_k/J^*$', 'Interpreter','latex');
legend(NumColumns=1,Location="northeast");
set(gca,"XScale","linear","YScale","lin")
ylim([0.9 1.9])
box on

for i = 1:numel(labels)

p=plot([lastPoint(i,1) lastPoint(i,1)],[1.00 lastPoint(i,2)],Color=colorPlot(i,:), ...
    Marker='pentagram', MarkerSize=8,LineWidth=2,DisplayName='',LineStyle='none');
set( get(get(p,'Annotation'), 'LegendInformation'),'IconDisplayStyle','off')

% p=plot(lastPoint(i,1),lastPoint(i,2),Color=colorPlot(i,:), ...
%     );
% set( get(get(p,'Annotation'), 'LegendInformation'),'IconDisplayStyle','off')



end
% ylim([0.95 3])
% txtL = buildTikzCoordinatesFromTable( ...
%     Tall, 'DOFs_total',   'L_opt', 'formulation', orderList, styleMap, 4,false);
% 
% txtW = buildTikzCoordinatesFromTable( ...
%     Tall, 'DOFs_total',   'b_opt', 'formulation', orderList, styleMap, 4,false);
