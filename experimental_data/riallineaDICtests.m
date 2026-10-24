function ExpMNodes=riallineaDICtests()
addpath(genpath('..\'))

%% Experimental data
[ExpV,ExpModal_freqs,ExpNodes] = readuff_modalshapes('real_modeset.unv');

Punti = ExpNodes(:,2:4)*1000;


% Beam dimensions in mm
a = 40; % Length (L)
b = 360;  % Width (W)
c = 4;   % Thickness (T)

Punti(:,1)=Punti(:,1)+a/2;
Punti(:,2)=Punti(:,2)-b/2;

[ ~,resortingInd]=sort(Punti(:,1),"descend");
PuntiR=Punti(resortingInd,:);



% --- Calcolo del piano ---
[n, centro] = fitPianoMinimiQuadrati(PuntiR);
%%
fprintf('Vettore Normale: [%.2e, %.2e, %.2e]\n', n);
fprintf('Baricentro: [%.4f, %.4f, %.4f]\n', centro);
%%
figure(1),
% subplot(2,1,1),
cla,
plot(Punti(:,2),Punti(:,1),'*-')
% axis equal
text(Punti(:,2),Punti(:,1),string((1:size(Punti,1))'))



%%


ImRef = imread("referenceDIC-FE-UPDATING.tif");
ImRef2 = imread("referenceDic-grid.tif");
ImRef2 = ImRef2(:,:,1:3);
% figure(13),cla,
% imagesc(ImRef2)
% ar=drawpolygon();




figure(2),
% subplot(2,1,1),cla,imshow(ImRef2(:,:,1:3));
Im=adapthisteq(ImRef,"ClipLimit",0.01 ,"NBins",64);
% subplot(2,1,2),
cla,imshow(Im)
%
roiBeamPixel=[47.5197300103842	16.3504672897195
1807.5	40
1799.	227
59.4823468328141	204.429387331256];

% roiBeamPixel=[42.5353063343718	15.0212876427829;
%            48.5166147455867	135.976635514019;
%      1177.32242990654	151.594496365524;
%      1179.98078920042	32.3006230529594];



hold on,
plot(roiBeamPixel(:,1),roiBeamPixel(:,2),'*r')

roiBeamPhysical=[40 -360;
                 40  0;
                 0   0;
                 0   -360];


%%

% Calcolo della matrice omografica
H = calcolaOmografia(roiBeamPhysical,roiBeamPixel);


PuntiDemo = applicaOmografia(PuntiR(:,1:2), H);

plot(PuntiDemo(:,1),PuntiDemo(:,2),'r*')
PuntiDemoOffs = PuntiDemo+repmat([-8 +6],size(PuntiDemo,1),1);
plot(PuntiDemoOffs(:,1),PuntiDemoOffs(:,2),'y+')
%%
roiRefs=[42.0430463576159	16.0741710296684
1182.22004816376	30.9606525911708
1175.40367248645	151.585514834206
49.4121011438892	135.743746364165];

Hrefs = calcolaOmografia(roiBeamPhysical,roiRefs);

PuntiRefs = applicaOmografia(PuntiR(:,1:2), Hrefs);

figure(3),
cla
imshow(ImRef)
% plot(PuntiDemo(:,1),PuntiDemo(:,2),'b*')
hold on,

PGrid=grigliaVincolata(70,41,15);
% plot(PGrid(:,1),PGrid(:,2),'.r')
% nq =size(PGrid,1);
% quiver(PuntiDemo(1:nq,1),PuntiDemo(1:nq,2), ...
%     PGrid(:,2)-PuntiDemo(1:nq,1), ...
%     PGrid(:,1)-PuntiDemo(1:nq,2),"off")

Indices = orderdGrid();

H_final=H;

[H_final] = fineTuneHomography(PuntiR(:,1:2), PGrid, H);

% PGridH = applicaOmografia(Punti(:,1:2), H_final);




Hinv = calcolaOmografia(roiBeamPixel,roiBeamPhysical);
% PgridCorrected= applicaOmografia(PGridH, Hinv);


% PuntiC = ExpNodes(:,2:4)*1000;

% 
% % Beam dimensions in mm
% a = 40; % Length (L)
% b = 360;  % Width (W)
% c = 4;   % Thickness (T)
% PuntiC(:,1)=Punti(:,1)+a/2;
% PuntiC(:,2)=Punti(:,2)-b/2;

PGridH = applicaOmografia(Punti(:,1:2), H_final);
plot(PGridH(:,1),PGridH(:,2),'.m')

PgridCorrected= applicaOmografia(PGridH, Hinv);


figure(4),cla,
plot(roiBeamPhysical([1 2 3 4 1],2),roiBeamPhysical([1 2 3 4 1],1),'-r'), hold on
plot(PgridCorrected(:,2),PgridCorrected(:,1),'.b')


ExpMNodes=[ExpNodes(:,1) PgridCorrected];
ExpMNodes(:,4)=-c;

axis equal




%%
% % p1=detectHarrisFeatures(ImRef,"ROI",[roiBeam2D(1,[1 2])-5 11 11],"FilterSize",7);
% % plot(p1);
% % p1loc = p1.Location;
% % 
% % p4=detectHarrisFeatures(ImRef,"ROI",[roiBeam2D(4,[1 2])-5 11 11],"FilterSize",7);
% % plot(p4)
% % p4loc = p4.Location;


end

function Q_stima = applicaOmografia(P, H)
    % P: Matrice Nx2 dei punti da trasformare
    N = size(P, 1);
    % Aggiungi la colonna di 1 per le coordinate omogenee
    P_omogenei = [P, ones(N, 1)];
    
    % Trasformazione: Q = H * P'
    res = (H * P_omogenei')';
    
    % Ritorno a coordinate euclidee (divisione per la terza coordinata w)
    Q_stima = res(:, 1:2) ./ res(:, 3);
end

function H = calcolaOmografia(P, Q)
    % P: Matrice Nx2 dei punti sorgente [x, y]
    % Q: Matrice Nx2 dei punti target [u, v]
    
    N = size(P, 1);
    if N < 4
        error('Servono almeno 4 coppie di punti.');
    end
    
    % Costruzione della matrice A per il sistema Ah = 0
    A = zeros(2*N, 9);
    for i = 1:N
        x = P(i,1); y = P(i,2);
        u = Q(i,1); v = Q(i,2);
        
        % Ogni coppia di punti fornisce due equazioni
        A(2*i-1, :) = [-x, -y, -1,  0,  0,  0, u*x, u*y, u];
        A(2*i, :)   = [ 0,  0,  0, -x, -y, -1, v*x, v*y, v];
    end
    
    % Risoluzione tramite SVD (il vettore nullo di A)
    [~, ~, V] = svd(A);
    h = V(:, end); % L'ultima colonna di V è la soluzione
    
    % Rimodella in matrice 3x3
    H = reshape(h, 3, 3)';
    
    % Normalizzazione (imposta l'ultimo elemento a 1)
    H = H / H(3,3);
end

function PGrid=grigliaVincolata(u0,v0,s,Indices)
arguments
    u0 (1,1) double=65;
    v0 (1,1) double=40;
    s (1,1) double=15;
    Indices =orderdGrid();
end
PGrid = [u0+s*Indices(:,2) v0+s*Indices(:,1)];

end
function Indices=orderdGrid()
Indices=[];
redStripe=(22:-1:0);
Indices=[Indices; 0*redStripe'  redStripe'];

orangeStripe=(115-24:-1:0);
Indices=[Indices; 1+0*orangeStripe'  orangeStripe']; %#ok<*NASGU>

fullStripe=(114:-1:0);
Indices=[Indices; 2+0*fullStripe'  fullStripe']; 
Indices=[Indices; 3+0*fullStripe'  fullStripe']; 
Indices=[Indices; 4+0*fullStripe'  fullStripe']; 
Indices=[Indices; 5+0*fullStripe'  fullStripe']; 
 
fullStripe=(114:-1:1);
Indices=[Indices; 6+0*fullStripe'  fullStripe']; 
Indices=[Indices; 7+0*fullStripe'  fullStripe']; 
Indices=[Indices; 8+0*fullStripe'  fullStripe']; 
Indices=[Indices; 9+0*fullStripe'  fullStripe']; 
 
fullStripe=(113:-1:35);
Indices=[Indices; 10+0*fullStripe'  fullStripe']; 
 
fullStripe=(113:-1:109);
Indices=[Indices; 11+0*fullStripe'  fullStripe']; 
% 




end


function [H_final] = fineTuneHomography(P_phys, Pgrid, H_guess)
    % P_phys: [N x 2] Coordinate sul piano fisico (es. metri o mm)
    % Indices: [N x 2] Indici interi della griglia (es. [0,0; 0,1; 1,0...])
    % H_guess: [3 x 3] Stima iniziale dell'omografia
    
    % 1. Preparazione parametri iniziali
    % Normalizziamo H_guess affinché H(3,3) = 1 e prendiamo i primi 8 elementi
    h_init = H_guess(1:8) / H_guess(9);
    
    % Stima iniziale per la griglia (u0, v0, s)
    % Proviamo a proiettare i punti con H_guess per avere un'idea della scala s
    Q_stima = applicaOmografia(P_phys, H_guess./H_guess(9));
    u_init = Q_stima(:,1);
    v_init = Q_stima(:,2);
    

    u0_guess = min(u_init);
    v0_guess = min(v_init);
    s_guess = 15;
    
    % Vettore x = [h1...h8, u0, v0, s]
    x0 = h_init(:)';
    % % lb = [-inf(1,8),0,0,5];
    % ub = [inf(1,8),500,500,50];
    % 2. Ottimizzazione non lineare
    options = optimoptions('lsqnonlin', ...
        'Display', 'iter', ...
        'Algorithm', 'levenberg-marquardt', ...
        'FunctionTolerance', 1e-9, ...
        'MaxIterations', 400);
    
    x_opt = lsqnonlin(@(x) residui_vincolati(x, P_phys, Pgrid), x0, [],[], options);
    
    % 3. Estrazione dei risultati finali
    H_final = [x_opt(1:3)', x_opt(4:6)',[ x_opt(7:8)'; 1]];
    % grid_params.u0 = x_opt(9);
    % grid_params.v0 = x_opt(10);
    % grid_params.s  = x_opt(11);
end

function R = residui_vincolati(x, P_phys, P_target)
    % Ricostruzione H (8 DoF)
    H = [x(1:3)', x(4:6)',[ x(7:8)'; 1]];
    
    % Parametri griglia (3 DoF)
    % u0 = x(9); v0 = x(10); s = x(11);
    
    % Proiezione proiettiva
     Q_stima = applicaOmografia(P_phys, H);
    u_p = Q_stima(:,1);
    v_p = Q_stima(:,2);
   
    
    % Calcolo posizioni teoriche sulla griglia equispaziata
    % u = u0 + i*s, v = v0 + j*s
    u_target = P_target(:,1) ;
    v_target = P_target(:,2);
    
    % Residui (differenza tra proiezione e griglia ideale)
    R = [u_p - u_target; v_p - v_target];
end
