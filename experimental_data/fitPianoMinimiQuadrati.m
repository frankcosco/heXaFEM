function [normale, centroide] = fitPianoMinimiQuadrati(P)
    % P è una matrice Nx3 dove ogni riga è un punto (x, y, z)
    
    % 1. Calcola il baricentro (centroide)
    centroide = mean(P, 1);
    
    % 2. Centra i punti sottraendo il baricentro
    P_centrati = P - centroide;
    
    % 3. Calcola la SVD
    % In MATLAB, [U, S, V] = svd(A) restituisce A = U*S*V'
    [~, ~, V] = svd(P_centrati, 0);
    
    % 4. Il vettore normale è l'ultima colonna di V
    % Corrisponde alla direzione di minima varianza
    normale = V(:, end);
end