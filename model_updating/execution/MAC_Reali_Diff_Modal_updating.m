function [Err_pi,Err_lam,Err_V,ErrorV]= MAC_Reali_Diff_Modal_updating(...
    modeShape1,lambda1,... 
    modeShape_ref,lambda_ref ) %%coordinate punti di controllo





numberOfModes = numel(lambda_ref);
numberOfModeshapes = size(modeShape1,2);
sizeOfModesShapes = size(modeShape1,1);
%% Reali metric 
modeShape_cn     = 0*modeShape1;
modeShape_cn_ref = 0*modeShape_ref;

% Err_V = 0*ErrMatrixRel;
% Err_V = 0*modeShape1cn;
Err_V      = zeros(numberOfModes,1);
Err_lam    = zeros(numberOfModes,1);
Err_pi  = zeros(numberOfModes,1);
ErrorV   =zeros(sizeOfModesShapes,numberOfModes);
for ij=1:numberOfModeshapes
    % modeShape_ref_scale = abs(modeShape_ref(:,ij));
    % modeShape_ref_scale=modeShape_ref_scale/norm(modeShape_ref_scale);

    modeShape_cn_ref(:,ij) = modeShape_ref(:,ij)/norm(modeShape_ref(:,ij));
    modeShape_cn(:,ij) = modeShape1(:,ij)/norm(modeShape1(:,ij));
    [~,ind_maxAbsVal]=max(abs(modeShape_cn_ref(:,ij)));
    v1 = modeShape_cn(ind_maxAbsVal,ij);
    v2 = modeShape_cn_ref(ind_maxAbsVal,ij);
    if v1*v2<0
        modeShape_cn(:,ij) = -modeShape_cn(:,ij) ;
    end


    
%     Err_V(i,:) = max((modeShape1cn(:,i)- modeShape2cn(:,i)).^2);

    Err_V_ij = ( modeShape_cn(:,ij) - modeShape_cn_ref(:,ij));%.*modeShape_ref_scale;
    % Err_V_ij_neg = ( -modeShape_cn(:,ij) - modeShape_cn_ref(:,ij));

    
    ErrorV(:,ij)=Err_V_ij;
    Err_V(ij) = norm( Err_V_ij).^2;

    rel_lam_ij = (lambda1(ij)-lambda_ref(ij))/lambda_ref(ij);     
    Err_lam(ij)= sqrt(abs(rel_lam_ij))*sign(rel_lam_ij);
    % Err_lam(ij)= rel_lam_ij;


     Err_pi(ij)= (sum(Err_V_ij.^2) + abs(rel_lam_ij));
%    MAC_reali(ij)=  abs(Err_lam(ij));




end

% Err_V_ = sum((modeShape_cn-modeShape_cn_ref).^2,1);

end
