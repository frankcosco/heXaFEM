function stop = outFun(x,optimValues,state,startedTimer,res2use,namefile_str, ...
    x0,var2optimize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
persistent  counter oldresidual

% if isempty(oldresidual)
%     oldresidual=[0 0 0 0]';% optimValues.residual*0;
% end
x_sc_optimal = 0*var2optimize;
x_sc_optimal(var2optimize) = x;
x_sc_optimal(~var2optimize)=1;
x_optimal = x0.*x_sc_optimal;



stop = false;
currentTimer = toc(startedTimer);
switch state
    case 'iter'
        % Make updates to plot or guis as needed
        % currentTimer = toc(startedTimer);
        % fprintf('elapsed Time %1.2f [s]. ',currentTimer);


        numModes = sum(res2use(1,:));
        if size(res2use,1)>1

            numV     = sum(res2use(2,:));
        else
            numV = numModes;
        end
        res_lam = optimValues.residual(1:numModes)';
        if numel(optimValues.residual)>numModes
            res_V   = reshape(optimValues.residual(numModes+1:end),[],numV);


            if numV<numModes
                sqrV = 0*res_lam;
                sqrV(res2use(2,:)) = sum(res_V.^2);
            else
                sqrV =  sum(res_V.^2);
            end



            Err_PI = (sqrV+abs(res_lam).^2);
            fprintf('Err_PI: [');
            fprintf(' %1.4f ', Err_PI);
            fprintf(']. ResNorm (sum(ErrP)=)%1.7f. RMS_pi=%1.7f\n',optimValues.resnorm, rms(Err_PI));
            % % % oldresidual = [oldresidual;counter currentTimer Err_PI optimValues.resnorm];
        else
            Err_PI = (res_lam).^2;
            fprintf('Err_Pi: [');
            fprintf(' %1.4f ', Err_PI);
            % fprintf(']. ResNorm (sum(ErrP)=)%1.7f. RMS_pi=%1.7f\n',optimValues.resnorm, rms(Err_PI));
            fprintf(']. sum(ErrP)=%1.7f. ResNorm=%1.7f. RMS_pi=%1.7f\n',sum(Err_PI),optimValues.resnorm, rms(Err_PI));
            % % % oldresidual = [oldresidual; counter currentTimer res_lam optimValues.resnorm];
        end
        % % % counter=counter+1;
        fprintf('Xscaled: [');
        fprintf(' %1.4f ',x);
        fprintf(']\n');

        % fprintf('elapsed Time %1.2f [s]. Output:%s\n',currentTimer,state);
        fprintf('[%s] elapsed Time %1.2f [s]. Output:%s\n', ...
            string(datetime('now', 'Format', 'dd-MM-yyyy HH:mm')), ...
            currentTimer, state);

        if var2optimize(1),fprintf('Current density: %1.14e [kg/mm³]\n', x_optimal(1));       end
        if var2optimize(2),fprintf('Current young''s Modulus: %1.14e [kPa]\n', x_optimal(2)); end
        if var2optimize(3),fprintf('Current poisson''s Ratio: %.14f \n', x_optimal(3));       end
        if var2optimize(4),fprintf('Current thikness        : %.3f [mm]\n', x_optimal(4));    end
        if var2optimize(5),fprintf('Current width           : %.3f [mm]\n', x_optimal(5));    end
        if var2optimize(6),fprintf('Current lenght          : %.3f [mm]\n', x_optimal(6));    end
        % fprintf('Current density         : %1.14e [kg/mm³]\n', x_optimal(1));
        % fprintf('Current young''s Modulus: %1.14e [kPa]\n'   , x_optimal(2));
        % fprintf('Current poisson''s Ratio: %.14f \n'         , x_optimal(3));
        % fprintf('Current thikness        : %.3f [mm]\n'      , x_optimal(4));
        % fprintf('Current width           : %.3f [mm]\n'      , x_optimal(5));
        % fprintf('Current lenght          : %.3f [mm]\n'      , x_optimal(6));

    case 'interrupt'

    case 'init'
        % % % oldresidual =zeros(0,6);
        % % % counter=0;
    case 'done'
        % % % namefile = sprintf("%s_details.mat",namefile_str);
        % % % save(namefile,"oldresidual");

    otherwise
        % fprintf('elapsed Time %1.2f [s]. Output:%s\n',currentTimer,state);
end

end