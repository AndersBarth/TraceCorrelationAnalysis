function [w_res,C] = FCS_two_state_kinetics_fFCS_oneexp(fit_parameters,t,cor,sem,method)
    % fit parameters are the two relaxation times and 6 ampltidues (two each
    % for the autocorrelations and two for the crosscorrelations)
    % further, the cross-correlation starts at -1
    % 
    tau1= fit_parameters(1);    
    A11 = fit_parameters(2);
    A22 = fit_parameters(3);   
    if numel(fit_parameters) > 3 % offset
        offset = fit_parameters(4:end);
    else
        offset = zeros(1,4);
    end
    % compute ideal FCS curves
    C = zeros(numel(t),4);
    C(:,1) = A11*exp(-t./tau1) + offset(1); %DxD
    C(:,2) = A22*exp(-t./tau1) + offset(2); %AxA
    C(:,3) = -exp(-t./tau1) + offset(3); %DxA
    C(:,4) = -exp(-t./tau1) + offset(4); %AxD
    
    % compute residuals
    w_res = (C-cor)./sem; w_res(isnan(w_res)) = 0;
    if method == 1 % lsqnonlin, return array of residuals
        w_res = w_res(:);
    elseif method == 2 %%% return chi2
        w_res = sum(w_res(:).^2);
    end
end