function [w_res,C] = FCS_two_state_kinetics_fFCS(rates,t,cor,sem,method)
    k12 = rates(1);
    k21 = rates(2);
    % two state filtered FCS model function
    if numel(rates) > 2
        % rest of the input parameters are offsets
        offset = rates(3:end);
    else
        offset = zeros(1,4);
    end
    
    % compute ideal FCS curves
    C = zeros(numel(t),4);
    C(:,1) = (k21/k12)*exp(-(k12+k21).*t)+offset(1); %DxD
    C(:,2) = (k12/k21)*exp(-(k12+k21).*t)+offset(2); %AxA
    C(:,3) = -exp(-(k12+k21).*t)+offset(3); %DxA
    C(:,4) = -exp(-(k12+k21).*t)+offset(4); %AxD
    
    % compute residuals
    w_res = (C-cor)./sem; w_res(isnan(w_res)) = 0;
    if method == 1 % lsqnonlin, return array of residuals
        w_res = w_res(:);
    elseif method == 2 %%% return chi2
        w_res = sum(w_res(:).^2);
    end
end