function [w_res,C] = FCS_two_state_kinetics_colorFCS(rates,E1,E2,t,cor,sem,method)
    k12 = rates(1);
    k21 = rates(2);
    % two state filtered FCS model function
    if numel(rates) > 2 % also fitting FRET efficiencies
        % rest of the input parameters are offsets
        E1 = rates(3);
        E2 = rates(4);
    end
    
    % compute ideal FCS curves
    C = zeros(numel(t),4);
    dE2 = (E1-E2).^2; %deltaE^2
    mE = (k21*E1+k12*E2)./(k12+k21); %meanE
    xd1xd2 = (k12*k21)./((k12+k21).^2);
    C(:,1) = (dE2./((1-mE)^2))*xd1xd2*exp(-(k12+k21).*t);
    C(:,2) = (dE2./(mE^2))*xd1xd2*exp(-(k12+k21).*t);
    C(:,3) = -(dE2./(mE*(1-mE)))*xd1xd2*exp(-(k12+k21).*t);
    C(:,4) = -(dE2./(mE*(1-mE)))*xd1xd2*exp(-(k12+k21).*t);
    
    % compute residuals
    w_res = (C-cor)./sem; w_res(isnan(w_res)) = 0;
    if method == 1 % lsqnonlin, return array of residuals
        w_res = w_res(:);
    elseif method == 2 %%% return chi2
        w_res = sum(w_res(:).^2);
    end
end