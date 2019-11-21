function [w_res,C] = FCS_two_state_kinetics_fFCS_threeexp(fit_parameters,t,cor,sem,method)
    % fit parameters are the three relaxation times and 9 ampltidues (three each
    % for the autocorrelations and two for the crosscorrelations)
    % it is assumed that cross-correlation is symmetrics, i.e. A_12 = A_21
    % further, since the cross-correlation starts at -1, we can set
    % A1_12+A2_12+A3_12+A4_12 = 1
    % tau1, tau2, tau3, tau4, A1_11, A2_11, A3_11, A4_11, A1_22, A2_22, A3_22,
    % A4_22, A1_12, A2_12, A3_12
    % 
    tau1 = fit_parameters(1);
    tau2 = fit_parameters(2);
    tau3 = fit_parameters(3);
    tau4 = fit_parameters(4);
    A1_11= fit_parameters(5);
    A2_11= fit_parameters(6);
    A3_11= fit_parameters(7);
    A4_11= fit_parameters(8);    
    A1_22= fit_parameters(9);
    A2_22= fit_parameters(10);
    A3_22= fit_parameters(11);
    A4_22= fit_parameters(12);
    A1_12= fit_parameters(13);
    A2_12= fit_parameters(14);
    A3_12= fit_parameters(15);
    % compute ideal FCS curves
    C = zeros(numel(t),4);
    C(:,1) = A1_11*exp(-t./tau1)+A2_11*exp(-t./tau2)+A3_11*exp(-t./tau3)+A4_11*exp(-t./tau4); %DxD
    C(:,2) = A1_22*exp(-t./tau1)+A2_22*exp(-t./tau2)+A3_22*exp(-t./tau3)+A4_22*exp(-t./tau4); %AxA
    C(:,3) = -(A1_12*exp(-t./tau1)+A2_12*exp(-t./tau2)+A3_12*exp(-t./tau3)+(1-A1_12-A2_12-A3_12)*exp(-t./tau4)); %DxA
    C(:,4) = -(A1_12*exp(-t./tau1)+A2_12*exp(-t./tau2)+A3_12*exp(-t./tau3)+(1-A1_12-A2_12-A3_12)*exp(-t./tau4)); %AxD
    
    % compute residuals
    w_res = (C-cor)./sem; w_res(isnan(w_res)) = 0;
    if method == 1 % lsqnonlin, return array of residuals
        w_res = w_res(:);
    elseif method == 2 %%% return chi2
        w_res = sum(w_res(:).^2);
    end
end