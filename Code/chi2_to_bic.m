function BIC = chi2_to_bic(chi2,n_params,n_data_points)
% converts chi2 and n_params and n_data_points to BIC
BIC = log(n_data_points)*n_params-2*(-0.5*chi2);
end

