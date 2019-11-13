E1 = 0.3;
E2 = 0.7;
dE = E2-E1;
k12 = 1; k21 = 2/3;

% theoretical relation
dE_possible = linspace(E1,E2,1000);
k21_possible = (2*E2-1)./(2*dE_possible)*(k21+k12);

%% symbolic math
% variables
syms E1 E2 k12 k21 A_DD A_AA lambda

% equations
eq_lambda = lambda == k12+k21;
dE2 = (E2-E1).^2; % delta E squared
eq_DD = A_DD == (k21*k12*dE2)./((k21*(1-E1)+k12*(1-E2)).^2);
eq_AA = A_AA == (k21*k12*dE2)./((k21*E1+k12*E2).^2);

S = solve([eq_DD,eq_AA,eq_lambda],k21);

