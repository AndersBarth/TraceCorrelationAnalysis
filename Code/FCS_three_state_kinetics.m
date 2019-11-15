function [w_res,C] = FCS_three_state_kinetics(rates,t,cor,sem,method)
    % rates are the kinetic rates
    % rates = [k12,k13,k21,k23,k31,32]
    % inputs: K and time axis
%     t = 0.1:0.1:100;
%     K = [-0.32	0.32	0;... %k11 k12 k13
%     0.11	-0.5	0.39;... %k21 k22 k23
%     0.06	0.4	-0.46]; %k31 k32 k33
    if any(rates < 0)
        % can only happen during MCMC sampling
        w_res = Inf;
        return;
    end
    K = [-(rates(1)+rates(2)), rates(1), rates(2);...
        rates(3),-(rates(3)+rates(4)),rates(4);...
        rates(5),rates(6),-(rates(5)+rates(6))];
    if numel(rates) > 6
        % rest of the input parameters are offsets
        offset = rates(7:end);
    else
        offset = zeros(1,9);
    end
    K = K';
    lambda = eig(K); %eigenvalues
    
    %check for complex eigenvalues
    if any(~isreal(lambda))
        % not a valid solution, return inf
        switch method
            case 1
                % lsqnonlin, return array
                w_res = Inf(size(cor));
            case 2
                w_res = Inf;
        end
        return;
    end
    
    [U,~] = eig(K);
    [~,idx] = sort(abs(lambda));
    lambda = lambda(idx);
    U = U(:,idx);

    % get equilibrium fraction from first eigenvector
    mu = U(:,1)./sum(U(:,1));
    % adjacent eigenvectors
    V = inv(U);

    % discard first eigenvalue (not needed)
    U = U(:,2:end);
    V = V(2:end,:);

    left = U';
    right = V;
    lambda = lambda(2:end);
    
    % calculate correlation amplitudes
    G = zeros(numel(K),numel(lambda));
    s1 = [1,0,0]'; s2 = [0,1,0]'; s3 = [0,0,1]';
    G(1,:) = correlation_amplitude(s1,s1,left,right,mu); %1x1
    G(2,:) = correlation_amplitude(s2,s2,left,right,mu); %2x2
    G(3,:) = correlation_amplitude(s3,s3,left,right,mu); %3x3
    G(4,:) = correlation_amplitude(s1,s2,left,right,mu); %1x2
    G(5,:) = correlation_amplitude(s1,s3,left,right,mu); %1x3
    G(6,:) = correlation_amplitude(s2,s1,left,right,mu); %2x1
    G(7,:) = correlation_amplitude(s2,s3,left,right,mu); %2x3
    G(8,:) = correlation_amplitude(s3,s1,left,right,mu); %3x1
    G(9,:) = correlation_amplitude(s3,s2,left,right,mu); %3x2
    
    % calculate time evolution
    correct_integration_window = false;
    if correct_integration_window
        % construct a finer time axis
        sampling = 10;
        t_extended = [0;t];
        t_fine = [];
        bin = [];        
        for i = 1:(numel(t_extended)-1)
            time_sampled = linspace(t_extended(i),t_extended(i+1),sampling+1);
            t_fine = [t_fine,time_sampled(1:end-1)];
            bin = [bin,i*ones(1,sampling)];
        end
        decay = exp(repmat(lambda,[1,numel(t_fine)]).*repmat(t_fine,[numel(lambda),1]));
        C_fine = zeros(numel(t_fine),size(G,1));
        for i = 1:size(G,1)
            C_fine(:,i) = G(i,1).*decay(1,:) + G(i,2).*decay(2,:) + offset(i);
        end
        % return to coarse sampling in C
        C = zeros(numel(t),9);
        for i = 1:numel(t)
            C(i,:) = mean(C_fine(bin==i,:),1);
        end
    else
        decay = exp(repmat(lambda,[1,numel(t)]).*repmat(t',[numel(lambda),1]));
        C = zeros(numel(t),size(G,1));
        for i = 1:size(G,1)
            C(:,i) = G(i,1).*decay(1,:) + G(i,2).*decay(2,:) + offset(i);
        end
    end
    % compute residuals
    w_res = (C-cor)./sem; w_res(isnan(w_res)) = 0;
    if method == 1 % lsqnonlin, return array of residuals
        w_res = w_res(:);
    elseif method == 2 %%% return chi2
        w_res = sum(w_res(:).^2);
    end
end

function G = correlation_amplitude(s1,s2,left,right,mu)
    G = ((left*s1).*(right*(mu.*s2)))./(dot(mu,s1).*dot(mu,s2));
end
