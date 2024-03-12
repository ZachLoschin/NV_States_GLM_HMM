% Zachary Loschinskey
% Dr. Brian Depasquale Rotation
% E-step for Gaussian HMM

function [xi, gamma, log_lik] = E_step(HMM, X)
    % HMM is the HMM struct with fields A, Pi, E
    % X is 2xN matrix of observations
    
    % Define number of states and data points
    K = length(HMM.A(:,1));
    N = length(X);

    % Forward and Backward with Scaling Factors
    [a_norm, c] = forward(HMM, X);
    b_norm = backward(HMM, X, c);

    % Calculate Gamma and normalize
    gamma = a_norm .* b_norm;
    gamma = gamma ./ sum(gamma,2);


    % Calculate Xi and normalize
    xi = zeros(K,K,N-1);
    
    for t = 1:N-1
        xi_unnormalized = zeros(K, K);

        for i = 1:K
            for j = 1:K
                xi_unnormalized(i, j) = c(t) * a_norm(t, i) * HMM.A(i,j) * mvnpdf(X(t+1,:), HMM.U{j}, HMM.Sigma{j}) * b_norm(t+1,j);
            end
        end

        % Normalize the xi values
        xi(:,:,t) = xi_unnormalized / sum(xi_unnormalized, 'all');
    end
    log_lik= sum(log(c));
end

