function [HMM] = M_step(HMM, xi, gamma, X)
    % Calculate the pi_k initial state probability updates
    % Pi_k = gamma(Z_1k) / Sum{j}[gamma(Z1j)]
    K = length(HMM.A(1,:));
    N = length(X);

    % Calculate new initial state probabilities
    HMM.Pi = gamma(1,:) ./ sum(gamma(1,:));

    % Calculate new transition matrix
    A_new = zeros(K,K);
    for i = 1:K
        for j = 1:K
            %A_new(i, j) = sum(xi(i,j,:)) / (sum(xi(i,i,:)) + sum(xi(i,j,:)));
            A_new(i, j) = sum(xi(i, j, :)) / sum(gamma(i, 1:end-1));
        end
    end

    % Normalize rows of the transition matrix
    HMM.A = A_new ./ sum(A_new, 2);
    
    % Update the gaussian emission model according to Bishop pg 618
    for i = 1:K
         t = (X' * gamma(:, i)) ./ sum(gamma(:, i));
         HMM.U{i} = t';
    end

    for i = 1:K
        numerator = 0;
        for n = 1:N
            numerator = numerator + (gamma(n,i)* (X(n,:)-HMM.U{i})' * (X(n,:)-HMM.U{i}));
        end

        denominator = sum(gamma(:,i));

        HMM.Sigma{i} = 0.5 * ((numerator / denominator) + (numerator / denominator)');
    end
end

