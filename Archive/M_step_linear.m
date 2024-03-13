function [Pi_new, A_new, E_new, weightDiff, state_prediction, weights1, weights2] = M_step_linear(input, output, xi, gamma, E)
    % Calculate the pi_k initial state probability updates
    % Pi_k = gamma(Z_1k) / Sum{j}[gamma(Z1j)]
    Pi_new = [];
    
    % Calculate the new initiation probabilities
    denom = sum(exp(gamma(1,:)));
    Pi_new = exp(gamma(1,:)) ./ denom;
    
    % Calculate new transition probabilities
    
    % -- Fit the new IRFs
    [K, l] = size(E);
    N = length(input);

    A_new = zeros(2,2);
    for i = 1:K
        for j = 1:K
            A_new(i, j) = exp(log(sum(exp(xi(i,j,:)))) - log(sum(exp(gamma(1:N-1, i)))));
        end
    end

    % Normalize rows of the transition matrix
    sum_term = sum(A_new, 2)

    A_new(1,:) = A_new(1,:) ./ sum_term(1);
    A_new(2,:) = A_new(2,:) ./ sum_term(2);

    
    [N, del] = size(input);

    % Normalize the weights
    weights = exp(gamma);
    sum_term = sum(weights, 2);
    weights = weights ./ sum_term;

    % Track diff in weights, should diverge!
    weightDiff = mean(abs(weights(:,1) - weights(:,2)));
    
    input_l = input';
    output_l = output';

    weights = repelem(weights, 70, 1);
    
    % Define the linear model
    linearModel = fittype('a*x + b', 'coefficients', {'a', 'b'});

    fitObject1 = fit(input_l(:), output_l(:), linearModel, 'Weights', weights(:,1));
    fitObject2 = fit(input_l(:), output_l(:), linearModel, 'Weights', weights(:,2));
    
    slope1 = fitObject1.a;
    intercept1 = fitObject1.b;
    
    slope2 = fitObject2.a;
    intercept2 = fitObject2.b;
    
    E_new = [slope1, intercept1;
             slope2, intercept2];

    % -- Calculate each states prediction based on max weight
    w1 = zeros(length(weights(:,1)), 1) + 2;
    w1_ind = find(weights(:,1)>=weights(:,2));
    w1(w1_ind) = 1;

    weights1 = weights(:,1);
    weights2 = weights(:,2);

    state_prediction = w1;
end












