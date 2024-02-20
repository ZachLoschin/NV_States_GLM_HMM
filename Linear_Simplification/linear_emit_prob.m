% Zachary Loschinskey
% Dr. Brian Depasquale Rotation
% January 2024
% Neurovascular States Project

% Function to calculate the probability of a hemoglobin output time
% series based on the Ca input and the current linear function parameters.
% Used in the forward and backward algorithms for alpha function
% Expecation maximization.

function likelihood = linear_emit_prob(calcium, real_Hb, E)
    % Function to calculate the log probability of emitting Hb given Ca
    % and the current state IRF

    % E is the m,b of the current state guess

    % Predict the hemoglobin wave from the current IRF and Ca input
    pred_Hb = E(1).*calcium + E(2);

    % Calculate the log likelihood -> normalized residuals
    log_likelihood = -mean(abs(pred_Hb - real_Hb));
    likelihood = exp(log_likelihood);

    % -- Plot for visualization --%
    % figure();
    % plot(pred_Hb, 'r');
    % hold on;
    % plot(real_Hb, 'b');
    % legend("Predicted Hb", "Current Guess IRF")
    
end



