% Zachary Loschinskey
% Anna Devor Rotation
% January 2024
% Neurovascular States Project

% Function to calculate the probability of a hemoglobin output time
% series based on the Ca input and the current alpha function parameters.
% Used in the forward and backward algorithms for alpha function
% Expecation maximization.
function likelihood = alpha_emit_prob(calcium, real_Hb, IRF)
    % Function to calculate the log probability of emitting Hb given Ca
    % and the current state IRF

    % Predict the hemoglobin wave from the current IRF and Ca input
    pred_Hb = ifft(fft(calcium).*fft(IRF));

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



