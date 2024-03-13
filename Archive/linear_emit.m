
function log_likelihood = linear_emit(x, u, E)

    % Extract the slope and intercept from E
    slope = E(1);
    intercept = E(2);

    % Calculate the expected value of x based on u
    u_expected = slope*x + intercept;

    % Calculate residuals as negative log likelihood
    log_likelihood = -mean(abs(u - u_expected));

end