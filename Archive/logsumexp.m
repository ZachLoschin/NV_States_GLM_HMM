function result = logsumexp(x)
    % Function to compute the log(sum(exp(x))) avoiding overflow/underflow

    c = max(x); % Find the maximum value
    result = c + log(sum(exp(x - c)));
end