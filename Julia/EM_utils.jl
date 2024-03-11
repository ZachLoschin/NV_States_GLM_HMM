#=
Zachary Loschinskey
Brian Depasquale Rotation
February 2024
=#
using Statistics
using DSP
using FFTW


function conv_emit_prob(input, output, IRF)
    
    # Takes in Ca, Hb, and a state IRF and gives emission probabilities
    # Convolve input and IRF to get pred_output
    pred_output = ifft(fft(input).*fft(IRF));

    # Calculate residuals
    residuals = broadcast(abs, output .- pred_output);

    # Calculate neg log likelihood
    neg_log_likelihood = -mean(residuals);

    # Convert to a likelihood
    likelihood = exp(neg_log_likelihood);

    return likelihood
end


function forward(INPUT, OUTPUT, IRF, A, π)
    # Input     - NxL matrix of input variable wehere N is number of timesteps and L is lenght of timestep
    # Output    - NxL matrix of output variable 
    # A         - KxK matrix of transition probabilities w/ K states
    # IRFS      - KxN matrix of current state specific IRFs
    # π         - Kx1 matrix of initial state probabiliites
    
    # Get parameters
    N,l = size(INPUT);
    k = length(A[1,:]);

    # NxK storage for posterior state probabilities
    F_norm = zeros(Float64, N, k);
    F_current = zeros(Float64, 1, k);
    C = zeros(1, N);

    # First forward probability
    for state = 1:k
        F_current[state] = π[state] * conv_emit_prob(INPUT[1,:], OUTPUT[1,:], IRF[:, state]);
    end

    C[1] = sum(F_current);
    F_norm[1,:] = F_current ./ C[1];

    # For the rest of the timesteps
    for n = 2:N
        F_current = zeros(Float64, 1, k);

        for state = 1:k
            for prev_state = 1:k
                F_current[prev_state] = F_norm[n-1, prev_state] * A[prev_state, state];
            end
            F_norm[n, state] = conv_emit_prob(INPUT[n,:], OUTPUT[n,:], IRF[:,state]) * sum(F_current);
        end

        # Normalizing factor
        C[n] = sum(F_norm[n,:]);
        F_norm[n,:] = F_norm[n,:] ./ C[n];
    end

    return F_norm, C
end


function backward(INPUT, OUTPUT, IRF, A, π, C)
    # Input     - NxL matrix of input variable wehere N is number of timesteps and L is lenght of timestep
    # Output    - NxL matrix of output variable 
    # A         - KxK matrix of transition probabilities w/ K states
    # IRFS      - KxN matrix of current state specific IRFs
    # π         - Kx1 vector of initial state probabiliites
    # C         - 1xN vector of scaling factors from forward recursion

    # Get parameters
    N,l = size(INPUT);
    k = length(A[1,:]);

    # NxK storage for posterior state probabilities
    B_norm = zeros(Float64, N, k);
    B_current = zeros(Float64, 1, k);

    # Initialize backward probabiliites
    for state = 1:k
        B_norm[N, state] = 0.5
    end

    for n = N-1:-1:1
        for state = 1:k
            B_current = zeros(Float64, 1, k);

            for next_state = 1:k
                B_current[next_state] = B_norm[n+1, next_state] * A[state, next_state] * conv_emit_prob(INPUT[n+1,:], OUTPUT[n+1,:], IRF[:, next_state])
            end
            # Normalize using scaling factors from forward algorithm
            B_norm[n, state] = sum(B_current) / C[n+1];
        end
    end

    return B_norm
end


function E_step(INPUT, OUTPUT, A, IRF, π, LATENT)
    # Calculate Forward/Backward Probabilities
    α_norm, C = forward(INPUT, OUTPUT, IRF, A, π);
    β_norm = backward(INPUT, OUTPUT, IRF, A, π, C);

    # Calculate poserior probabilities of the states and Normalize
    γ1 = α_norm[:,1] .* β_norm[:,1];
    γ2 = α_norm[:,2] .* β_norm[:,2];

    # Normalize the γ values
    sum_term = γ1 .+ γ2;
    γ1 = γ1 ./ sum_term;
    γ2 = γ2 ./ sum_term;

    # Caclulate ξ values
    N, l = size(INPUT);

    ξ_11 = zeros(Float64, N-1, 1);
    ξ_12 = zeros(Float64, N-1, 1);
    ξ_21 = zeros(Float64, N-1, 1);
    ξ_22 = zeros(Float64, N-1, 1);

    for n = 1:N-1
        ξ_11[n] = (α_norm[n, 1] * A[1,1] * conv_emit_prob(INPUT[n+1,:], OUTPUT[n+1,:], IRF[:,1]) * β_norm[n+1, 1]) / C[n+1];
        ξ_12[n] = (α_norm[n, 1] * A[1,2] * conv_emit_prob(INPUT[n+1,:], OUTPUT[n+1,:], IRF[:,2]) * β_norm[n+1, 2]) / C[n+1];
        ξ_21[n] = (α_norm[n, 2] * A[2,1] * conv_emit_prob(INPUT[n+1,:], OUTPUT[n+1,:], IRF[:,1]) * β_norm[n+1, 1]) / C[n+1];
        ξ_22[n] = (α_norm[n, 2] * A[2,2] * conv_emit_prob(INPUT[n+1,:], OUTPUT[n+1,:], IRF[:,2]) * β_norm[n+1, 2]) / C[n+1];
    end

    # Normalize ξ values
    sum_term = ξ_11 .+ ξ_12 .+ ξ_21 .+ ξ_22;

    # Normalize after rounding
    ξ_11 = ξ_11 ./ sum_term
    ξ_12 = ξ_12 ./ sum_term
    ξ_21 = ξ_21 ./ sum_term
    ξ_22 = ξ_22 ./ sum_term


    # Calculate log_likelihood
    log_likelihood = log(prod(C));

    return ξ_11, ξ_12, ξ_21, ξ_22, γ1, γ2, log_likelihood
end


function M_step(INPUT, OUTPUT, ξ_11, ξ_12, ξ_21, ξ_22, γ1, γ2)

    # Update π
    π_new = zeros(Float64, 1, 2);
    sum_term = γ1[1] + γ2[1];
    π_new[1] = γ1[1] / sum_term;
    π_new[2] = γ2[1] / sum_term;

    # Update A
    A_new = zeros(Float64, 2, 2);
    A_new[1,1] = sum(ξ_11) / (sum(ξ_11) + sum(ξ_12));
    A_new[1,2] = sum(ξ_12) / (sum(ξ_11) + sum(ξ_12));
    A_new[2,1] = sum(ξ_21) / (sum(ξ_21) + sum(ξ_22));
    A_new[2,2] = sum(ξ_22) / (sum(ξ_21) + sum(ξ_22));

    # Update the IRFs
    N,l = size(INPUT);  # N is timesteps, l is samples in a timestep
    PRED_IRF = zeros(Float64, N, l);

    for idx = 1:N
        PRED_IRF[idx, :] = real(ifft(fft(OUTPUT[idx,:]) ./ fft(INPUT[idx,:])));
    end
    
    IRF_NEW = zeros(Float64, l, 2);

    IRF_NEW[:,1] = PRED_IRF' * γ1;
    IRF_NEW[:,2] = PRED_IRF' * γ2;

    return π_new, A_new, IRF_NEW

end
