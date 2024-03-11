#=
Zachary Loschinskey
Brian Depasquale Rotation
February 2024
Code for generating synthetic Ca conv IRF = Hb dataset
=# 

# alt+j, alt+r to clear workspace variables
using Plots
using SpecialFunctions
using Distributions
using FFTW


function gen_IRF(τ)

    x = 0:0.1:6.9;
    t0 = 0;
    
    # Construct the IRFs
    α1_1 = ((x .- t0) ./ τ[1, 1]) .^ 3 .* exp.(-(x .- t0) ./ τ[1, 1])
    α1_2 = -((x .- t0) ./ τ[1, 2]) .^ 3 .* exp.(-(x .- t0) ./ τ[1, 2])
    irf_1 = α1_1 + α1_2

    α2_1 = ((x .- t0) ./ τ[2, 1]) .^ 3 .* exp.(-(x .- t0) ./ τ[2, 1])
    α2_2 = -((x .- t0) ./ τ[2, 2]) .^ 3 .* exp.(-(x .- t0) ./ τ[2, 2])
    irf_2 = α2_1 + α2_2

    IRF = hcat(irf_1, irf_2);

    return IRF

end

function gen_alpha_data(A, π, τ, num_samples)
    # Define range of each timestep
    x = 0:0.1:6.9;

    # Construct time lag
    t0 = 0;

    # Construct the IRFs
    α1_1 = ((x .- t0) ./ τ[1, 1]) .^ 3 .* exp.(-(x .- t0) ./ τ[1, 1])
    α1_2 = -((x .- t0) ./ τ[1, 2]) .^ 3 .* exp.(-(x .- t0) ./ τ[1, 2])
    irf_1 = α1_1 + α1_2

    α2_1 = ((x .- t0) ./ τ[2, 1]) .^ 3 .* exp.(-(x .- t0) ./ τ[2, 1])
    α2_2 = -((x .- t0) ./ τ[2, 2]) .^ 3 .* exp.(-(x .- t0) ./ τ[2, 2])
    irf_2 = α2_1 + α2_2

    IRF = hcat(irf_1, irf_2)

    # Initialize storage for Calcium, Hemoglobin, and Latent variables
    CA = zeros(Float64, num_samples, 70)
    HB = zeros(Float64, num_samples, 70)
    LATENT = zeros(Int, num_samples)

    noise_std = 0.02;

    # Sample π for starting State
    sample = rand(Multinomial(1, π))
    latent = findfirst(x -> x == 1, sample)

    # Generate initial emission and store State
    calcium_sample = 3 .* rand(70) .- 1.5
    noise = noise_std .* randn(70)
    calcium_sample_noise = calcium_sample .+ noise

    # Generate hemoglobin values from current state IRF and noisy calcium
    hb = real.(ifft(fft(IRF[:, latent]) .* fft(calcium_sample_noise)))

    # Add values to storage
    CA[1, :] = calcium_sample
    HB[1, :] = hb
    LATENT[1] = latent

    # Loop through num_samples - 1
    for idx = 2:num_samples
        # Sample from transition multinomial Distributions
        prev_latent = LATENT[idx-1]
        sample = rand(Multinomial(1, A[prev_latent, :]))
        latent = findfirst(x -> x==1, sample)
        
        calcium_sample = 3 .* rand(70) .- 1.5
        noise = noise_std .* randn(70)
        calcium_sample_noise = calcium_sample .+ noise

        # Generate hemoglobin values from current state IRF and noisy calcium
        hb = real.(ifft(fft(IRF[:, latent]) .* fft(calcium_sample_noise)))

        # Add values to storage
        CA[idx, :] = calcium_sample
        HB[idx, :] = hb
        LATENT[idx] = latent
    end

    return CA, HB, LATENT, IRF
end

function plot_IRFS(x, α1_1, α1_2, irf_1, α2_1, α2_2, irf_2)
    # Create two subplots
    plot1 = plot(x, α1_1, label="α1_1", color=:blue)
    plot!(x, α1_2, label="α1_2", color=:red)
    plot!(x, irf_1, label="IRF1", color=:green)
    plot2 = plot(x, α2_1, label="α2_1", color=:blue)
    plot!(x, α2_2, label="α2_2", color=:red)
    plot!(x, irf_2, label="IRF2", color=:green)

    # Display the plots
    plot(plot1, plot2, layout=(2, 1), size=(800, 600))
    title!("State Specific IRF Construction")
    xlabel!("Time Point")
    ylabel!("IRF")

end

function plot_data(CA, HB, LATENT)
    ca = vec(CA');
    hb = vec(HB');

    repeated_latent = repeat(LATENT, inner=length(CA[1,:]))


    plot1 = plot(ca, legend=false)
    title!("Calcium")

    plot2 = plot(hb, legend=false)
    title!("Hemoglobin")

    plot3 = plot(repeated_latent, legend=false)
    title!("Latent")
    xlabel!("Samples")

    plot(plot1, plot2, plot3, layout=(3,1), size=(1500, 1000))
end

