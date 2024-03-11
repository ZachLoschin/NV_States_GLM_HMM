#=
Zachary Loschinskey
Brian Depasquale Rotation
February 2024
Code for generating synthetic Ca conv IRF = Hb dataset
=# 

using CSV
using DataFrames
include("synthetic_data_utils.jl")
include("EM_utils.jl")

# -- Define parameters for synthetic data generation -- #
# Initial state Probabilities
π_data = [0.5, 0.5];

# Transition Probabilities
A_data = [0.75 0.25;
    0.12 .88];

# Define alpha function parameters
# tau1 tau2; tau1, tau2
τ_data = [0.2 0.6;
    0.6 0.8];

num_samples = 100;

# Generate the synthetic dataset
CA, HB, LATENT, IRF = gen_alpha_data(A_data, π_data, τ_data, num_samples);
#plot_data(CA, HB, LATENT)

# Storage 
num_iterations = 5;


# Generate initial guesses
A_guess = [0.8 0.2; 0.3 0.7];
π_guess = [0.5 0.5];
τ_guess = [0.3 0.4; 0.31 0.41];

IRF_guess = gen_IRF(τ_guess)

# Struct for each iteration of EM
struct Parameters
    IRF::Array{Float64, 2}
    A::Array{Float64, 2}
    π::Array{Float64, 2}
    log_lik::Float64
end

# Define vector of structs for storage
parameter_set = []


# Put initial parameters inside
push!(parameter_set, Parameters(IRF_guess, A_guess, π_guess, 0))

# EM loop
for chim=1:num_iterations

    # Get current parameters from parameter_set
    A_g = parameter_set[end].A;
    IRF_g = parameter_set[end].IRF;
    π_g = parameter_set[end].π;

    ξ_11_g, ξ_12_g, ξ_21_g, ξ_22_g, γ1_g, γ2_g, log_likelihood_g = E_step(CA, HB, A_g, IRF_g, π_g, LATENT)

    π_g, A_g, IRF_g = M_step(CA, HB, ξ_11_g, ξ_12_g, ξ_21_g, ξ_22_g, γ1_g, γ2_g)

    push!(parameter_set, Parameters(IRF_g, A_g, π_g, log_likelihood_g));

end

log_lik = [param.log_lik for param in parameter_set]
A = [param.A for param in parameter_set]
