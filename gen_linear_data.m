% Zachary Loschinskey
% Dr. Anna Devor Rotation
% January 2024
% Generation of a test dataset for Input-Output HMM Training

% Zachary Loschinskey
% Dr. Anna Devor Rotation
% January 2024
% Generation of a test dataset for Input-Output HMM Training

%--------------------------------------------------------------------------
% Input variable x and output variable y will be used
%
% Two different dual alpha function relationships between x and y will be used
% to generate the test dataset
%
% The final dataset will be a timeseries of the input and output variables
% with defined stretches of each different linear relationship
%
% The goal is to use this dataset to help me write, verify, and fully
% understand input-output HMM expectation-maximization so that I can
% implement the IOHMM for neurovascular states relating Hb to Ca2+ through
% alpha functions
%--------------------------------------------------------------------------


function [Ca, Hb, latent, transition_probabilities, E] = gen_linear_data()
    % transition_probabilities = [0.75, 0.25;  % Probability of staying in same state or transitioning
                            %        0.12, 0.88]; % Probability of transitioning to the other state or staying
    transition_probabilities = [0.75 .25; 
        0.12 0.88];
    % Define range of x values -> constructed so that there are 70 points
    x = 0:1:100;

    % Define the two linear relationships
    E = [2.5 3; 2 1];  % m1 b1; m2 b2

    
    
    %% Generating the data set
    % Create the dataset for the IOHMM
    num_samples = 200;
    
    Ca = [];
    Hb = [];
    latent = [];
    
    % Generate the dataset
    current_state = 1; % Start in state 1
    for sample = 1:num_samples
        % Record the current state
        latent = [latent; current_state];
        
        % Generate 70 calcium sample random number between -5 and 5
        calcium_values = zeros(1, 70);
        for chim =1:70
            calcium_values(chim) = -5 + 10 * rand();
        end
    
        % Set the standard deviation for the noise
        noise_std = 0.02;
        
        % Generate random noise with the same size as calcium_values
        noise = noise_std * randn(size(calcium_values));
        
        % Add noise to each element of calcium_values
        % % % calcium_values_with_noise = calcium_values + noise;
        calcium_values_with_noise = calcium_values + noise;

        % Generate data point based on the current state and the randomly selected position
        if current_state == 1
            Hb_values = calcium_values_with_noise.*E(1,1) + E(1,2);
    
            Ca = [Ca; calcium_values]; % Record input
            Hb = [Hb; Hb_values]; % Record output
        else
            % Generate data from y1 with noise
            % Hb_values = conv(calcium_values_with_noise, dual_alpha2, 'full');
            % Hb_values = Hb_values(1:70) + offset;
            Hb_values = calcium_values_with_noise.*E(2,1) + E(2,2);

            Ca = [Ca; calcium_values]; % Record input
            Hb = [Hb; Hb_values]; % Record output
        end
        
        % Transition to next state based on transition probabilities
        if rand < transition_probabilities(current_state, current_state)
            % Stay in the current state
            continue;
        else
            % Transition to the other state
            current_state = 3 - current_state; % Flip between 1 and 2
        end
    end
    
    
    
    %% Plotting the synthetic dataset
    figure()
    
    t = linspace(0, 7 * num_samples, 70*num_samples);
    Hb_t = Hb';
    Ca_t = Ca';
    
    subplot(3, 1, 1)
    plot(t, Ca_t(:))
    xlabel('Time')
    ylabel('Calcium Conc')
    title('Calcium Timeseries')
    
    subplot(3, 1, 2)
    plot(t, Hb_t(:))
    xlabel('Time')
    ylabel('Hemoglobin Conc')
    title('Hemoglobin Timeseries')
    
    subplot(3, 1, 3)
    plot(t, repelem(latent, 70))
    xlabel('Time')
    ylabel('State')
    title('Latent Timeseries')
    ylim([0.5, 2.5])

    writematrix(Ca, "Calcium.csv")
    writematrix(Hb, "Hb.csv")
    writematrix(latent, "latent.csv")

end