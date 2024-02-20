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


function [Ca, Hb, latent, y1, y2, transition_probabilities, E_real] = gen_linear_data()
    % transition_probabilities = [0.75, 0.25;  % Probability of staying in same state or transitioning
                            %        0.12, 0.88]; % Probability of transitioning to the other state or staying
    transition_probabilities = [0.7 .3; 
                                0.3 0.7];
    % Define range of x values -> constructed so that there are 70 points
    x = 0.1:0.1:7;
    
    % Define parameters of two linear classess
    m1 = 5;
    b1 = 0.5;

    m2 = -3;
    b2 = 0;

    E_real = [m1 b1; m2 b2];

    % Define classes
    y1 = m1.*x + b1;
    y2 = m2.*x + b2;

    
    % --- Plotting to verify synthetic data creation --- %
    figure()
    plot(x, y1, 'r')
    hold on;
    plot(x, y2, 'b')
    title("Linear Class Functions")
    xlabel("x")
    ylabel("y")
    legend("Class 1", "Class 2")
    
    
    %% Generating the data set
    % Create the dataset for the IOHMM
    num_samples = 100;
    
    Ca = [];
    Hb = [];
    latent = [];
    
    % Generate the dataset
    current_state = 1; % Start in state 1
    for sample = 1:num_samples
        % Record the current state
        latent = [latent; current_state];
        
        % Generate 70 calcium sample random number between -1.5 and 1.5
        calcium_values = zeros(1, 70);
        for chim =1:70
            calcium_values(chim) = -1.5 + 3 * rand();
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
            % Generate data from y1 with noise
            Hb_values = m1.*calcium_values_with_noise+b1;
    
            Ca = [Ca; calcium_values]; % Record input
            Hb = [Hb; Hb_values]; % Record output
        else
            % Generate data from y1 with noise
            Hb_values = m2.*calcium_values_with_noise+b2;

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
    
    t = linspace(0, 7 * num_samples, 7000);
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

end
    
    
    
    

