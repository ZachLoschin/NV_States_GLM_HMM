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


function [Ca, Hb, latent, alpha1, alpha2] = gen_alpha_data()
    transition_probabilities = [0.75, 0.25;  % Probability of staying in same state or transitioning
                                    0.12, 0.88]; % Probability of transitioning to the other state or staying
    
    % Define range of x values -> constructed so that there are 70 points
    x = 0.1:0.1:7;
    
    % Define parameters of the two dual alpha functions
    alpha1_t0 = 0;
    alpha1_tau1 = 0.5;
    alpha1_tau2 = 0.7;

    % alpha1_t0 = 0;
    % alpha1_tau1 = 0.5;
    % alpha1_tau2 = 0.6;
    
    alpha1_1 = (((x-alpha1_t0) ./ alpha1_tau1) .^3) .* exp(-(x-alpha1_t0) ./ alpha1_tau1);
    alpha1_2 = -(((x-alpha1_t0) ./ alpha1_tau2) .^3) .* exp(-(x-alpha1_t0) ./ alpha1_tau2);
    dual_alpha1 = alpha1_1 + alpha1_2;
    alpha1 = dual_alpha1;
    
    % Define parameters of the two dual alpha functions
    alpha2_t0 = 0;
    alpha2_tau1 = 0.1;
    alpha2_tau2 = 0.8;
    % 
    % alpha2_t0 = 0;
    % alpha2_tau1 = 0.4;
    % alpha2_tau2 = 0.55;
    
    alpha2_1 = (((x-alpha2_t0) ./ alpha2_tau1) .^3) .* exp(-(x-alpha2_t0) ./ alpha2_tau1);
    alpha2_2 = -(((x-alpha2_t0) ./ alpha2_tau2) .^3) .* exp(-(x-alpha2_t0) ./ alpha2_tau2);
    dual_alpha2 = alpha2_1 + alpha2_2;
    alpha2 = dual_alpha2;
    
    
    figure()
    subplot(2, 1, 1)
    plot(x, alpha1_1, 'r')
    hold on;
    plot(x, alpha1_2, 'b')
    plot(x, dual_alpha1, 'k')
    title("Construction of Class 1 Dual Alpha Function")
    legend("Component1", "Component2", "Dual Alpha 1")
    
    subplot(2, 1, 2)
    plot(x, alpha2_1, 'r')
    hold on;
    plot(x, alpha2_2, 'b')
    plot(x, dual_alpha2, 'k')
    title("Construction of Class 2 Dual Alpha Function")
    legend("Component1", "Component2", "Dual Alpha 2")
    
    
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
        noise_std = 0.0;
        
        % Generate random noise with the same size as calcium_values
        noise = noise_std * randn(size(calcium_values));
        
        % Add noise to each element of calcium_values
        % % % calcium_values_with_noise = calcium_values + noise;
        calcium_values_with_noise = calcium_values + noise;
        
        % if sample ==1
        %     offset = 0;
        % else
        %     offset = Hb(end);
        % end

        offset = 0;


        % Generate data point based on the current state and the randomly selected position
        if current_state == 1
            % Generate data from y1 with noise
            % Hb_values = conv(calcium_values_with_noise, dual_alpha1, 'full');
            % Hb_values = Hb_values(1:70) + offset;
            Hb_values = ifft(fft(dual_alpha1) .* fft(calcium_values_with_noise));
    
            Ca = [Ca; calcium_values]; % Record input
            Hb = [Hb; Hb_values]; % Record output
        else
            % Generate data from y1 with noise
            % Hb_values = conv(calcium_values_with_noise, dual_alpha2, 'full');
            % Hb_values = Hb_values(1:70) + offset;
            Hb_values = ifft(fft(dual_alpha2) .* fft(calcium_values_with_noise));

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
    
    
    
    

