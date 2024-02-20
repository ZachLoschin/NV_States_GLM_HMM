function plot_fit_linear(ECell, real_E)
    % Number of classes
    numIterations = numel(ECell);
    currentIndex = 1; % Initialize current index
    
    % Create figure
    fig = figure;
    
    % Plot the first set of fits
    plot_iteration(ECell{currentIndex}, real_E);
    
    % Add arrow buttons for navigation
    nextPosition = [80 20 50 20];
    prevPosition = [20 20 50 20];
    nextButton = uicontrol('Style', 'pushbutton', 'String', '>', ...
        'Position', nextPosition, 'Callback', @next_callback);
    prevButton = uicontrol('Style', 'pushbutton', 'String', '<', ...
        'Position', prevPosition, 'Callback', @prev_callback);
    
    % Callback functions for navigation buttons
    function next_callback(~, ~)
        currentIndex = min(currentIndex + 1, numIterations);
        plot_iteration(ECell{currentIndex}, real_E);
        recreate_buttons();
    end

    function prev_callback(~, ~)
        currentIndex = max(currentIndex - 1, 1);
        plot_iteration(ECell{currentIndex}, real_E);
        recreate_buttons();
    end

    function plot_iteration(currentClass, real_E)
        % Extract Linear Fit Parameters for the current class
        y1 = currentClass(1,:);
        y2 = currentClass(2,:);
        
        x = 0.1:0.1:7;

        y1_plot = y1(1).*x + y1(2);
        y2_plot = y2(1).*x + y2(2);

        y1_real = real_E(1,1).*x + real_E(1,2);
        y2_real = real_E(2,1).*x + real_E(2,2);


        % Clear current axes
        clf(fig, 'reset');
        
        % Plot the lines and the data for both classes
        plot(y1_plot, "b");
        hold on;
        plot(y2_plot, 'r')
        plot(y1_real, "b--")
        plot(y2_real, 'r--')
        title(sprintf("Fits after Iteration %i", currentIndex));
        legend("Class 1 Fit", "Class 2 Fit", "Real Class 1", "Real Class 2");
    end

    function recreate_buttons()
        % Delete existing buttons
        delete(nextButton);
        delete(prevButton);
        
        % Recreate arrow buttons
        nextButton = uicontrol('Style', 'pushbutton', 'String', '>', ...
            'Position', nextPosition, 'Callback', @next_callback);
        prevButton = uicontrol('Style', 'pushbutton', 'String', '<', ...
            'Position', prevPosition, 'Callback', @prev_callback);
        
        % Update figure
        drawnow;
    end
end
