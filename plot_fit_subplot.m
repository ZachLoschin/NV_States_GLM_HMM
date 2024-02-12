function plot_fit_subplot(ECell, real_E)
    % Number of classes
    numIterations = numel(ECell);

    % Create a figure with subplots
    figure();

    for classIndex = 1:numIterations
        % Extract Linear Fit Parameters for the current class
        currentClass = ECell{classIndex};

        % Extract slopes and intercepts for both classes
        alpha1 = currentClass(1,:);
        alpha2 = currentClass(2,:);

        % Create subplot
        subplot(numIterations, 1, classIndex);

        % Plot the lines and the data for both classes in the current subplot
        plot(alpha1, "b");
        title(sprintf("Fits after Iteration %i", classIndex));
        hold on;
        plot(alpha2, 'r')
        plot(real_E(1,:), "b--")
        plot(real_E(2,:), 'r--')
        legend("Class 1 Fit", "Class 2 Fit", "Real Class 1", "Real Class 2");
    end
end


