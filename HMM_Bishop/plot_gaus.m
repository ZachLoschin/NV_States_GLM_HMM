function plot_gaus(u, sigma)
    % Define grid for plotting
    x1_range = linspace(min(u(1) - 3*sqrt(sigma(1,1)), u(1) + 3*sqrt(sigma(1,1))), ...
                        max(u(1) - 3*sqrt(sigma(1,1)), u(1) + 3*sqrt(sigma(1,1))), 100);
    x2_range = linspace(min(u(2) - 3*sqrt(sigma(2,2)), u(2) + 3*sqrt(sigma(2,2))), ...
                        max(u(2) - 3*sqrt(sigma(2,2)), u(2) + 3*sqrt(sigma(2,2))), 100);
    [X1, X2] = meshgrid(x1_range, x2_range);
    X = [X1(:) X2(:)];
    
    % Calculate PDF values
    pdf_values = mvnpdf(X, u, sigma);
    pdf_values = reshape(pdf_values, length(x2_range), length(x1_range));
    
    % Plot contour lines
    contour(X1, X2, pdf_values, 'LineWidth', 2);
    hold on;
    
    % Plot mean
    plot(u(1), u(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
end