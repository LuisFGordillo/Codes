%% SIMULATE_BRANCHING_DISPERSION  
%   Multi-generation branching process with
%   dispersal according to the stopped subdiffusive fBm radial density g(r).
%
%   - Starts with 5 initial "parents" (Generation 0) sampled from g(r).
%   - Each individual produces a Poisson(1.2) number of offspring.
%   - Each offspring disperses from its parent according to g(r).

function simulate_branching_dispersion()

    close all; clc;
    % Parameters
    n_initial = 5;          
    lambda_poisson = 1.2;      
    num_generations = 20;     
    a = 0.5;                 
    b = 0.5;                 
    sigma = 1.0;             
    H = 0.35;                
    T_max = 0.05;            % <<< KEY PARAMETER: controls overall spatial scale
                             %    Smaller T_max → smaller spread 

    % Seed for reproducibility
    %rng(42);

    % Storage
    all_X = [];
    all_Y = [];
    all_Gen = [];

    % ====================== Generation 0 ======================
    [X0, Y0] = sample_stopped_positions(n_initial, a, b, sigma, H, T_max);
    all_X = [all_X; X0];
    all_Y = [all_Y; Y0];
    all_Gen = [all_Gen; zeros(n_initial, 1)];

    current_X = X0;
    current_Y = Y0;

    % ====================== Generations 1–num_generations ======================
    for gen = 1:num_generations
        n_current = length(current_X);
        new_X = [];
        new_Y = [];
        
        for i = 1:n_current
            n_children = poissrnd(lambda_poisson);
            if n_children > 0
                [dX, dY] = sample_stopped_positions(n_children, a, b, sigma, H, T_max);
                new_X = [new_X; current_X(i) + dX];
                new_Y = [new_Y; current_Y(i) + dY];
            end
        end
        
        all_X = [all_X; new_X];
        all_Y = [all_Y; new_Y];
        all_Gen = [all_Gen; gen * ones(length(new_X), 1)];
        
        current_X = new_X;
        current_Y = new_Y;
    end

    % ====================== Plotting (Blue Gradient) ======================
    figure('Position', [100 100 900 700]);
    hold on;
    
    n_colors = num_generations + 1;
    dark_blue  = [0.0, 0.0, 0.6];
    light_blue = [0.7, 0.85, 1.0];
    
    colors = zeros(n_colors, 3);
    for i = 1:n_colors
        t = (i-1)/(n_colors-1);
        colors(i, :) = (1-t)*dark_blue + t*light_blue;
    end
    
    for g = 0:num_generations
        idx = (all_Gen == g);
        if any(idx)
            scatter(all_X(idx), all_Y(idx), 25, colors(g+1, :), 'filled', ...
                    'MarkerEdgeColor', [0.2 0.2 0.2], 'LineWidth', 0.4, ...
                    'DisplayName', sprintf('Generation %d', g));
        end
    end
    hold off;
    
    axis equal;
    grid off;
    title(sprintf('Branching Dispersal Process'), ...
          'FontSize', 14);
    xlabel('x');
    ylabel('y');
    legend('Location', 'bestoutside');
    
    fprintf('Simulation complete (T_max = %.4f):\n', T_max);
    for g = 0:num_generations
        n = sum(all_Gen == g);
        fprintf('  Generation %d: %d individuals\n', g, n);
    end
end

% ====================== Helper Function ======================
function [X, Y] = sample_stopped_positions(n, a, b, sigma, H, T_max)
% Sample n 2D positions consistent with the stopped fBm model.
% T_max rescales the stopping time to match desired spatial scale.
    if nargin < 6
        T_max = 0.05;               
    end
    
    lambda = gamrnd(a, b, n, 1);
    tau    = exprnd(1 ./ lambda, n, 1);
    
    tau_scaled = tau * T_max;                    
    scale_t    = sigma .* (tau_scaled .^ H);
    
    R = raylrnd(scale_t, n, 1);
    theta = 2 * pi * rand(n, 1);
    
    X = R .* cos(theta);
    Y = R .* sin(theta);
end