%% FFT-Accelerated Subdiffusive fBm in the Plane with Gamma-Mixed Exponential Stopping
% Each trajectory has its own λ ~ Gamma(α, β), then τ ~ Exp(λ)
% Stopping locations and radial distances are recorded; histogram produced.
% Modified to include heavy-tail diagnostics for the distribution of radial distances r_stop,
% following recommended methods: log-log survival plot, mean excess plot, Hill plot,
% and Hill estimator for tail index α.

clear; close all; clc;
rng(42); % For reproducibility

% Simulation parameters
H = 0.35; % Hurst exponent (0 < H < 0.5 for subdiffusion)
N = 2^17; % Number of time steps (N+1 points); power-of-2 friendly
T = 50000; % Total simulation horizon (should be >> typical 1/λ)
M = 5000; % Number of independent trajectories

% Gamma prior parameters for the exponential rate λ
alpha = 0.5; % shape parameter (α > 0)
beta = 0.5; % rate parameter (mean(λ) = α/β, var(λ) = α/β²)

% Preallocate results
stopping_locations = zeros(M, 2);
r_stop = zeros(M, 1);
lambda_realized = zeros(M, 1); % for diagnostics
tau_realized = zeros(M, 1);

fprintf('Generating %d trajectories of subdiffusive fBm (H = %.2f)\n', M, H);
fprintf('Exponential stopping with λ ~ Gamma(α=%.1f, β=%.1f) → E[λ]=%.3f\n', ...
        alpha, beta, alpha/beta);

for i = 1:M
    % Sample rate parameter λ from Gamma(α, β)
    lambda_i = gamrnd(alpha, 1/beta); % MATLAB uses scale = 1/rate
    
    % Draw stopping time τ ~ Exp(λ_i)
    tau = -log(rand(1)) / lambda_i;
    if tau > T
        tau = T; % rare truncation when tail is heavy
    end
    
    % Generate full 1D paths (only once per trajectory)
    [X_full, t] = simulate_fbm_daviesharte(H, N, T);
    [Y_full, ~] = simulate_fbm_daviesharte(H, N, T);
    
    % Linear interpolation → position at τ
    X_tau = interp1(t, X_full, tau, 'linear');
    Y_tau = interp1(t, Y_full, tau, 'linear');
    
    stopping_locations(i, :) = [X_tau, Y_tau];
    r_stop(i) = sqrt(X_tau^2 + Y_tau^2);
    lambda_realized(i) = lambda_i;
    tau_realized(i) = tau;
end

fprintf('Completed.\n');
fprintf(' Mean λ realized = %.4f (theoretical = %.4f)\n', mean(lambda_realized), alpha/beta);
fprintf(' Mean τ realized = %.4f (theoretical = %.4f)\n', mean(tau_realized), beta/alpha);
fprintf(' Mean stopping distance = %.4f\n', mean(r_stop));

%% Visualization
% Figure 1: Scatter of all stopping positions
figure('Position', [100 100 700 600]);
scatter(stopping_locations(:,1), stopping_locations(:,2), 8, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('X(\tau)');
ylabel('Y(\tau)');
title(sprintf('Stopping positions (M = %d, H = %.2f, λ ∼ Gamma(%.1f, %.1f))', ...
              M, H, alpha, beta));
axis equal;
grid on; box on;set(gca, 'FontSize', 12);

% Figure 2: Histogram of radial distances at stopping
figure('Position', [100 100 800 500]);
histogram(r_stop, 40, 'Normalization', 'pdf', 'FaceColor', [0 0.4470 0.7410], ...
'EdgeColor', 'none');
xlabel('Distance r from origin at stopping time \tau');
ylabel('Probability density');
title(sprintf('Histogram of stopping distances (M = %d)', M));
grid on; box on;set(gca, 'FontSize', 12);

% Optional extra figure: distribution of realized λ and τ
figure('Position', [100 100 900 400]);
subplot(1,2,1);
histogram(lambda_realized, 30, 'Normalization','pdf');
xlabel('Realized λ'); ylabel('Density'); title('Distribution of λ ∼ Gamma(α,β)');
subplot(1,2,2);
histogram(tau_realized, 40, 'Normalization','pdf');
xlabel('Realized τ'); ylabel('Density'); title('Resulting mixture-of-exponentials τ');

%% Heavy-Tail Diagnostics for Radial Distances r_stop
% Sort the data in descending order for tail analysis
r_sorted = sort(r_stop, 'descend'); % Upper order statistics

% 1. Log-log survival plot (empirical CCDF on log-log scale)
figure('Position', [100 100 800 500]);
loglog(r_sorted, (1:M)/M, 'b.', 'MarkerSize', 10);
xlabel('log(r)');
ylabel('log(P(R > r))');
title('Log-log Survival Plot for Radial Distances');
grid on; box on;set(gca, 'FontSize', 12);

% Interpretation: Straight line in the right tail suggests power-law heavy tails.

% 2. Mean excess plot
thresholds = r_sorted(2:floor(M/10)); % Use upper 10% as potential thresholds, skipping max
mean_excess = zeros(length(thresholds), 1);
for idx = 1:length(thresholds)
    u = thresholds(idx);
    exceedances = r_stop(r_stop > u) - u;
    mean_excess(idx) = mean(exceedances);
end
figure('Position', [100 100 800 500]);
plot(thresholds, mean_excess, 'r-',  'LineWidth', 1.4);
xlabel('Threshold u');
ylabel('Mean Excess e(u)');
title('(e) Mean Excess Plot for Radial Distances');
grid on; box on;set(gca, 'FontSize', 12);
% Interpretation: Upward trend suggests heavy tails (Pareto-like).

% 3. Hill plot for tail index estimation
k_max = floor(M/2); % Max number of upper order stats to consider
alpha_hat = zeros(k_max, 1);
for k = 2:k_max
    alpha_hat(k) = k / sum(log(r_sorted(1:k)) - log(r_sorted(k)));
end
figure('Position', [100 100 800 500]);
plot(2:k_max, alpha_hat(2:end),'Color', [0 0.6 0], 'LineWidth', 1.4);
xlabel('Number of upper order statistics k');
ylabel('Hill estimator \alpha(k)');
title('Hill Plot for Tail Index');
grid on; box on; set(gca, 'FontSize', 12);
% Interpretation: Stable plateau indicates reliable \alpha estimate; small \alpha (<2-3) suggests heavy tails.

% Compute and display a simple Hill estimate (using k = sqrt(M) as heuristic)
k_heuristic = floor(sqrt(M));
alpha_est = k_heuristic / sum(log(r_sorted(1:k_heuristic)) - log(r_sorted(k_heuristic)));
fprintf('Estimated tail index \alphâ (using k=%d): %.4f\n', k_heuristic, alpha_est);
% Rule-of-thumb: \alphâ < 2 → very heavy tails (infinite variance); \alphâ < 1 → infinite mean.

% Additional heuristic: sample kurtosis (if finite)
if all(isfinite(r_stop))
    kurt = kurtosis(r_stop);
    fprintf('Sample kurtosis: %.4f (>>3 suggests heavy tails)\n', kurt);
else
    fprintf('Kurtosis undefined (infinite values present).\n');
end

function [B, t] = simulate_fbm_daviesharte(H, N, T)
    dt = T / N;
    t = (0:N) * dt;
    n = N;
    m = 2^nextpow2(2 * n);
    r = zeros(n, 1);
    r(1) = 1;
    for k = 1:n-1
        r(k+1) = 0.5 * ((k+1)^(2*H) - 2*k^(2*H) + (k-1)^(2*H));
    end
    c = zeros(m, 1);
    c(1:n) = r;
    c(m-n+2:m) = r(n:-1:2);
    lambda = real(fft(c));
    lambda = max(lambda, 1e-12);
    Z = randn(m, 1) + 1i * randn(m, 1);
    fGn = real(ifft(sqrt(lambda) .* Z));
    fGn = fGn(1:n);
    B_unit = cumsum([0; fGn]);
    B = B_unit * (dt ^ H);
end