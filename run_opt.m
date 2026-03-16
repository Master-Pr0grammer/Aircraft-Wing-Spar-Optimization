clear all; %#ok<CLALL>
close all;

global LENGTH; %#ok<*GVMIS>
global DENSITY;
global MIN_DIST;
global YOUNGS_MODULUS;
global AIRCRAFT_MASS;
global ULTIMATE_STRENGTH;

global N_QUAD_POINTS;

% Model hyper parameters
LENGTH = 7.5; % m
DENSITY = 1600; % kg/m^3
MIN_DIST = 0.0025; % m
YOUNGS_MODULUS = 70E9; % Pa
AIRCRAFT_MASS = 500; % kg (including spar)
ULTIMATE_STRENGTH = 600E6;

N_QUAD_POINTS = 3;

% Main optimization function
function [optimized_params, mass, iterInfo] = optimize(params)
    global MIN_DIST

    % Compute constraints
    N_nodes = size(params, 1)/2;
    N_elem = N_nodes-1;
    Aineq = zeros(N_nodes, 2*(N_elem + 1));
    bineq = zeros(N_nodes, 1);

    lb = zeros(2*N_nodes,1);         % will fill the lower bounds
    ub = zeros(2*N_nodes,1);         % will fill the upper bounds
    for i = 1:N_nodes
        % min separation inequaity constraints
        Aineq(i, i*2 - 1) = 1.0;
        Aineq(i, i*2 - 0) = -1.0;

        bineq(i) = -MIN_DIST; % the upper bound value

        % --- bounds for inner radius (ri) ---------------
        lb(2*i-1)  = 0.01;         
        ub(2*i-1)  = 0.0475;       

        % --- bounds for outer radius (ro) ---------------
        lb(2*i)    = 0.0125;          
        ub(2*i)    = 0.05;     
    end

    % function to collect fmincon optimization iteration data
    iterInfo = struct('iteration', [], 'fval', [], 'constrViolation', [], 'firstOrderOpt', []); % Struct to hold data
    function stop = myOutputFcn(~, optimvalues, state)
        stop = false;
        if strcmp(state, 'iter')       % only collect during iterations
            iterInfo.iteration(end+1) = optimvalues.iteration;
            iterInfo.fval(end+1) = optimvalues.fval;

            % feasibility (constraint violation)
            if isfield(optimvalues, 'constrviolation')
                iterInfo.constrViolation(end+1) = optimvalues.constrviolation;
            else
                iterInfo.constrViolation(end+1) = NaN;
            end

            % first–order optimality
            if isfield(optimvalues, 'firstorderopt')
                iterInfo.firstOrderOpt(end+1) = optimvalues.firstorderopt;
            else
                iterInfo.firstOrderOpt(end+1) = NaN;
            end
        end
    end

    options = optimset('GradObj', 'on', 'GradConstr', 'on', 'Display', 'iter', 'Algorithm','sqp', 'OutputFcn', @myOutputFcn);
    [optimized_params, mass] = fmincon(@objective, params, Aineq, bineq, [], [], lb, ub, @noncon, options);
    fprintf("final mass: %f\n", mass);
end

%% ---------- NOMINAL DESIGN ----------
nominal_design_vars = [
    0.0415; 0.05;
    0.0415; 0.05;
    0.0415; 0.05;
    0.0415; 0.05;
    0.0415; 0.05;
    0.0415; 0.05;
    0.0415; 0.05;
    0.0415; 0.05;
    0.0415; 0.05;
    0.0415; 0.05;
    0.0415; 0.05;
    0.0415; 0.05;
    0.0415; 0.05;
    0.0415; 0.05;
    0.0415; 0.05;
];

% ---------- Plot nominal design ----------
figure
hold on
[mass, ~] = objective(nominal_design_vars);
ri = nominal_design_vars(1:2:end);   % inner radii at faces   (size nFaces)
ro = nominal_design_vars(2:2:end);   % outer radii at faces   (size nFaces)
x = linspace(0, 7.5, size(nominal_design_vars, 1)/2);
plot(x, ri, 'Color', 'red');
plot(x, ro, 'Color', 'blue');
plot(x, -ri, 'Color', 'red');
plot(x, -ro, 'Color', 'blue');
legend('inner radius', 'outer radius')
title(sprintf('Nominal Design (mass = %.2f kg)', mass))
ylabel('Distance from spar centerline (m)')
xlabel('Distance from spar root (m)')
hold off

% ---------- Plot Nominal Design force distribution ----------
figure
hold on
[stress_means, stress_stds] = spar_stress_analysis(nominal_design_vars, N_QUAD_POINTS);
N_nodes = size(nominal_design_vars, 1)/2;
x = linspace(0,7.5, N_nodes);
plot(x, stress_means, '-r');
plot(x, stress_means + 6*stress_stds, "--c");
plot(x, stress_means - 6*stress_stds, ":c");
ylabel('Stress (Pa)');
xlabel('Spar x Position (m)');
title(sprintf('Spar Nominal Design Stress Distribution (mass = %.2f kg)', mass))
legend('mean stress', 'mean + 6 std', 'mean - 6 std');
hold off





% ---------- Mesh convergence study ----------
N_nodes_test = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36];

stress_mean_values = zeros(size(N_nodes_test));
stress_std_values = zeros(size(N_nodes_test));
for i = 1:size(N_nodes_test, 2)
    N_nodes = N_nodes_test(i);

    %initialize design vars
    design_vars = zeros(N_nodes*2, 1);
    for j = 1:N_nodes
        design_vars(2*j-1) = 0.0415; 
        design_vars(2*j) = 0.05;
    end

    [stress_means, stress_stds] = spar_stress_analysis(design_vars, N_QUAD_POINTS);
    stress_mean_values(i) = stress_means(end);
    stress_std_values(i) = stress_stds(end);
end


% Plot mesh convergence
figure
hold on;
plot(N_nodes_test, stress_mean_values, '-r', 'LineWidth', 2, 'MarkerSize', 8);
plot(N_nodes_test, stress_std_values, ':c', 'LineWidth', 2, 'MarkerSize', 8);

% Add delta error bounds based on last two points (stress mean)
stress_converged = stress_mean_values(end);
total_delta = stress_mean_values(end) - stress_mean_values(end-1);
yline(stress_converged + 0.5 * total_delta, 'w--');
yline(stress_converged - 0.5 * total_delta, 'w--', 'LineWidth', 1.5);

% Add convergence line
yline(stress_converged, 'm--', 'LineWidth', 1.5, 'FontSize', 12);


xlabel('N_{nodes}');
ylabel('Stress at end of spar (Pa)');
title('Nominal Gemoetry Mesh Convergence Study');
legend('stress mean values', 'stress std values')
grid on; 
box on;

% ---------- Quadurature point convergence study ----------
N_points_test = [1, 2, 3, 4, 5, 6, 7, 8];
stress_mean_values = zeros(size(N_points_test));
stress_std_values = zeros(size(N_points_test));
for i = 1:size(N_points_test, 2)
    N_pts = N_points_test(i);

    [stress_means, stress_stds] = spar_stress_analysis(nominal_design_vars, N_pts);
    stress_mean_values(i) = stress_means(floor(size(stress_means, 1)/2));
    stress_std_values(i) = stress_stds(floor(size(stress_means, 1)/2));
end


% Plot convergence
figure
hold on;
plot(N_points_test, stress_mean_values, '-r', 'LineWidth', 2, 'MarkerSize', 8);

% Add delta error bounds based on last two points (stress mean)
stress_converged = stress_mean_values(end);
total_delta = stress_mean_values(end) - stress_mean_values(end-1);
delta_percent = abs(total_delta / stress_converged) * 100;
yline(stress_converged + 0.5 * total_delta, 'w--', sprintf('Δ%.1f%%', delta_percent), 'LineWidth', 1.5, 'FontSize', 12);
yline(stress_converged - 0.5 * total_delta, 'w--', 'LineWidth', 1.5);

% Add convergence line
yline(stress_converged, 'm--', sprintf('Converged = %.4e               ', stress_converged), 'LineWidth', 1.5, 'FontSize', 12);

xlabel('N_{quad points}');
ylabel('Mean stress at center of spar (Pa)');
title('Nominal Gemoetry Number of Quadrature Points Vs Mean Stress Analysis');
grid on; 
box on;
hold off;



% Add delta error bounds based on last two points (stress standard deviation)
figure
hold on;
plot(N_points_test, stress_std_values, ':c', 'LineWidth', 2, 'MarkerSize', 8);

stress_converged = stress_std_values(end);
total_delta = stress_std_values(end) - stress_std_values(end-1);
delta_percent = abs(total_delta / stress_converged) * 100;
yline(stress_converged + 0.5 * total_delta, 'w--', sprintf('Δ%.1f%%', delta_percent), 'LineWidth', 1.5, 'FontSize', 12);
yline(stress_converged - 0.5 * total_delta, 'w--', 'LineWidth', 1.5);

% Add convergence line
yline(stress_converged, 'm--', sprintf('Converged = %.4e               ', stress_converged), 'LineWidth', 1.5, 'FontSize', 12);

xlabel('N_{quad points}');
ylabel('Stress std. dev. at center of spar (Pa)');
title('Nominal Gemoetry Number of Quadrature Points Vs Stress Standard Deviation Analysis');
grid on; 
box on;
hold off;


%% ---------- OPTIMIZED DESIGN ---------- 
% Initialize design vars
design_vars = [
    0.015; 0.045;
    0.015; 0.045;
    0.015; 0.045;
    0.015; 0.045;
    0.015; 0.045;
    0.015; 0.045;
    0.015; 0.045;
    0.015; 0.045;
    0.015; 0.045;
    0.015; 0.045;
    0.015; 0.045;
    0.015; 0.045;
    0.015; 0.045;
    0.015; 0.045;
    0.015; 0.045;
];

% ---------- Plot initial geometry ----------
figure
hold on
[mass, dm] = objective(design_vars);
ri = design_vars(1:2:end);   % inner radii at faces   (size nFaces)
ro = design_vars(2:2:end);   % outer radii at faces   (size nFaces)
x=linspace(0, 7.5, size(design_vars, 1)/2);
plot(x, ri, 'Color', 'red');
plot(x, ro, 'Color', 'blue');
plot(x, -ri, 'Color', 'red');
plot(x, -ro, 'Color', 'blue');
legend('inner radius', 'outer radius')
title(sprintf('Initial Unoptimized Design (mass = %.2f kg)', mass))
ylabel('Distance from spar centerline (m)')
xlabel('Distance from spar root (m)')
hold off

[params, mass, optLog] = optimize(design_vars);
disp(params);

% ---------- Optimiztion convergence study ----------
if exist('optLog','var') && isfield(optLog,'iteration')
    epsv = eps;                     % tiny offset so log(0) is avoided
    figure;

    % 1st subplot – first–order optimality
    subplot(3,1,1);
    semilogy(optLog.iteration, abs(optLog.firstOrderOpt)+epsv, 'b-');
    title('First-order optimality VS Iteration');
    ylabel('Optimality');
    xlabel('Iteration');

    % 2nd subplot – feasibility
    subplot(3,1,2);
    semilogy(optLog.iteration, abs(optLog.constrViolation)+epsv, 'r-');
    title('Feasibility VS Iteration');
    ylabel('Feasibility');
    xlabel('Iteration');

    % 3rd subplot – objective value
    subplot(3,1,3);
    semilogy(optLog.iteration, abs(optLog.fval)+epsv, 'g-');
    title('Objective Value VS Iteration');
    ylabel('Objective Value')
    xlabel('Iteration');
end

% ---------- Plot final optimized geometery ---------- 
figure
hold on
ri = params(1:2:end);   % inner radii at faces   (size nFaces)
ro = params(2:2:end);   % outer radii at faces   (size nFaces)
x=linspace(0, 7.5, size(params, 1)/2);
plot(x, ri, 'Color', 'red');
plot(x, ro, 'Color', 'blue');
plot(x, -ri, 'Color', 'red');
plot(x, -ro, 'Color', 'blue');
legend('inner radius', 'outer radius')
title(sprintf('Final Optimized Design (mass = %.2f kg)', mass))
ylabel('Distance from spar centerline (m)')
xlabel('Distance from spar root (m)')
hold off

% ---------- Plot Optimal Design force distribution ----------
figure
hold on
[stress_means, stress_stds] = spar_stress_analysis(params, N_QUAD_POINTS);
N_nodes = size(params, 1)/2;
x = linspace(0,7.5, N_nodes);
plot(x, stress_means, '-r');
plot(x, stress_means + 6*stress_stds, "--c");
plot(x, stress_means - 6*stress_stds, ":c");
ylabel('Stress (Pa)');
xlabel('Spar x Position (m)');
title(sprintf('Spar Optimal Design Stress Distribution (mass = %.2f kg)', mass))
legend('mean stress', 'mean + 6 std', 'mean - 6 std');
hold off