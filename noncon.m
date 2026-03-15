% nonlinear constraint function for fmincon to enforce stress constraints
function [cineq, ceq, jineq, jeq] = noncon(x)
    global ULTIMATE_STRENGTH %#ok<*GVMIS>
    global N_QUAD_POINTS;

    % calc constraint Helper function
    function [cineq] = calc_constraint(x, n)
        [stress_means, stress_stds] = spar_stress_analysis(x, n);
        obj = stress_means + 6 * stress_stds;
        cineq = obj./ULTIMATE_STRENGTH - 1;
    end

    N_nodes = size(x, 1)/2;
    cineq = calc_constraint(x, N_QUAD_POINTS);

    % Calclulate partial derivatives
    jineq = zeros(size(x, 1), N_nodes);

    h = 10e-30;
    for i = 1:size(x, 1)
        temp = x;
        temp(i) = temp(i) + 1i*h;
        jineq(i, :) = (imag(calc_constraint(temp, N_QUAD_POINTS))/h)'; % complex step
    end

    ceq = [];
    jeq = [];
end