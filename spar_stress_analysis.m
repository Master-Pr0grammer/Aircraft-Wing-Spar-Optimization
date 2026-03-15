%calculate stresses along beam during convergence studies
function [stress_means, stress_stds] = spar_stress_analysis(params, N_quad_points)
    N_nodes = size(params, 1)/2;

    %Calculate Force Distribution
    x = linspace(0,7.5, N_nodes);
    [means, stds] = compute_stress_stats(x, params, N_quad_points);

    stress_means = means;
    stress_stds = stds;
end

function [stress_means, stress_stds] = compute_stress_stats(x, params, N_quad_points)

    [xi, wts] = GaussHermite_Locations_Weights(N_quad_points);
    wts = wts./sqrt(pi);


    % Compute Sigma values
    mu = [0; 0; 0; 0];
    sigma = [0; 0; 0; 0];
    for n = 1:4
        sigma(n) = fnom(0)/(10*n);
    end

    Ef = 0;
    Ef2 = 0;
    for i1 = 1:size(xi, 1)
        pt1 = sqrt(2) * sigma(1) * xi(i1) + mu(1);
        for i2 = 1:size(xi, 1)
            pt2 = sqrt(2) * sigma(2) * xi(i2) + mu(2);
            for i3 = 1:size(xi, 1)
                pt3 = sqrt(2) * sigma(3) * xi(i3) + mu(3);
                for i4 = 1:size(xi, 1)
                    pt4 = sqrt(2) * sigma(4) * xi(i4) + mu(4);

                    c = [pt1; pt2; pt3; pt4];
                    stresses = compute_stresses(x, params, c);

                    Ef = Ef + wts(i1).*wts(i2).*wts(i3).*wts(i4).* stresses;
                    Ef2 = Ef2 + wts(i1).*wts(i2).*wts(i3).*wts(i4).* stresses.^2;
                end
            end
        end
    end

    stress_means = Ef;
    stress_stds = sqrt(Ef2 - stress_means.^2);

end

function [val] = fnom(x)
    global LENGTH %#ok<*GVMIS>
    global AIRCRAFT_MASS

    weight = AIRCRAFT_MASS * 9.81; % <--------------------------------------- NOTE: MAY NEED TO REMOVE 0.5!!!!!

    val = (2.5*weight)/LENGTH .* (1-x./LENGTH);
end

function [stresses] = compute_stresses(x, params, c)
    global YOUNGS_MODULUS
    global LENGTH

    N_nodes = size(params, 1)/2;
    N_elem = N_nodes - 1;

    %--- split radii -----------------------------------------------
    ri = params(1:2:end);   % inner radii at faces   (size nFaces)
    ro = params(2:2:end);   % outer radii at faces   (size nFaces)
    Iyy = (pi*( (2*ro).^4 - (2*ri).^4 ))/64;

    tol = 2.7028741483E-9;
    for i = 1:size(Iyy, 1)
        if real(Iyy(i)) < tol
            Iyy(i) = tol;
        end
    end

    Fx = f(x, c);
    u = CalcBeamDisplacement(LENGTH, YOUNGS_MODULUS, Iyy, Fx, N_elem);
    stresses = CalcBeamStress(LENGTH, YOUNGS_MODULUS, ro, u, N_elem);
end

function [val] = f(x, c)
    val = fnom(x) + delta_f(x, c);
end

function [val] = delta_f(x, c)
    global LENGTH
    sum = 0;
    for n = 1:4
        sum = sum + (c(n) .* cos( ((2*n - 1).*pi.*x) ./ (2*LENGTH)) );
    end
    val = sum;
end