function [mass, dm] = objective(design_vars)

    mass = calc_mass(design_vars);

    % Calclulate partial derivatives
    dm = zeros(size(design_vars));
    h = 10e-30;
    for i = 1:size(design_vars, 1)
        temp = design_vars;
        temp(i) = temp(i) + 1i*h;
        dm(i) = imag(calc_mass(temp))/h; % complex step
    end
end

% Helper function to help calculate mass (makes complex step differentiation simpler)
function [mass] = calc_mass(design_vars)
    global LENGTH %#ok<*GVMIS>
    global DENSITY

    N_elem = (size(design_vars, 1)/2) - 1;
    dx = LENGTH/N_elem;

    vol = 0;
    for i = 3:2:size(design_vars, 1)
        ri1 = design_vars(i-2);
        ro1 = design_vars(i-1);

        ri2 = design_vars(i-0);
        ro2 = design_vars(i+1);
        
        vol = vol + segment_vol(ro1, ri1, ro2, ri2, dx);
    end

    mass = vol*DENSITY;
end

% Helper function to find volume of individual segment/element
function [vol] = segment_vol(ro1, ri1, ro2, ri2, dx)
    m1 = (ro2 - ro1)/dx;
    b1 = ro1;

    m2 = (ri2 - ri1)/dx;
    b2 = ri1;

    vol = pi * ((1/3)*(m1^2 - m2^2)*dx^3 + (m1*b1 - m2*b2)*dx^2 + (b1^2 - b2^2)*dx);
end