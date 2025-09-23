function [xi, exit_flag] = multi_newton(fun,x_guess,solver_params)
    %unpack values from struct (if fields in struct have been set)
    dxtol = 1e-14; %terminate early if |x_{i+1}-x_{i}|<dxtol
    if isfield(solver_params,'dxtol')
        dxtol = solver_params.dxtol; 
    end

    ftol = 1e-14; %terminate early if |f(x_i)|<ftol
    if isfield(solver_params,'ftol')
        ftol = solver_params.ftol; 
    end

    dxmax = 1e8; %terminate early if |x_{i+1}-x_{i}|>dxmax
    if isfield(solver_params,'dxmax')
        dxmax = solver_params.dxmax; 
    end

    numerical_diff = 1; %1 = numeric diff, 0 = analytical
    if isfield(solver_params,'numerical_diff')
        numerical_diff = solver_params.numerical_diff; 
    end
    
    % ACTUAL NEWTON FUNCTION CODE
    max_iter = 250;
    iter = 1;
    xi_prev = x_guess;
    
    xi = x_guess;
    % calculate your first 2 inital points (second one needed to see if
    % change in x is already too small)
    if numerical_diff == 0
        [fx, J] = fun(xi_prev);
    else
        fx = fun(xi_prev);
        J = approximate_jacobian_for_newton(fun, xi_prev);
    end
    
    delta_x = -J\fx;
    
    if det(J * J.') ~= 0
        xi = xi_prev - J\fx;
    else
        xi = xi_prev + ftol*10; % here so program doesn't error
    end
    
    % keep finding "next point" until either change too small, too many
    % iter, or find the root
    while iter < max_iter && norm(xi - xi_prev) > dxtol && norm(fx) > ftol && norm(xi - xi_prev) < dxmax && det(J * J.') ~= 0
        xi_prev = xi;
        
        if numerical_diff == 0
            [fx, J] = fun(xi);
        else
            fx = fun(xi);
            J = approximate_jacobian_for_newton(fun, xi);
        end

        if det(J * J.') ~= 0
            xi = xi - J\fx;
        end
        
        iter = iter + 1;
    end
    
    distance_from_zero = norm(fx)
    
    % different exit flags
    if iter == max_iter
        exit_flag = "too many iterations";
    elseif norm(xi - xi_prev) < dxtol && norm(fx) < 100*ftol
        exit_flag = "close enough to success";
    elseif norm(xi - xi_prev) > dxmax
        exit_flag = "too big of a change in x";
    elseif det(J * J.') == 0
        exit_flag = "determinant is 0";
    elseif norm(fx) < ftol    
        exit_flag = "success"; % yippee
    else
        exit_flag = 'termination occured for another reason... how';
    end
    

end



% Numerical Approximation of the Jacobian
function J = approximate_jacobian_for_newton(fun,x)
    J = [];
    h = 1e-6;

    for j = 1:length(x)
        basis_j = zeros(length(x), 1);
        basis_j(j) = 1;
        column = (fun(x + h*basis_j) - fun(x - h*basis_j)) / (2*h);
        J = [J, column];
    end
end