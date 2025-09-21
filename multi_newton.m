function [x, exit_flag] = multi_newton(fun,x_guess,solver_params)
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

    % numerical_diff = 0 -> fun is assumed to return [fval,J]
    max_iter = 200;
    iter = 1;
    
    
    if numerical_diff == 0
        xi = x_guess;
        [fx, J] = fun(xi);

        if round(det(J * J.'), 6) ~= 0
            xi_new = xi - (J\fx).';
        else
            iter = max_iter;
        end

        while iter < max_iter && norm(xi_new - xi) > dxtol && norm(fx) > ftol && norm(xi_new - xi) < dxmax && round(det(J * J.'), 6) ~= 0
            xi = xi_new;
            [fx, J] = fun(xi);
            if(round(det(J * J.'), 6) ~= 0)
                xi_new = xi - (J\fx).';
            end

            iter = iter + 1;
        end

        if norm(fx) > ftol
            exit_flag = 0; % success
        else
            exit_flag = 1; % fail
        end
    
    else
        % numerical_diff = 1 -> fun is assumed to only return fval        
        xi = x_guess;
        fx = fun(xi);
        J = approximate_jacobian(fun, xi);
        if(round(det(J * J.'), 6) ~= 0)
            xi_new = xi - (J\fx).';
        end

        while iter < max_iter && norm(xi_new - xi) > dxtol && norm(fx) > ftol && norm(xi_new - xi) < dxmax && round(det(J * J.'), 6) ~= 0
            xi = x_guess;
            fx = fun(xi);
            J = approximate_jacobian(fun, xi);
            if(round(det(J * J.'), 6) ~= 0)
                xi_new = xi - (J\fx).';
            end

            iter = iter + 1;
        end

        if norm(fx) > ftol
            exit_flag = 0; % success
        else
            exit_flag = 1; % fail
        end

    end

    x = xi;
    

    % NUMERICAL APPROX
    function J = approximate_jacobian(fun,x)
        J = [];
        h = 1e-6;
    
        for j = 1:length(x)
            basis_j = zeros(length(x), 1);
            basis_j(j) = 1;
            column = (fun(x + h*basis_j) - fun(x - h*basis_j)) / (2*h);
            J = [J column];
        end
    end
end


