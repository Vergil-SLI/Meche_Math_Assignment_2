function class_2025_09_23_reference()
    X0 =  randn(3,1);

    [~,J_analysis] = test_function01(X0);
    J_numerical = approximate_jacobian01(@test_function01,X0);

    solver_params = struct();
    solver_params.dxmin = 1e-10;
    solver_params.ftol  = 1e-10;
    solver_params.dxmax = 1e8;
    solver_params.max_iter = 200;
    solver_params.approx = 1;

    Xguess = randn(3,1);

    X_root = multivar_newton(fun,Xguess,solver_params);

    disp(X_root)

    f_root = test_function01(X_root);

    disp(f_root)

    % disp(J_analysis)
    % disp(J_numerical)
end

function X_root = multivar_newton(fun,X,solver_params)
    dxmin = solver_params.dxmin;
    ftol = solver_params.ftol;
    dxmax = solver_params.dxmax;
    max_iter = solver_params.max_iter;
    approx = solver_params.approx;


        delta_x = 1;

        count = 0;

    while count<max_iter && norm(delta_x) > dxmin && norm(fval) > ftol && norm(delta_x) < dxmax
        count = count+1;
        if approx
            fval = fun(X);
            J = approximate_jacobian01(fun,X);
        else
            [fval,J] = fun(X);
        end

        delta_x = -J/fval;

        X = X + delta_x;

    end
        X_root = X;

end

function J = approximate_jacobian01(fun,X)
    f0 = fun(X);

    J = zeros(length(f0),length(X)); %height is outputs of f width is # of inputs of f

    en = zeros(length(X),1);

    delta_X = 1e-6;

    for n = 1:length(X)
        en(n) = 1;
        
        fplus = fun(X+en*delta_X);
        fminus = fun(X-en*delta_X);

        J(:,n) = (fplus-fminus)/(2*delta_X);

        en(n) = 0;
    end
end

function [fval,J] = test_function01(X)
    x1 = X(1); x2 = X(2); x3 = X(3);
    
    f1 = x1^2 + x2^2 -6-x3^5;
    f2 = x1*x3+x2-12;
    f3 = sin(x1+x2+x3);

    fval = [f1;f2;f3];

    J = [[2*x1,2*x2,-5*x3^4];[x3, 1, x1];[cos(x1+x2+x3),cos(x1+x2+x3),cos(x1+x2+x3)]];
end 