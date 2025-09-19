% [f_val, J] = test_func1([1, 2, 3]);
% [1, 2, 3] - J\f_val;

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