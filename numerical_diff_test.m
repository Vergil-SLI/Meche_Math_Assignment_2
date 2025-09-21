%function to test numerical differentiation function
function test_numerical_jacobian()
    %number of tests to perform
    num_tests = 100;
    fail = -1;

    for n = 1:num_tests
        %generate a randomized input and output dimension
        input_dim = randi([1,15]);
        output_dim = randi([1,15]);

        %generate a input_dim x input_dim x output_dim matrix stack A
        A = randn(input_dim,input_dim,output_dim);
        %generate a matrix, B of dimension output_dim x input_dim
        B = randn(output_dim,input_dim);
        %generate a column vector, C of height output_dim
        C = randn(output_dim,1);

        %create a new test function
        %this is essentially a random second-order (quadratic) function
        %with input dimension input_dim and output dimension output_dim
        test_fun = @(X) jacobian_test_function(X,A,B,C);
        X_guess = randn(input_dim,1);

        %evaluate numerical Jacobian of test_fun
        %use whatever your function name was here!
        J_numerical = approximate_jacobian(test_fun,X_guess);
        %compute the analytical jacobian of jacobian_test_function
        J_analytical = B;
        for n = 1:output_dim
            J_analytical(n,:)=J_analytical(n,:)+X_guess'*A(:,:,n);
            J_analytical(n,:)=J_analytical(n,:)+X_guess'*A(:,:,n)';
        end

        %compare with Jacobian of A
        largest_error = max(max(abs(J_numerical-J_analytical)));
        %if J is not close to A, print fail.
        if largest_error>1e-7
            fail = 1';
        end
    end

    if fail == -1
        disp('yippeeee')
    else
        disp('womp womp')
    end
end

%computes a quadratic function on input X
function f_val = jacobian_test_function(X,A,B,C)
    output_length = length(C);
    f_val = B*X+C;
    for n = 1:output_length
    f_val(n)=f_val(n)+(X'*A(:,:,n)*X);
    end
end

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