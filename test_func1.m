function [f_val,J] = test_func1(X)
    f1 = X(1)^2 + X(2)^2 - X(3)^5 -6;
    f2 = X(1)*X(3) + X(2) - 12;
    f3 = sin(X(1) + X(2) + X(3));

    f_val = [f1; f2; f3];
    J = [2*X(1), 2*X(2), 5*X(3)^4;
        X(3), 1, X(1);
        cos(X(1) + X(2) + X(3)), cos(X(1) + X(2) + X(3)), cos(X(1) + X(2) + X(3))]   ;
end