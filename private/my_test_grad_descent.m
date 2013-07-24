function [] = my_test_grad_descent(n)

    A = randn(n);
    x0 = randn(n, 1);
    function [fval, dx] = fun(x)
        fval = 0.5 * norm(A*x(:)).^2;
        dx = A' * A * x(:);
    end

    [x, fval, hist] = my_grad_descent(@fun, x0);

    semilogy(hist, '.-')

end
