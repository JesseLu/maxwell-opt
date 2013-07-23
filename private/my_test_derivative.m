function my_test_derivative(fun, df, f0, x0, delta)

  
    dx = delta * randn(size(x0));
    for k = 1 : 5
        dx = zeros(size(x0));
        dx(k) = delta;

        df0(k) = fun(x0+dx) - f0;
        df1(k) = real(df.' * dx);
    end
    [df0(:), df1(:)]

    norm(df0 - df1) / norm(df0);
