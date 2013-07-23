function test_structgrad()

    % Create a problem
    A0 = randn(10) + 1i * randn(10);
    A = @(z) A0 + diag(z(:));
    b = randn(10, 1) + 1i * randn(10, 1);
    B = @(x) diag(x(:));
    d = @(x) b - A0 * x;

    % Test equivalence of (A, b) and (B, d).
    x = randn(10, 1) + 1i * randn(10, 1);
    z = randn(10, 1);

    c1 = A(z)*x - b;
    c2 = B(x)*z - d(x);

    norm(c1 - c2) / norm(c1)

    % Fitness function.
    fx = @(x) 0.5 * norm(x)^2;
    f = @(z) fx(A(z)\b);
    dfx = @(x) x;
   
    grad_test(fx, conj(dfx(x)), x)


    % Okay, now let's try to find df/dz.
    delta = 1e-6;

    % Brute force.
    z0 = z;
    f0 = f(z0);
    for k = 1 : length(z)
        z1 = z0;
        z1(k) = z1(k) + delta;
        f1(k) = f(z1);
        dfdz(k) = (f1(k) - f0) / delta;
        dxdz(:,k) = 1/delta * ((A(z1)\b) - (A(z0)\b));
    end
    dfdz;
    grad_test(f, dfdz(:), z0)
        
    % Gradient test for dxdz.
    dxdz;
    grad_test2(@(z)A(z)\b, dxdz', z0)
    x0 = A(z0) \ b; 
    grad_test2(@(z)A(z)\b, (-A(z0)\B(x0))', z0)
    % Now using adjoint method.
    dfdz_adj = real((-A(z0)' \ dfx(x0))' * B(x0));

    [dfdz(:), dfdz_adj(:)]

function grad_test(fun, df, x0)
    
    delta = 1e-6;
    f0 = fun(x0); 
    for k = 1 : 10
        dx = delta * randn(size(x0));
        actual = fun(x0+dx) - f0;
        grad = real(df' * dx);
        err(k) = norm(actual - grad) / norm(actual);
    end

    fprintf('Max error: %e\n', max(err));

function grad_test2(fun, df, x0)
    
    delta = 1e-6;
    f0 = fun(x0); 
    for k = 1 : 10
        dx = delta * randn(size(x0));
        actual = fun(x0+dx) - f0;
        grad = (df' * dx);
        err(k) = norm(actual(:) - grad(:)) / norm(actual(:));
    end

    fprintf('Max error: %e\n', max(err));
