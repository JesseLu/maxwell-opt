function my_gradient_test(fun, grad_x, x0, use_real, text)
   
    delta = 1e-6;
    f0 = fun(x0); 
    for k = 1 : 10
        dx = delta * (randn(size(x0)) + 0i * randn(size(x0)));
        diff1 = fun(x0+dx) - f0;
        diff2 = (grad_x * dx);
        if use_real
            diff2 = real(diff2);
        end
        err(k) = norm(diff1(:) - diff2(:)) / norm(diff1(:));
    end
    fprintf('Max error (%s): %e\n', text, max(err));
