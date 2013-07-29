function my_gradient_test(fun, grad_x, x0, type, text)
   
    delta = 1e-6;
    f0 = fun(x0);
    for k = 1 : 10
        dx = delta * (randn(size(x0)));
        f1 = fun(x0+dx);
        diff1 = f1 - f0;
        diff2 = (grad_x * dx);
        switch type
            case 'real'
                diff2 = real(diff2);
            case 'complex'
                diff2 = diff2;
            otherwise
                error('Invalid type.');
        end
        err(k) = norm(diff1(:) - diff2(:)) / norm(diff1(:));
    end
    fprintf('Max error (%s): %e\n', text, max(err));
