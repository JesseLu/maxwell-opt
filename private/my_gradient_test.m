function my_gradient_test(fun, grad_x, x0, type, text)
   
    delta = 1e-6;
    f0 = fun(x0);
    for k = 1 : 10
        dx = delta * randn(size(x0));
        if strcmp(type, 'real_with_imag')
            dx = dx + 1i * delta * randn(size(x0));
        end
        f1 = fun(x0+dx);
        diff1 = f1 - f0;
        diff2 = (grad_x * dx);
        switch type
            case {'real', 'real_with_imag'}
                diff2 = real(diff2);
            case 'complex'
                diff2 = diff2;
            otherwise
                error('Invalid type.');
        end
        err(k) = norm(diff1(:) - diff2(:)) / norm(diff1(:));
    end
    fprintf('Max error (%s): %e\n', text, max(err));
