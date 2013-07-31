function [max_err, err] = my_gradient_test(fun, grad_x, x0, type, text, varargin)
% Just checks gradients.
   
    f0 = fun(x0);

    delta = mean(abs(x0(:))) / 1e6;
    if delta == 0
        delta = 1e-6;
    end

    for k = 1 : 10
        dx = delta * randn(size(x0));
        if strcmp(type, 'real_with_imag')
            dx = dx + 1i * delta * randn(size(x0));
        end
        f1 = fun(x0+dx);
        diff1 = f1 - f0;
        diff2 = (grad_x' * dx);
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

    max_err = max(err);

    if length(text) > 0
        fprintf('Max error (%s): %e\n', text, max_err);
    end

