%% my_grad_descent
% Simple gradient descent algorithm 

function [x, fval, hist] = my_grad_descent(fun, x0, varargin)

    
    % Parse options.
    options = my_parse_options(struct(  'init_step', 1, ...
                                        'max_delta', 1, ...
                                        'step_shrink', 0.5, ...
                                        'step_grow', 1.1, ...
                                        'max_iters', 100, ...
        'vis_progress', @(hist) fprintf('%d: %e\n', length(hist), hist(end))), ...
                                varargin, mfilename);

    x = x0;
    [f, dx] = fun(x);
    step_size = options.init_step;

    first_success_step = true;

    for k = 1 : options.max_iters
        hist(k) = f;
        options.vis_progress(hist);

        if any(abs(dx) > options.max_delta)
            dx = dx * (options.max_delta / max(abs(dx)));
        end

        x1 = x - step_size * dx;
        [f1, dx] = fun(x1);

        if f1 < f % Grow.
            step_size = step_size * options.step_grow;
        else % Shrink.
            step_size = step_size * options.step_shrink;
        end

        f = f1;
        x = x1;


    end
    fval = f;

        



