%% maxopt_gradient_descent
% Simple gradient descent optimization algorithm.

%%% Syntax
%
% * |maxopt_gradient_descent(fun, x0)|
%   attempts to minimize the function |fun| via gradient-descent
%   starting at the point |x0|.
%   |fun| is a function handle that returns both 
%   the value of the fitness function (figure of merit) and
%   the the gradient of said function with respect to |x|.
%   
% * |[x, fval, hist] = maxopt_gradient_descent(...)|
%   outputs |x|, the final 
%

%%% Description
% |maxopt_gradient_descent| is a simple function that performs a simple 
% gradient-descent optmization algorithm.
% |maxopt_gradient_descent| keeps a step even if it increases 
% (instead of decreasing) the figure of merit.
% However, in this case, the step-length will be decreased.
%

function [x, fval, hist] = maxopt_gradient_descent(fun, x0, varargin)

        %
        % Validate and parse inputs.
        %

    validateattributes(fun, {'function_handle'}, {}, 'fun', mfilename);

    validateattributes(x0, {'numeric'}, ...
                        {'vector', 'real', 'finite', 'nonnan'}, ...
                        'x0', mfilename);
        
    % Optional parameters.
    options = my_parse_options(struct(  'init_step', 1, ...
                                        'max_delta', 1, ...
                                        'step_shrink', 0.5, ...
                                        'step_grow', 1.1, ...
                                        'max_iters', 100, ...
        'vis_progress', @(hist) fprintf('%d: %e\n', length(hist), hist(end))), ...
                                varargin, mfilename);

    simple_check = @(var, var_name) ...
        validateattributes(var, {'numeric'}, ...
                        {'positive', 'scalar', 'finite', 'nonnan'}, ...
                        var_name, mfilename);

    simple_check(options.init_step, 'init_step');
    simple_check(options.max_delta, 'max_delta');
    simple_check(options.step_shrink, 'step_shrink');
    simple_check(options.step_grow, 'step_grow');
    simple_check(options.max_iters, 'max_iters');

    validateattributes(options.vis_progress, {'function_handle'}, {}, ...
                        'vis_progress', mfilename);


        %
        % Perform minimization.
        %

    x = x0;
    [f, dx] = fun(x);
    step_size = options.init_step;

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

        



