%% maxopt_freq_gradient
% Calculate structural gradients for frequency of the eigenmode.

function [param_grad, eps_grad] = maxopt_freq_gradient(grid, E, omega, fitness_fun, ...
                                                        params0, create_eps, ...
                                                        varargin)
                                                        
        %
        % Validate and parse inputs.
        %

    my_validate_grid(grid, mfilename);
    my_validate_field(E, grid.shape, 'E', mfilename);
    validateattributes(omega, {'numeric'}, {'scalar', 'nonnan', 'finite'}, 'omega', mfilename);
    validateattributes(params0, {'numeric'}, {'real', 'nonnan'}, 'params0', mfilename);

    % Check fitness_fun.
    validateattributes(fitness_fun, {'function_handle'}, {}, 'fitness_fun', mfilename);
    [fval, grad_omega] = fitness_fun(omega);
    validateattributes(fval, {'numeric'}, {'scalar', 'real'}, ...
                        'fval (from fitness_fun)', mfilename);
    validateattributes(grad_omega, {'numeric'}, {'scalar', 'nonnan', 'finite'}, ...
                        'grad_omega (from fitness_fun)', mfilename);

    % Check fitness_fun's gradient.
    err = my_gradient_test(fitness_fun, grad_omega, omega, 'real_with_imag', '');
    if err > 1e-3
        warning('Error in fitness_fun gradient is large (%e).', err);
    end

    validateattributes(grad_omega, {'numeric'}, ...
                        {'nonnan', 'finite', 'scalar'}, ...
                        'grad_omega', mfilename);

    validateattributes(create_eps, {'function_handle'}, {}, 'create_eps', mfilename);

    % Optional arguments.
    options = my_parse_options(struct(  'delta_p', 1e-6, ...
                                        'solver', [], ...
                                        'check_gradients', false), ...
                                varargin, mfilename);

    validateattributes(options.delta_p, {'numeric'}, ...
                        {'nonnan', 'finite', 'real'}, ...
                        'delta_p', mfilename);

    validateattributes(options.check_gradients, {'logical'}, {'scalar'}, ...
                        'check_gradients', mfilename);


        %
        % Convert to linear algebra language.
        %

    [vec, unvec] = my_vec(grid.shape); % Helper functions.
    N = 3 * prod(grid.shape);

    params0 = params0(:);

    x0 = vec(E);
    y0 = my_left_eigenvector(grid, x0);

    p2z = @(p) vec(create_eps(p));
    z0 = p2z(params0);

    lambda = omega^2;
    df_dw = grad_omega';


        %
        % Find the derivative.
        %

    % Find the df/dlambda derivative.
    lambda2w = @sqrt;
    df_dlambda = df_dw * 0.5 * lambda^(-1/2);

    % Find the dlambda/dz derivative.
    dlambda_dz = -(lambda / (y0' * (z0 .* x0))) * (conj(y0) .* x0).';

    % Find the dz/dp derivative.
    dz_dp = my_parameter_gradient(p2z, params0, options.delta_p); % Find the dz/dp derivative.

    % Put everything together the df/dp derivative that we're looking for.
    df_dz = df_dlambda * dlambda_dz; % Get the df/dz derivative.
    df_dp = df_dz * dz_dp; % Final df/dp derivative.

    % Output parameters (transfer everything to gradients).
    param_grad = real(df_dp).';
    eps_grad = df_dz';


        %
        % Check the derivatives, if desired..
        %

    if options.check_gradients % Check results.
        my_gradient_test(@(lambda) fitness_fun(sqrt(lambda)), df_dlambda, lambda, 'real', 'df/dlambda');
        my_gradient_test(@(z) (options.solver(unvec(z))), dlambda_dz', z0, 'complex', 'dlambda/dz');
        my_gradient_test(@(z) fitness_fun(lambda2w(options.solver(unvec(z)))), df_dz', z0, 'real', 'df/dz');
        my_gradient_test(p2z, dz_dp', params0, 'complex', 'dz/dp'); % Test grad_p (dz/dp).
        my_gradient_test(@(p) fitness_fun(lambda2w(options.solver(unvec(p2z(p))))), df_dp', params0, 'real', 'df/dp');
    end

