%% maxopt_field_gradient
% Calculate structural gradients for the E-field of a simulation.

function [param_grad, eps_grad] = maxopt_field_gradient(grid, E, fitness_fun, ...
                                                        params0, create_eps, ...
                                                        varargin)
                                                        

        %
        % Validate and parse inputs.
        %

    my_validate_grid(grid, mfilename);
    my_validate_field(E, grid.shape, 'E', mfilename);
    validateattributes(params0, {'numeric'}, {'real', 'nonnan'}, 'params0', mfilename);
    validateattributes(create_eps, {'function_handle'}, {}, 'create_eps', mfilename);

    % Check fitness_fun.
    validateattributes(fitness_fun, {'function_handle'}, {}, 'fitness_fun', mfilename);
    [fval, grad_E] = fitness_fun(E);
    validateattributes(fval, {'numeric'}, {'scalar', 'real'}, ...
                        'fval (from fitness_fun)', mfilename);
    my_validate_field(grad_E, grid.shape, 'grad_Efval (from fitness_fun)', mfilename);

    % Check fitness_fun's gradient.
    [vec, unvec] = my_vec(grid.shape); % Helper function.
    err = my_gradient_test(@(x) fitness_fun(unvec(x)), vec(grad_E), vec(E), 'real_with_imag', '');
    if err > 1e-3
        warning('Error in fitness_fun gradient is large (%e).', err);
    end


    % Optional arguments.
    options = my_parse_options(struct(  'delta_p', 1e-6, ...
                                        'solver_args', {{}}, ...
                                        'solver_fun', [], ...
                                        'check_gradients', false), ...
                                varargin, mfilename);

    validateattributes(options.delta_p, {'numeric'}, ...
                        {'nonnan', 'finite', 'real'}, ...
                        'delta_p', mfilename);

    validateattributes(options.solver_args, {'cell'}, {}, ...
                        'solver_args', mfilename);

    validateattributes(options.check_gradients, {'logical'}, {'scalar'}, ...
                        'check_gradients', mfilename);


        %
        % Convert to linear algebra language.
        %

    [vec, unvec] = my_vec(grid.shape); % Helper functions.
    N = 3 * prod(grid.shape);

    params0 = params0(:);

    x0 = vec(E);
    grad_x0 = vec(grad_E); 

    p2z = @(p) vec(create_eps(p));
    z0 = p2z(params0);

    solve_A_dagger = @(z, b) maxopt_solve_adjoint(grid, unvec(z), ...
                                            unvec(b ./ (-i * grid.omega)), ...
                                            options.solver_args{:});

    B = spdiags(-grid.omega^2 * vec(E), 0, N, N);


        %
        % Compute the adjoint and the dz/dp derivative while waiting.
        %

    % Initiate adjoint solve.
    fprintf('[start adjoint solve] '); 
    cb = solve_A_dagger(z0, grad_x0);

    % Find the dz/dp derivative.
    dz_dp = my_parameter_gradient(p2z, params0, options.delta_p); 

    % Complete adjoint solve.
    while ~cb(); end 
    [~, y] = cb();
    y = vec(y);


        %
        % Form the derivative.
        %
    
    % Form the df/dz derivative.
    df_dz = -y' * B; 

    % Form df/dp the parameter derivative.
    df_dp = df_dz * dz_dp; 

    % Output parameters.
    param_grad = real(df_dp');
    eps_grad = df_dz';


        %
        % Check gradients, if desired.
        %

    if options.check_gradients
        % Check result of the dagger solve.
        A = maxwell_axb(grid, unvec(z0), E, E);
        fprintf('Error from A_dagger solve: %e\n', norm(A'*y - grad_x0));

        my_gradient_test(p2z, dz_dp', params0, 'real', 'dz/dp'); % Test dz/dp.

        if ~isempty(options.solver_fun)
            my_gradient_test(@(z) fitness_fun(options.solver_fun(unvec(z))), df_dz', z0, 'real', 'df/dz');

            my_gradient_test(@(p) fitness_fun(options.solver_fun(unvec(p2z(p)))), df_dp', params0, 'real', 'df/dp');
        end

%             % Check equivalence of Ax-b and Bz-d.
%             z0 = p2z(params0);
%             x0 = vec(E);
%             multA = maxwell_axb(grid, unvec(z0), E, E, 'functional', true);
%             b = randn(N, 1);
%             d = b - (multA(x0) + grid.omega^2 * (z0 .* x0));
%     
%             res1 = multA(x0) - b;
%             res2 = B * z0 - d;
%     
%             fprintf('Error between Ax-b and Bz-d: %e\n', norm(res1 - res2)/norm(res1));
    end
