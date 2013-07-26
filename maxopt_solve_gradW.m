%% maxopt_solve_gradW
% Calculate structural gradients for frequency of the eigenmode.

function [struct_grad] = maxopt_solve_gradient( grid, omega, E, grad_w, params0, ...
                                                create_eps, varargin)
                                                        

        %
        % Validate and parse inputs.
        %

    my_validate_grid(grid, mfilename);
    my_validate_field(E, grid.shape, 'E', mfilename);

    validateattributes(grad_w, {'numeric'}, ...
                        {'nonnan', 'finite', 'scalar'}, ...
                        'grad_w', mfilename);

    validateattributes(create_eps, {'function_handle'}, {}, 'create_eps', mfilename);

    % Optional arguments.
    options = my_parse_options(struct(  'delta_p', 1e-6, ...
                                        'solver_args', {{}}, ...
                                        'fitness', [], ...
                                        'solver', [], ...
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
    y0 = my_left_eigenvector(grid, x0);

    p2z = @(p) vec(create_eps(p));
    z0 = p2z(params0);

    B = spdiags(-omega^2 * vec(E), 0, N, N);
    lambda = omega^2;


        %
        % Find the derivative.
        %

    % Find the dw/dlambda derivative.
    grad_l = grad_w * 0.5 * lambda^(-1/2);
    % my_gradient_test(@(omega) options.fitness(omega), grad_w, omega, true, 'df/dw');
    my_gradient_test(@(lambda) options.fitness(sqrt(lambda)), grad_l, omega^2, true, 'df/dlambda');

    % Find the dlambda/dz derivative.
    grad_z = (1 / (y0' * x0)) * (y0' * B); % Using B is not right here. Need to use "F-field".
    % A = maxwell_axb(grid, unvec(z), unvec(x0), unvec(x0));
    lambda
    my_gradient_test(@(z) (options.solver(unvec(z))), grad_z, z0, false, 'df/dz');

    % Find the dz/dp derivative.
    progress_text = '';
    for k = 1 : numel(params0) 
        p = params0;
        p(k) = p(k) + options.delta_p;
        z = p2z(p);
        grad_p(:, k) = sparse((z - z0) ./ options.delta_p);

        fprintf(repmat('\b', 1, length(progress_text)));
        progress_text = sprintf('[%d/%d gradients computed] ', k, numel(params0));
        fprintf(progress_text);

%         % Debug.
%         maxwell_view(grid, unvec((z - z0) ./ options.delta_p), [], 'y', [nan nan 0], 'clims', [-1 1]);
%         pause
    end

    grad_z = -y' * B; % Form the df/dz derivative.
    df_dp = grad_z * grad_p; % Form the structural gradient.
    grad = real(df_dp).';


    if options.check_gradients
        % Check result of the dagger solve.
        A = maxwell_axb(grid, unvec(z0), E, E);
        fprintf('Error from A_dagger solve: %e\n', norm(A'*y - grad_x0));

        my_gradient_test(p2z, grad_p, params0, false, 'dz/dp') % Test grad_p (dz/dp).

        if ~isempty(options.fitness)
            % Check grad_z.
            my_gradient_test(@(z) options.fitness(unvec(z)), grad_z, z0, true, 'df/dz');

            % Check struct_grad.
            my_gradient_test(@(p) options.fitness(unvec(p2z(p))), df_dp, params0, true, 'df/dp');
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

