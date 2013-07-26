%% maxopt_solve_gradient
% Calculate structural gradients

function [struct_grad, omega_grad] = ...
            maxopt_solve_gradient(  grid, E, grad_E, params0, ...
                                    create_eps, varargin)
                                                        

        %
        % Validate and parse inputs.
        %

    my_validate_grid(grid, mfilename);
    my_validate_field(E, grid.shape, 'E', mfilename);
    my_validate_field(grad_E, grid.shape, 'grad_E', mfilename);

    validateattributes(create_eps, {'function_handle'}, {}, 'create_eps', mfilename);

    % Optional arguments.
    options = my_parse_options(struct(  'delta_p', 1e-6, ...
                                        'solver_args', {{}}, ...
                                        'fitness', [], ...
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
%     solve_A_dagger = @(z, b) maxopt_solve_adjoint(grid, unvec(z), ...
%                                             unvec(b ./ (-i * grid.omega)), ...
%                                             options.solver_args{:});
% 
    solve_A_dagger = @(z, b) maxopt_solve_adjoint(grid, unvec(z), ...
                                            unvec(b), ...
                                            options.solver_args{:});

    B = spdiags(-grid.omega^2 * vec(E), 0, N, N);


        %
        % Find the derivative.
        %

    % Initiate solve.
    fprintf('[start adjoint solve] ');
    cb = solve_A_dagger(z0, grad_x0);

%     % Find the dz/dp derivative.
%     progress_text = '';
%     for k = 1 : numel(params0) 
%         p = params0;
%         p(k) = p(k) + options.delta_p;
%         z = p2z(p);
%         grad_p(:, k) = sparse((z - z0) ./ options.delta_p);
% 
%         fprintf(repmat('\b', 1, length(progress_text)));
%         progress_text = sprintf('[%d/%d gradients computed] ', k, numel(params0));
%         fprintf(progress_text);
%     end

    % Obtain result from dagger solve.
    while ~cb(); end
    [~, y] = cb();
    y = vec(y);

        [A, x, b] = maxwell_axb(grid, unvec(z0), unvec(y), unvec(grad_x0));
        norm(A'*x-b) / norm(b)
        fprintf('Error from A_dagger solve: %e\n', norm(A'*y - grad_x0));
    grad_z = -y' * B; % Form the df/dz derivative.
    df_dp = grad_z * grad_p; % Form the structural gradient.
    struct_grad = real(df_dp).';
    omega_grad = nan;

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

%         % Check equivalence of Ax-b and Bz-d.
%         z0 = p2z(params0);
%         x0 = vec(E);
%         multA = maxwell_axb(grid, unvec(z0), E, E, 'functional', true);
%         b = randn(N, 1);
%         d = b - (multA(x0) + grid.omega^2 * (z0 .* x0));
% 
%         res1 = multA(x0) - b;
%         res2 = B * z0 - d;
% 
%         fprintf('Error between Ax-b and Bz-d: %e\n', norm(res1 - res2)/norm(res1));
    end

