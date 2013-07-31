%% maxopt_case_L3
% Sets up the optimization for an L3 photonic crystal resonator.

function [fun, x0] = maxopt_case_L3(type, varargin)

        %
        % Parse inputs.
        %

    validateattributes(type, {'char'}, {'vector'}, 'type', mfilename);

    options = my_parse_options(struct(  'flatten', false), ...
                                varargin, mfilename);


        %
        % Return recommended starting parameters.
        %

    pc_size = [6 5];
    x0 = zeros(2*prod(pc_size), 1); % Start with no shifts.


        %
        % Helper functions for solving for the eigenmode.
        %

    % Save previously used values.
    omega_cache = []; 
    E_cache = [];

    function [fval, grad_f, omega, E, H, grid, eps] = cached_solve(varargin)
    % This is a cached solve which uses the most recent omega and E values
    % as initial guesses for the eigenmode solve.

        [fval, grad_f, omega, E, H, grid, eps] = ...
            solve_resonator(varargin{:}, omega_cache, E_cache); % Solve.

        % Update cached values.
        omega_cache = omega;
        E_cache = E;
    end

    function [omega, E, H, grid, eps] = get_fields(varargin)
    % Just return the resulting eigenmode field.
        [~, ~, omega, E, H, grid, eps] = cached_solve(varargin{:});
    end


        %
        % Return appropriate function handle.
        %

    switch type
        case 'get_fields'
            fun = @(x) get_fields(pc_size, x, options.flatten, false);
        case 'fval'
            fun = @(x) cached_solve(pc_size, x, options.flatten, false);
        case 'grad_f'
            fun = @(x) cached_solve(pc_size, x, options.flatten, true);
        otherwise
            error('Invalid type.');
    end
end



% Simulate an L3 photonic crystal resonator and return a gradient.
function [fval, grad_f, omega, E, H, grid, eps] = ...
                    solve_resonator(pc_size, shifts, flatten, calc_grad, omega_guess, E_guess)

        %
        % Sanity check for shifts.
        %

    if any(abs(shifts) > 3.2) % Check for shifts which are way too large.
        fval = 1e9;
        grad_f = -1e9 * shifts;
        return
    end


        %
        % Get the initial guess field.
        % If 2D, use a central excitation at the guess frequency.
        % If 3D, use the 2D eigenmode as as the excitation at the central plane
        % to obtain the initial guess E-field.
        %

    [grid, eps, J] = make_resonator_structure(pc_size, shifts, flatten);
    [x, y, z] = maxwell_pos2ind(grid, 'Ey', [0 0 0]); % Get central Ey component.

    if isempty(omega_guess) % Don't need if we already have a cached guess.

        % Obtain the excitation that we'll use for obtaining the guess field.
        if flatten % 2D case.
            x = x-1; % Slight adjustment.
            y = y;
            J{2}(x+[0 1], y, z) = 1; % Use a central point excitation.
        
        else % 3D case.
            [~, ~, ~, E] = solve_resonator(pc_size, shifts, true, false, [], []); % 2D solve
            J{1}(:,:,z) = E{1}; % Use it as the excitation at the central plane.
            J{2}(:,:,z) = E{2}; 
        end

        % Run a simulation to obtain the initial guess field.
        [E, H] = maxwell_solve(grid, eps, J);

    else % We have a cached guess, no need to simulate to get one.
        [grid, eps, J] = make_resonator_structure(pc_size, shifts, flatten);
        E = E_guess;
    end


        %
        % Solve for the eigenmode.
        %

    [omega, E, H] = maxwell_solve_eigenmode(grid, eps, E, 'eig_max_iters', 3);
    figure(1); maxwell_view(grid, eps, E, 'y', [nan nan 0]); % Visualize.


        % 
        % Compose and evaluate fitness function (function to minimize).
        %

    [vec, unvec] = my_vec(grid.shape);
    function [fval, grad_w] = fitness(w)
    % Calculates figure of merit (fitness function) and its derivative.
        fval = 0.5 * norm(w)^2;
        grad_w = w;
        fval = imag(w);
        grad_w = 1i;
    end
        
    [fval, grad_w] = fitness(omega);


        % 
        % Calculate structural gradient needed for gradient descent optimization.
        %

    if ~calc_grad % Skip if not needed.
        grad_f = nan;
        return
    end

    function [eps] = make_eps(params)
    % Function handle for creating the structure.
        [~, eps] = make_resonator_structure(pc_size, params, flatten);
    end

    function [lambda] = solver(eps)
    % Function that evaluates the fitness based on eps.
    % Only used for gradient checking.
        [omega_fit, E_fit, H_fit] = maxwell_solve_eigenmode(grid, eps, E, 'err_thresh', 1e-2);
        lambda = omega_fit^2;
    end

    % Calculate the structural gradient.
    if ~flatten; figure(3); end
    grad_f = maxopt_freq_gradient(grid, E, omega, @fitness, shifts, @make_eps, ...
                'solver', @solver, ...
                'check_gradients', false);
end



function [grid, eps, J] = make_resonator_structure(pc_size, shifts, flatten, varargin)
% Function to create a square lattice photonic crystal structure.

        %
        % Create grid.
        % Note that we keep frequencies constant in order to have constant PML
        % scale factors.
        %

    % Make a grid for a wavelength of 1550 nm.
    d = 0.025;
    x = -2.6:d:2.6;
    y = -2:d:2;
    z = -1.5:d:1.5;
    omega = 2*pi/1.583;

    if flatten
        omega = 2*pi/1.911;
        z = 0;
    end

    [grid, eps, ~, J] = maxwell_grid(omega, x, y, z);


        %
        % Setup the structure.
        %

    % Structure constants.
    height = 0.24;
    radius = 0.12;
    a = 0.4;
    a1 = [1 0 0];
    a2 = [-0.5 sqrt(3)/2 0];
    si_eps = 13;
    air_eps = 1;

    % Draw slab.
    eps = maxwell_shape(grid, eps, si_eps, ...
                        maxwell_box([0 0 0], [inf inf height]));

    % Determine the positions of the etched holes.
    % Things get a little complicated here because we force x- and y-symmetry.
    % In particular, note that many shift parameters go ignored.
    shifts = reshape(shifts, [round(numel(shifts)/2) 2]);
    pos = {};
    cnt = 0;
    for i = 0 : pc_size(1)-1
        for j = 0 : pc_size(2)-1
            cnt = cnt + 1;
            p{1} = a * ((i+floor(j/2))*a1 + j*a2);

            % Skip holes that have negative x-position.
            if p{1}(1) < 0 
                continue
            end

            % Skip the 3 holes to be removed.
            if (i == 0 || i == 1) && j == 0 
                continue
            end
            
            % Add x-shift.
            if p{1}(1) ~= 0
                p{1}(1) = p{1}(1) + shifts(cnt, 1);
            end

            % Add y-shift.
            if p{1}(2) ~= 0
                p{1}(2) = p{1}(2) + shifts(cnt, 2);
            end

            % Reflect across x.
            if p{1}(1) == 0
                p{2} = [];
            else
                p{2} = p{1} .* [-1 1 1];
            end

            % Reflect across y.
            if p{1}(2) == 0
                p{3} = [];
            else
                p{3} = p{1} .* [1 -1 1];
            end

            % Reflect across x and y.
            if p{1}(1) == 0 || p{1}(2) == 0
                p{4} = [];
            else
                p{4} = p{1} .* [-1 -1 1];
            end

            pos = [pos, p]; % Add to position list.
        end
    end

    % Create the holes in the L3 structure.
    for k = 1 : length(pos)
        if ~isempty(pos{k})
            eps = maxwell_shape(grid, eps, air_eps, ...
                                maxwell_cyl_smooth(pos{k}, radius, 2*height, ...
                                                    'smooth_dist', 2*d));
        end
    end
end
