%% maxopt_case_squarepc
% Sets up optimization problem for a square photonic crystal array.

function [fun, x0] = maxopt_case_squarepc(type, varargin)

        %
        % Parse inputs.
        %

    validateattributes(type, {'char'}, {'vector'}, 'type', mfilename);

    options = my_parse_options(struct(  'flatten', false), ...
                                varargin, mfilename);


        %
        % Return recommended starting parameters.
        %

    pc_size = [3 3]; % Actually turns out to be a symmetric 6x6 crystal.
    x0 = zeros(2*prod(pc_size), 1); % Start with no shifts.


        %
        % Return appropriate function handle.
        %
        
    function [E, H, grid, eps] = get_fields(varargin)
        [~, ~, E, H, grid, eps] = solve_resonator(varargin{:});
    end

    switch type
        case 'get_fields'
            fun = @(x) get_fields(pc_size, x, options.flatten, false);
        case 'fval'
            fun = @(x) solve_resonator(pc_size, x, options.flatten, false);
        case 'grad_f'
            fun = @(x) solve_resonator(pc_size, x, options.flatten, true);
        otherwise
            error('Invalid type.');
    end
end

function [fval, grad_f, E, H, grid, eps] = ...
                    solve_resonator(pc_size, shifts, flatten, calc_grad)
% Simulate a square photonic crystal array with shifts, 
% also return the structural gradient.

        %
        % Sanity check for shifts.
        %

    if any(abs(shifts) > 3.2) % Check for shifts which are way too large.
        fval = 1e9;
        grad_f = -1e9 * shifts;
        return
    end


        %
        % Simulate with a central current source.
        %

    [grid, eps, J] = make_resonator_structure(pc_size, shifts, flatten);

    % Use a central point excitation.
    [x, y, z] = maxwell_pos2ind(grid, 'Ey', [0 0 0]); % Get central Ey component.
    x = x-1; % Slight adjustment.
    y = y;
    J{2}(x+[0 1], y, z) = 1;

    if ~flatten; figure(3); end
    [E, H] = maxwell_solve(grid, eps, J); % Solve.

    figure(1); maxwell_view(grid, eps, E, 'y', [nan nan 0], 'field_phase', nan); % Visualize.


        % 
        % Measure power reflected back to the center (figure of merit).
        %

    function [fval, grad_E] = fitness(E)
    % Calculates figure of merit and its derivative.
        % Figure of merit.
        E_meas = [E{2}(x, y, z); E{2}(x+1, y, z)];
        fval = -sum(abs(E_meas)); % This is the figure of merit.

        % Field gradient.
        grad_E = my_default_field(grid.shape, 0); 
        grad_E{2}(x, y, z) = -E{2}(x, y, z) / abs(E{2}(x, y, z));
        grad_E{2}(x+1, y, z) = -E{2}(x+1, y, z) / abs(E{2}(x+1, y, z));
    end
        
    [fval, grad_E] = fitness(E);


        % 
        % Calculate structural gradient.
        %

    if ~calc_grad % Skip if not needed.
        grad_f = nan;
        return
    end

    function [eps] = make_eps(params)
    % Function handle for creating the structure.
        [~, eps] = make_resonator_structure(pc_size, params, flatten);
    end

    % Calculate the structural gradient.
    if ~flatten; figure(3); end
    grad_f = maxopt_field_gradient(grid, E, @fitness, shifts, @make_eps, ...
                'solver_fun', @(eps) maxwell_solve(grid, eps, J), ...
                'check_gradients', false);
end



function [grid, eps, J] = make_resonator_structure(pc_size, shifts, flatten)
% Function to create a square lattice photonic crystal structure.

        %
        % Create grid.
        %

    % Make a grid for a wavelength of 1550 nm.
    d = 0.05;
    if flatten
        [grid, eps, ~, J] = maxwell_grid(2*pi/1.55, -2:d:2, -2:d:2, 0); % 2D.
    else
        [grid, eps, ~, J] = maxwell_grid(2*pi/1.55, -2:d:2, -2:d:2, -1.5:d:1.5);
    end


        %
        % Setup the structure.
        %

    % Structure constants.
    height = 0.4;
    radius = 0.15;
    a = 0.5;
    si_eps = 13;
    air_eps = 1;

    % Draw slab.
    eps = maxwell_shape(grid, eps, si_eps, ...
                        maxwell_box([0 0 0], [inf inf height]));

    % Determine positions of holes.
    shifts = reshape(shifts, [round(numel(shifts)/2) 2]);
    pos = {};
    cnt = 1;
    for i = 1 : pc_size(1)
        for j = 1 : pc_size(2)
            p{1} = a * [(i-0.5) (j-0.5) 0] + ...
                    [shifts(cnt, 1) shifts(cnt, 2) 0];
            p{2} = p{1} .* [-1 1 1];
            p{3} = p{1} .* [1 -1 1];
            p{4} = p{1} .* [-1 -1 1];
            pos = [pos, p];
            cnt = cnt + 1;
        end
    end

    % Draw photonic crystal.
    for k = 1 : length(pos)
        eps = maxwell_shape(grid, eps, air_eps, ...
                            maxwell_cyl_smooth(pos{k}, radius, 2*height, ...
                                                'smooth_dist', d));
    end
end
