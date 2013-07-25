%% case3_L3
% L3 photonic crystal resonator.

function [fun, x0] = case3_L3(type, varargin)

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
        % Return appropriate function handle.
        %
        
    function [omega, E, H, grid, eps] = get_fields(varargin)
        [~, ~, omega, E, H, grid, eps] = solve_resonator(varargin{:});
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


function [fval, df_dp, omega, E, H, grid, eps] = ...
                    solve_resonator(pc_size, shifts, flatten, calc_grad)
% Simulate a photonic crystal resonator.
% Also return the structural gradient.

        %
        % Sanity check for shifts.
        %

    if any(abs(shifts) > 3.2) % Check for shifts which are way too large.
        fval = 1e9;
        df_dp = -1e9 * shifts;
        return
    end


        %
        % Get the initial excitation.
        %

    [grid, eps, J] = make_resonator_structure(pc_size, shifts, flatten);

    if flatten

        % Use a central point excitation.
        [x, y, z] = maxwell_pos2ind(grid, 'Ey', [0 0 0]); % Get central Ey component.
        x = x-1; % Slight adjustment.
        y = y;
        J{2}(x+[0 1], y, z) = 1;
    
    else % If 3D recursively use 2D eigenmode as initial excitation.

        [~, ~, ~, E] = solve_resonator(pc_size, shifts, true, calc_grad);
        [x, y, z] = maxwell_pos2ind(grid, 'Ey', [0 0 0]); % Get central Ey component.
        J{2}(:,:,z) = E{2};
    end

    [E, H] = maxwell_solve(grid, eps, J);
    [omega, E, H] = maxwell_solve_eigenmode(grid, eps, E);

    figure(1); maxwell_view(grid, eps, E, 'y', [nan nan 0], 'field_phase', nan); % Visualize.


        % 
        % Measure power reflected back to the center (figure of merit).
        %

    function [fval] = fitness(E)
        E_meas = [E{2}(x, y, z); E{2}(x+1, y, z)];
        fval = -0.5 * norm(E_meas)^2; % This is the figure of merit.
        fval = -norm(E_meas); % This is the figure of merit.
    end
        
%     E_meas = [E{2}(x, y, z); E{2}(x+1, y, z)];
%     fval = -0.5 * norm(E_meas)^2; % This is the figure of merit.
    fval = fitness(E);


        % 
        % Calculate gradient.
        %

    if ~calc_grad % Skip if not needed.
        df_dp = nan;
        return
    end

    % Field gradient.
    grad_E = my_default_field(grid.shape, 0); 
    a = norm([E{2}(x,y,z); E{2}(x+1,y,z)]);
    grad_E{2}(x, y, z) = -E{2}(x, y, z) / a;
    grad_E{2}(x+1, y, z) = -E{2}(x+1, y, z) / a;

    % Function handle for creating the structure.
    function [eps] = make_eps(params)
        [~, eps] = make_resonator_structure(pc_size, params, flatten);
    end

    % Calculate the structural gradient.
    if ~flatten; figure(3); end
    df_dp = maxopt_gradients(grid, E, grad_E, shifts, @make_eps, ...
                'fitness', @(eps) fitness(maxwell_solve(grid, eps, J)), ...
                'check_gradients', false);
end



function [grid, eps, J] = make_resonator_structure(pc_size, shifts, flatten)
% Function to create a square lattice photonic crystal structure.

        %
        % Create grid.
        %

    % Make a grid for a wavelength of 1550 nm.
    d = 0.025;
    if flatten
        [grid, eps, ~, J] = maxwell_grid(2*pi/1.911, -2.6:d:2.6, -2:d:2, 0); % 2D.
    else
        [grid, eps, ~, J] = maxwell_grid(2*pi/1.55, -2.6:d:2.6, -2:d:2, -1.5:d:1.5);
    end


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

    % Draw photonic crystal.
    shifts = reshape(shifts, [round(numel(shifts)/2) 2]);
    pos = {};
    cnt = 1;
    for i = 0 : pc_size(1)-1
        for j = 0 : pc_size(2)-1
            p{1} = a * ((i+floor(j/2))*a1 + j*a2) + ...
                    [shifts(cnt, 1) shifts(cnt, 2) 0];

            if p{1}(1) < 0 % Skip.
                continue
            end

            if (i == 0 || i == 1) && j == 0 % Skip the 3 holes to be removed.
                continue
            end

            if p{1}(1) == 0
                p{2} = [];
            else
                p{2} = p{1} .* [-1 1 1];
            end

            if p{1}(2) == 0
                p{3} = [];
            else
                p{3} = p{1} .* [1 -1 1];
            end

            if p{1}(1) == 0 || p{1}(2) == 0
                p{4} = [];
            else
                p{4} = p{1} .* [-1 -1 1];
            end

            pos = [pos, p];
            cnt = cnt + 1;
        end
    end

    for k = 1 : length(pos)
        if ~isempty(pos{k})
            eps = maxwell_shape(grid, eps, air_eps, ...
                                maxwell_cyl_smooth(pos{k}, radius, 2*height, ...
                                                    'smooth_dist', d));
        end
    end
end
