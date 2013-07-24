%% example3_cavity
% Derivative-based optimization of an L3 cavity mode.


function [params, E, H, grid, eps] = example3_cavity(varargin)

        %
        % Parse inputs.
        %

    options = my_parse_options(struct(  'delta', 0.5 * [0 ones(1, 5)], ... 
                                        'width', 0.1 * ones(1, 6), ...
                                        'sim_only', false, ...
                                        'flatten', false), ...
                                varargin, mfilename);


        %
        % Set up the optimization problem.
        %

    pc_size = [6 5];
    x = zeros(2*prod(pc_size), 1); % Start with no shifts.
    fun = @(x) solve_resonator(pc_size, x, options.flatten, true);

    function vis_progress(hist)
        figure(2);
        plot(hist, '.-');
        xlabel('optimization iterations');
        ylabel('fval');
        title('structure optimization progress');
    end

    if ~options.sim_only
        search_options = optimset(  'Display', 'iter', ...
                                    'TolX', 1e-6, ...
                                    'Tolfun', 1e-6, ...
                                    'DerivativeCheck', 'off', ...
                                    'Diagnostics', 'on', ...
                                    'LargeScale', 'on', ...
                                    'MaxFunEvals', 1e3, ...
                                    'MaxIter', 1e3, ...
                                    'FunValCheck', 'on', ...
                                    'GradObj', 'off', ...
                                    'TypicalX', 1e-3 + ones(size(x)), ...
                                    'PlotFcns', {[]}, ...
                                    'OutputFcn', {[]});
        
        
        [x, fval, hist] = my_grad_descent(fun, x,  'init_step', 0.1, ...
                                                    'max_delta', 0.1, ...
                                                    'vis_progress', @vis_progress);
    end

    [~, ~, E, H, grid, eps] = solve_resonator(pc_size, x, options.flatten, false);
    params = reshape(x, [round(length(x)/2) 2]);

end

function [fval, df_dp, E, H, grid, eps] = ...
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
        % Simulate with a central current source.
        %

    [grid, eps, J] = make_resonator_structure(pc_size, shifts, flatten);

    maxwell_view(grid, eps, [], 'y', [nan nan 0]);
    pause
    % Use a central point excitation.
    [x, y, z] = maxwell_pos2ind(grid, 'Ey', [0 0 0]); % Get central Ey component.
    x = x-1; % Slight adjustment.
    y = y;
    J{2}(x+[0 1], y, z) = 1;

    % Solve.
    if ~flatten; figure(3); end
    [E, H] = maxwell_solve(grid, eps, J);

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
        [grid, eps, ~, J] = maxwell_grid(2*pi/1.55, -2.6:d:2.6, -2:d:2, 0); % 2D.
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
                                                    'smooth_dist', 3*d));
        end
    end
end
