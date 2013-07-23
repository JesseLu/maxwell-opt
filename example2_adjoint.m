%% example2_adjoint
% Derivative-based optimization of an H1 photonic crystal resonator.


function [fun, fun1, Emag, params] = example2_adjoint(varargin)

        %
        % Parse inputs.
        %

    options = my_parse_options(struct(  'delta', 0.5 * [0 ones(1, 5)], ... 
                                        'width', 0.1 * ones(1, 6), ...
                                        'flatten', false), ...
                                varargin, mfilename);


        %
        % Set up the optimization problem.
        %

    pc_size = [1 1];
    x0 = zeros(2*prod(pc_size), 1); % Start with no shifts.
    fun = @(x) solve_resonator(pc_size, x, options.flatten, true);
    fun1 = @(x) solve_resonator(pc_size, x, options.flatten, false);

    search_options = optimset(  'Display', 'iter', ...
                                'TolX', 1e-3, ...
                                'Tolfun', 1e-6, ...
                                'DerivativeCheck', 'off', ...
                                'Diagnostics', 'on', ...
                                'LargeScale', 'on', ...
                                'MaxFunEvals', 1e3, ...
                                'MaxIter', 1e3, ...
                                'FunValCheck', 'on', ...
                                'GradObj', 'on', ...
                                'PlotFcns', {[]}, ...
                                'OutputFcn', {[]});
    
    fun(x0);
    return 
    % Perform the optimization.
    [x, fval] = fminunc(fun, x0, search_options);


end

function [fval, df_dp, E, H, grid, eps] = ...
                    solve_resonator(pc_size, shifts, flatten, calc_grad)
% Simulate a photonic crystal resonator.
% Also return the structural gradient.

        %
        % Simulate with a central current source.
        %

    [grid, eps, J] = make_resonator_structure(pc_size, shifts, flatten);

    % Use a central point excitation.
    [x, y, z] = maxwell_pos2ind(grid, 'Ey', [0 0 0]); % Get central Ey component.
    x = x-1; % Slight adjustment.
    y = y-1;
    J{2}(x+[0 1], y, z) = 1;

    % Solve.
    if ~flatten; figure(2); end
    [E, H] = maxwell_solve(grid, eps, J);

    figure(1); maxwell_view(grid, eps, E, 'y', [nan nan 0]); % Visualize.


        % 
        % Measure power reflected back to the center (figure of merit).
        %

    function [fval] = fitness(E)
        E_meas = [E{2}(x, y, z); E{2}(x+1, y, z)];
        fval = -0.5 * norm(E_meas)^2; % This is the figure of merit.
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
    grad_E{2}(x, y, z) = -E{2}(x, y, z);
    grad_E{2}(x+1, y, z) = -E{2}(x+1, y, z);

    % Function handle for creating the structure.
    function [eps] = make_eps(params)
        [~, eps] = make_resonator_structure(pc_size, params, flatten);
    end

    % Calculate the structural gradient.
    df_dp = maxopt_gradients(grid, E, grad_E, shifts, @make_eps, ...
                'fitness', @(eps) fitness(maxwell_solve(grid, eps, J)), ...
                'check_gradients', true);
end



function [grid, eps, J] = make_resonator_structure(pc_size, shifts, flatten)
% Function to create a square lattice photonic crystal structure.

        %
        % Create grid.
        %

    % Make a grid for a wavelength of 1550 nm.
    d = 0.05;
    if flatten
        [grid, eps, ~, J] = maxwell_grid(2*pi/1.55, -3.2:d:3.2, -3.2:d:3.2, 0); % 2D.
    else
        [grid, eps, ~, J] = maxwell_grid(2*pi/1.55, -3.2:d:3.2, -3.2:d:3.2, -1.5:d:1.5);
    end


        %
        % Setup the structure.
        %

    % Structure constants.
    height = 0.3;
    radius = 0.15;
    a = 0.5;
    si_eps = 13;
    air_eps = 1;

    % Draw slab.
    eps = maxwell_shape(grid, eps, si_eps, ...
                        maxwell_box([0 0 0], [5.4 5.4 height]));

    % Draw photonic crystal.
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

    for k = 1 : length(pos)
        eps = maxwell_shape(grid, eps, air_eps, ...
                            maxwell_cyl_smooth(pos{k}, radius, 2*height, ...
                                                'smooth_dist', d));
    end
end
