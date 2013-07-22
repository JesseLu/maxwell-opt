%% example1_brutesearch
% Derivative-free optimization of a nanophotonic grating coupler.


function [efficiency, delta, width] = example1_brutesearch(varargin)

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

    % Get the input power
    fprintf('Calculating input power... ');
    P_in = solve_grating([], [], options.flatten);
    fprintf('Input power: %e\n', P_in);

    n = length(options.delta);
    fun = @(x)  -1/P_in * solve_grating(x(1:n), x(n+1:end), options.flatten);
    x0 = [options.delta(:); options.width(:)];
    search_options = optimset(  'Display', 'iter', ...
                                'TolX', 0.05, ...
                                'Tolfun', 1e-16, ...
                                'MaxFunEvals', 1e3, ...
                                'MaxIter', 1e3, ...
                                'FunValCheck', 'on', ...
                                'PlotFcns', {[]}, ...
                                'OutputFcn', {[]});

    % Perform the optimization.
    [x, fval] = fminsearch(fun, x0, search_options);


        %
        % Get the result.
        %

    delta = x(1:n);
    width = x(n+1:end);
    efficiency = -fval;
    
        



end


function [P, E, H, grid, eps] = solve_grating(delta, width, flatten)
% Simulate a grating coupler.

    if isempty(delta) && isempty(width)
        no_struct = true;
    else
        no_struct = false;
    end


        %
        % Create grid.
        %

    % Make a grid for a wavelength of 1550 nm.
    d = 0.05;
    if flatten
        [grid, eps] = maxwell_grid(2*pi/1.55, -4:d:4, 0, -1.5:d:1.5); % 2D.
    else
        [grid, eps] = maxwell_grid(2*pi/1.55, -4:d:4, -3:d:3, -1.5:d:1.5);
    end


        %
        % Setup the structure.
        %

    % Structure constants.
    height = 0.2;
    si_eps = 13;
    air_eps = 1;
    
    % Draw tapered waveguide.
    function [inside] = tapered_wg(x, y, z)
        wg_width = 0.4 + 3.6 * (x < 0) + ...
                        3.6 .* cos(pi*x./4).^2 .* (x >= 0 & x <= 2);
        inside =    x >= -3 & ...
                    abs(y) <= wg_width/2 & ...
                    abs(z) <= height/2;
    end
    eps = maxwell_shape(grid, eps, si_eps, @tapered_wg);

    % Draw half-trenches.
    x_curr = 0;
    for k = 1 : length(delta)
        x_curr = x_curr + abs(delta(k));
        eps = maxwell_shape(grid, eps, air_eps, ...
                    maxwell_box_smooth( [-x_curr 0 height/2], ...
                                        [abs(width(k)) 3.6 height], ...
                                        'smooth_dist', 0.025));
    end

%     if flatten
%         subplot 211; maxwell_view(grid, eps, [], 'y', [nan 0 nan]);
%     else
%         subplot 221; maxwell_view(grid, eps, [], 'y', [nan nan height/4]);
%         subplot 222; maxwell_view(grid, eps, [], 'y', [nan 0 nan]);
%     end


        %
        % Use Gaussian as initial excitation.
        %

    J = maxwell_gaussian(grid, eps, [-1.5 0 0.8], [4 4 -inf], ...
                        'y', 0.8, 2.0);

        
        %
        % Solve.
        %

    if no_struct
        for k = 1 : 3
            eps{k} = air_eps + 0*eps{k};
        end
    end

    if ~flatten; figure(2); end
    [E, H] = maxwell_solve(grid, eps, J);

    figure(1);
    maxwell_view(grid, eps, E, 'y', [nan 0 nan], 'field_phase', nan); % Visualize.


        % 
        % Measure power in reflected wave.
        %

    if no_struct
        P = maxwell_flux(grid, [E H], [0 0 0], [1e9 1e9 -inf]);
    else
        [~, E1, H1] = maxwell_wgmode(grid, eps, [2.5 0 0], [+inf 2 2]);
        P = maxwell_flux(grid, [E H], [E1 H1]);
    end
end
