%% maxopt_case_grating
% Sets up a grating coupler optimization problem.

function [fun, x0] = maxopt_case_grating(type, varargin)

        %
        % Parse inputs.
        %

    validateattributes(type, {'char'}, {'vector'}, 'type', mfilename);

    options = my_parse_options(struct(  'flatten', false), ...
                                varargin, mfilename);


        %
        % Return recommended starting parameters.
        %

    delta0 = 0.5 * [0 ones(1, 5)]; % Spacings.
    width0 = 0.1 * ones(1, 6); % Widths.
    x0 = [delta0(:); width0(:)];


        %
        % Calculate input power.
        % 

    P_in = abs(solve_grating([], options.flatten, 1, false));

        %
        % Return appropriate function handle.
        %
        
    function [E, H, grid, eps] = get_fields(varargin)
        [~, ~, E, H, grid, eps] = solve_grating(varargin{:});
    end

    switch type
        case 'get_fields'
            fun = @(x) get_fields(x, options.flatten, P_in, false);
        case 'fval'
            fun = @(x) solve_grating(x, options.flatten, P_in, false);
        case 'grad_f'
            fun = @(x) solve_grating(x, options.flatten, P_in, true);
        otherwise
            error('Invalid type.');
    end
end

function [fval, grad_f, E, H, grid, eps] = solve_grating(params, flatten, P_in, calc_grad)
% Simulate a grating coupler.

    n_2 = round(length(params)/2);
    delta = params(1:n_2);
    width = params(n_2+1:end);

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
    for k = 1 : n_2
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

    fval = -P / P_in;
    grad_f = [];
end
