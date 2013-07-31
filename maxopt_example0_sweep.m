%% maxopt_example0_sweep
% Sweep the parameters of a waveguide-coupled disk resonator structure.

%%% Description
% Note that this example is especially tuned for the 2D, flattened case.

function [gap, radius, P] = maxopt_example0_sweep(varargin)

        %
        % Parse inputs.
        %

    options = my_parse_options(struct(  'gap', 0.018:0.001:0.021, ... 
                                        'radius', 1.039:0.001:1.043, ...
                                        'flatten', true), ...
                                varargin, mfilename);


        % 
        % Run structures.
        %

    [gap, radius] = ndgrid(options.gap, options.radius);

    for k = 1 : numel(gap)
        fprintf('[%d/%d] ', k, numel(gap));
        [P(k), E{k}, H{k}, grid, eps{k}] = solve_ring(gap(k), radius(k), ...
                                                options.flatten);

        maxwell_view(grid, eps{k}, E{k}, 'y', [nan nan 0], 'field_phase', nan); 
        fprintf('\tgap: %1.3f, radius: %1.3f, reflected power: %1.3f\n', ...
                    gap(k), radius(k), P(k));
        drawnow
    end
    P = reshape(P, size(gap));


        %
        % Find best structure.
        %

    [P_max, ind] = max(P(:));
    gap_max = gap(ind);
    radius_max = radius(ind);
    fprintf('\nBest structure -- gap: %1.3f, radius: %1.3f, reflected power: %1.3f\n', ...
                gap_max, radius_max, P_max);
    maxwell_view(grid, eps{k}, E{k}, 'y', [nan nan 0], 'field_phase', inf); 



function [P, E, H, grid, eps] = solve_ring(gap, radius, flatten)
% Simulate a waveguide-coupled disk resonator.

        %
        % Create grid.
        %

    % Make a grid for a wavelength of 1550 nm.
    if flatten
        [grid, eps] = maxwell_grid(2*pi/1.55, -3:0.05:3, -3:0.05:3, 0); % Use this for 2D.
        m = 2;
    else
        [grid, eps] = maxwell_grid(2*pi/1.55, -3:0.05:3, -3:0.05:3, -1.5:0.05:1.5);
        m = 1;
    end


        %
        % Setup the structure.
        %

    % Structure constants.
    height = 0.2;
    si_eps = 13;
    air_eps = 1;
    wg_ypos = -2;
    wg_width = 0.4;

    % Draw waveguide
    eps = maxwell_shape(grid, eps, si_eps, ...
                        maxwell_box([0 wg_ypos 0], [inf wg_width height]));

    % Draw disk.
    eps = maxwell_shape(grid, eps, si_eps, ...
                        maxwell_cyl([0 wg_ypos+wg_width/2+radius+gap 0], radius, height));


        %
        % Solve for initial excitation.
        %

    % Excitation for the fundamental mode (of the ring's waveguide).
    J = maxwell_wgmode(grid, eps, [-1.5 wg_ypos 0], [+inf 2 2], 'mode_number', m);

    [E, H] = maxwell_solve(grid, eps, J);


        % 
        % Measure power in reflected wave.
        %

    [~, E1, H1] = maxwell_wgmode(grid, eps, [-1.7 wg_ypos 0], [-inf 2 2], 'mode_number', m);

    P = maxwell_flux(grid, [E H], [E1 H1]);
    % maxwell_view(grid, eps, E, 'y', [nan nan 0], 'field_phase', inf); % Visualize the excited waveguide.
