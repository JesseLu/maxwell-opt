%% example2_adjoint
% Derivative-based optimization of an H1 photonic crystal resonator.


function [Emag, params] = example2_adjoint(varargin)

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

    E = solve_resonator([5, 5], zeros(25, 2), options.flatten);
end

function [fval, df_dp, E, H, grid, eps] = solve_resonator(pc_size, shifts, flatten)
% Simulate a photonic crystal resonator.

        %
        % Build the structure.
        %

    [grid, eps, J] = make_resonator_structure(pc_size, shifts, flatten);


        %
        % Use a central point excitation.
        %

    [x, y, z] = maxwell_pos2ind(grid, 'Ey', [0 0 0]); % Get central Ey component.
    x = x-1; % Slight adjustment.
    y = y-1;
    J{2}(x+[0 1], y, z) = 1;

        
        %
        % Solve.
        %

    if ~flatten; figure(2); end
    [E, H] = maxwell_solve(grid, eps, J);

    figure(1);
    maxwell_view(grid, eps, E, 'y', [nan nan 0]); % Visualize.


        % 
        % Measure power reflected back to the center (figure of merit).
        %

    E_meas = [E{2}(x, y, z); E{2}(x+1, y, z)];
    fval = norm(E_meas).^2; % This is the figure of merit.


        % 
        % Calculate gradient.
        %

    Egrad = my_default_field(grid.shape, 0); 
    Egrad{2}(x, y, z) = E{2}(x, y, z);
    Egrad{2}(x+1, y, z) = E{2}(x+1, y, z);

    function [eps] = make_eps(params)
        [~, eps] = make_resonator_structure(pc_size, params, flatten);
    end

    df_dp = maxopt_structgrad(grid, eps, E, Egrad, @make_eps, shifts);
end

function [grid, eps, J] = make_resonator_structure(pc_size, shifts, flatten)

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
    pos = {};
    for i = 1 : pc_size(1)
        for j = 1 : pc_size(2)
            p{1} = a * [(i-0.5) (j-0.5) 0] + [shifts(1) shifts(2)*(i~=1) 0];
            p{2} = p{1} .* [-1 1 1];
            p{3} = p{1} .* [1 -1 1];
            p{4} = p{1} .* [-1 -1 1];
            pos = [pos, p];
        end
    end

    for k = 1 : length(pos)
        eps = maxwell_shape(grid, eps, air_eps, ...
                            maxwell_cyl_smooth(pos{k}, radius, 2*height, ...
                                                'smooth_dist', 0.025));
    end

%     if flatten
%         subplot 211; maxwell_view(grid, eps, [], 'y', [nan 0 nan]);
%     else
%         subplot 111; maxwell_view(grid, eps, [], 'y', [nan nan height/4]);
%         % subplot 222; maxwell_view(grid, eps, [], 'y', [nan 0 nan]);
%     end

end
