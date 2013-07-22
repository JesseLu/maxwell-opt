%% test_solveAT
% Test maxopt_solveAT.

function [] = test_solveAT(varargin)

    options = my_parse_options(struct(  'delta', 0.5 * [0 ones(1, 5)], ... 
                                        'width', 0.1 * ones(1, 6), ...
                                        'flatten', false), ...
                                varargin, mfilename);
    
    x = -2:0.05:2;
    y = -2:0.05:2;
    z = -2:0.05:2;
    m = 1;
    if options.flatten
        z = 0;
        m = 2;
    end

    [grid, eps, ~, J] = maxwell_grid(2*pi/1.55, x, y, z);

    eps = maxwell_shape(grid, eps, 13, maxwell_box([0 0 0], [inf .4 .2]));

    J = maxwell_wgmode(grid, eps, [0 0 0], [+inf 2 2], 'mode_number', m);

    % cb = maxwell_solve_async(grid, eps, J);
    cb = maxopt_solveAT(grid, eps, J);

    while ~cb(); end;

    [~, E] = cb();

    [A, x, b] = maxwell_axb(grid, eps, E, J);
    norm(A'*x-b) / norm(b)

    maxwell_view(grid, eps, E, 'y', [nan nan 0], 'field_phase', inf);    
