%% maxopt_structgrad
% Structural gradient.

%%% Syntax

%%% Description
%
%
% Note that we assume mu = 1 here.
%

function [dp] = maxopt_structgrad(grid, eps, E, Egrad, make_eps)

        %
        % Check and validate inputs.
        %

    my_validate_grid(grid, mfilename);

    my_validate_field(eps, grid.shape, 'eps', mfilename);

    my_validate_field(E, grid.shape, 'E', mfilename);
    
    my_validate_field(Egrad, grid.shape, 'Egrad', mfilename);

    validateattributes(make_eps, {'function_handle'}, {}, 'make_eps', mfilename);


        %
        %
