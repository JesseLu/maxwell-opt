%% maxopt_solve_adjoint
% Asynchronously solve for the adjoint of the operator.

%%% Description
% Really just a simple wrapper around maxwell_solve_async.
%
% This works via a symmetrization matrix S.
%


function [cb] = maxopt_solve_adjoint(grid, eps_mu, J, varargin)

        %
        % Validate and parse inputs.
        %

    my_validate_grid(grid, mfilename);

    [eps, mu] = my_split(eps_mu, grid.shape, {'eps', 'mu'}, mfilename);
    if isempty(mu)
        mu = my_default_field(grid.shape, 1);
    end
    my_validate_field(eps, grid.shape, 'eps', mfilename);
    my_validate_field(mu, grid.shape, 'mu', mfilename);

    my_validate_field(J, grid.shape, 'J', mfilename);


        %
        % Prepare elements for adjoint solve.
        %

    % Form elements of diagonal symmetrization matrix S.
    [spx, spy, spz] = ndgrid(   s2D(grid.s_prim{1}), ...
                                s2D(grid.s_prim{2}), ...
                                s2D(grid.s_prim{3}));

    [sdx, sdy, sdz] = ndgrid(   s2D(grid.s_dual{1}), ...
                                s2D(grid.s_dual{2}), ...
                                s2D(grid.s_dual{3}));

    s = {   sdx.*spy.*spz, ...
            spx.*sdy.*spz, ...
            spx.*spy.*sdz};

    t = {   spx.*sdy.*sdz, ...
            sdx.*spy.*sdz, ...
            sdx.*sdy.*spz};
     
    % Conjugate all the elements which compose the A matrix.
    for k = 1 : 3
        grid.s_prim{k} = conj(grid.s_prim{k});
        grid.s_dual{k} = conj(grid.s_dual{k});
        eps{k} = conj(eps{k});
        mu{k} = conj(mu{k});

        % Invert by conjugate of s.
        % Also, take care of the problem that omega needs to be conjugated
        % for matrix A, but not for b.
        J{k} = (grid.omega) / conj(grid.omega) * (J{k}) ./ conj(s{k}); 
    end
    grid.omega = conj(grid.omega);


        % 
        % Initiate asynchronous solve.
        %

    % Because of the conjugated terms, this actually solves A conjugate.
    cb_orig = maxwell_solve_async(grid, [eps mu], J);

    % Special callback...
    cb = @() my_adjoint_callback(cb_orig, s, t);

end 

function [is_done, E, err] = my_adjoint_callback(cb_orig, s, t)
    [is_done, E, ~, err] = cb_orig();
    if is_done
        for k = 1 : 3 % Re-transform by s.
            E{k} = conj(s{k}) .* E{k};
        end
    end
end

function [s] = s2D(s)
    if any(isinf(s))
        s = 1;
    end
end
