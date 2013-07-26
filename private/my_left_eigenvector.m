function [y] = my_left_eigenvector(grid, x)

    % Form elements of diagonal symmetrization matrix S.
    [spx, spy, spz] = ndgrid(   s2D(grid.s_prim{1}), ...
                                s2D(grid.s_prim{2}), ...
                                s2D(grid.s_prim{3}));

    [sdx, sdy, sdz] = ndgrid(   s2D(grid.s_dual{1}), ...
                                s2D(grid.s_dual{2}), ...
                                s2D(grid.s_dual{3}));

    my_diag = @(z) spdiags(z(:), 0, numel(z), numel(z));
    S = my_diag([   sdx(:).*spy(:).*spz(:); ...
                    spx(:).*sdy(:).*spz(:); ...
                    spx(:).*spy(:).*sdz(:)]);

    % Obtain right eigenvector.
    y = conj(S * x);
    y = y / norm(y);

function [s] = s2D(s)
    if any(isinf(s))
        s = 1;
    end
