function f = my_validate_grid(grid, name)

    my_validate = @(var, type, attr, var_name) ...
        validateattributes(var, type, attr, name, var_name);

    my_validate(grid, {'struct'}, {'nonempty'}, 'grid')

    if ~isfield(grid, 'omega') || ~isfield(grid, 'shape') || ...
        ~isfield(grid, 's_prim') || ~isfield(grid, 's_dual') || ...
        ~isfield(grid, 'origin') || length(fieldnames(grid)) ~= 5
        error(['grid variable must be a structure with these 5 fields: ', ...
                'omega, shape, origin, s_prim, s_dual.']);
    end

    my_validate(grid.omega, {'double'}, {'nonnan', 'finite', 'scalar'}, 'grid.omega');
    my_validate(grid.shape, {'numeric'}, {'real', 'integer', 'numel', 3}, 'grid.shape');
    my_validate(grid.origin, {'numeric'}, {'real', 'numel', 3}, 'grid.origin');
    my_validate(grid.s_prim, {'cell'}, {'numel', 3}, 'grid.s_prim');
    my_validate(grid.s_dual, {'cell'}, {'numel', 3}, 'grid.s_dual');
    for k = 1 : 3
        my_validate(grid.s_prim{k}, {'double'}, {'numel', grid.shape(k)}, sprintf('grid.s_prim{%d}', k));
        my_validate(grid.s_dual{k}, {'double'}, {'numel', grid.shape(k)}, sprintf('grid.s_dual{%d}', k));
    end
end
