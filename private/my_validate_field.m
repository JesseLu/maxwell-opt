function my_validate_field(F, shape, varname, filename)

    my_validate = @(var, type, attr, var_name) ...
        validateattributes(var, type, attr, filename, var_name);

    my_validate(F, {'cell'}, {'numel', 3}, varname);
    for k = 1 : 3
        my_validate(F{k}, {'double'}, {'size', shape}, sprintf('%s{%d}', varname, k));
    end
