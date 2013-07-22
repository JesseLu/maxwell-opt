function [F] = my_default_field(shape, val)
    for k = 1 : 3
        F{k} = val * ones(shape);
    end
