function [a, b] = my_split(a_b, shape, varnames, filename)
    if length(a_b) == 6
        a = a_b(1:3);
        b = a_b(4:6);   
    else
        a = a_b;
        b = [];
    end

    my_validate_field(a, shape, varnames{1}, filename)
    if ~isempty(b)
        my_validate_field(a, shape, varnames{2}, filename)
    end

