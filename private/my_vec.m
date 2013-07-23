
function [vec, unvec] = my_vec(shape)
% Here, A1 referes to curl on H-fields, and A2 refers to curl on E-fields.

    n = prod(shape);

    % Inline functions to vectorize and unvectorize.
    function [z] = vec_fun(z)
        z = [z{1}(:); z{2}(:); z{3}(:)];
    end

    function [z] = unvec_fun(z0)
        [z{1}, z{2}, z{3}] = deal(  reshape(z0(1:n), shape), ...
                                    reshape(z0(n+1:2*n), shape), ...
                                    reshape(z0(2*n+1:3*n), shape));
    end

    [vec, unvec] = deal(@vec_fun, @unvec_fun);
end
