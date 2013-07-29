function test_eiggrad(n)

    % Create a problem
    A0 = randn(n) + 1i * randn(n);
    A0 = A0' * A0;
    A = @(z) A0 + diag(z(:));
    B = @(x) diag(x(:));
    b = randn(n, 1) + 1i * randn(n, 1);
    z = randn(n, 1);

    % Choose a eigenvector.
    [V, D] = eig(A(z));
    ind = randi(n);
    v = V(:,ind);
    w = v;
    lambda = D(ind, ind);

    % Function to find nearest eigenvector.
    function [lambda, v] = my_eig(lambda, z)
        [v, dlambda] = eigs(A(z) - lambda * eye(n), 1, 'SM');
        lambda = lambda + dlambda;
    end

    function [v] = my_eig2(lambda, z)
        [~, v] = my_eig(lambda, z);
        v = v ./ v(1)
    end

    [lambda, v] = my_eig(lambda, z);

    % Test the eigenvector.
    norm((A(z) - lambda*eye(n)) * v) / norm(v)
    norm(w' * (A(z) - lambda*eye(n))) / norm(w)

    % Calculate dlambda/dz.
    dl_dz = w' * B(v) / (w' * v);
    my_gradient_test(@(z) my_eig(lambda, z), dl_dz, z, 'real', 'dlambda/dz');

%     % Calculate dv/dz (this doesn't scale in real cases).
%     dv_dz = (A(z) - lambda*eye(n)) \ ((eye(n) - (1/(w' * v)) * v * w') * B(v));
%     my_gradient_test(@(z) my_eig2(lambda, z), dv_dz, z, 'eig', 'dv/dz');
end 
