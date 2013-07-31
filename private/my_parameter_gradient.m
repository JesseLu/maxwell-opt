function [grad_p] = my_parameter_gradient(p2z, params0, delta_p)

    % Find the dz/dp derivative.
    z0 = p2z(params0);
    progress_text = '';
    for k = 1 : numel(params0) 
        p = params0;
        p(k) = p(k) + delta_p;
        z = p2z(p);
        grad_p(:, k) = sparse((z - z0) ./ delta_p);

        fprintf(repmat('\b', 1, length(progress_text)));
        progress_text = sprintf('[%d/%d gradients computed] ', k, numel(params0));
        fprintf(progress_text);

%         % Debug.
%         maxwell_view(grid, unvec((z - z0) ./ options.delta_p), [], 'y', [nan nan 0], 'clims', [-1 1]);
%         pause
    end


