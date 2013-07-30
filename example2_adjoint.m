%% example2_adjoint
% Form a cavity out of a square lattice using gradient optimization.


function [fval, x, f_vis] = example2_adjoint(case_name, varargin)

        %
        % Parse inputs.
        %

    options = my_parse_options(struct(  'iters', 100, ...
                                        'flatten', false), ...
                                varargin, mfilename);


        %
        % Set up the optimization problem.
        %

    switch case_name
        case 'square_pc'
            [f, x0] = case2_square('grad_f', 'flatten', options.flatten);
            [f_vis] = case2_square('get_fields', 'flatten', options.flatten);
        otherwise
            error('Invalid case_name.');
    end


    % Visualization function for optimization progress.
    function vis_progress(hist)
        figure(2);
        plot(hist, '.-');
        xlabel('optimization iterations');
        ylabel('fval');
        title('structure optimization progress');
    end

        
        %
        % Perform the optimization.
        %
        
    [x, fval, hist] = maxopt_gradient_descent(f, x0, ...
                                                'init_step', 0.1, ...
                                                'max_delta', 0.1, ...
                                                'max_iters', options.iters, ...
                                                'vis_progress', @vis_progress);

end
