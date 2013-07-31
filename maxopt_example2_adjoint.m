%% maxopt_example2_adjoint
% Form a cavity out of a square lattice using gradient optimization.


function [fval, x, f_vis] = maxopt_example2_adjoint(varargin)

    case_name = 'squarepc';

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
        case 'squarepc'
            [f, x0] = maxopt_case_squarepc('grad_f', 'flatten', options.flatten);
            [f_vis] = maxopt_case_squarepc('get_fields', 'flatten', options.flatten);
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
