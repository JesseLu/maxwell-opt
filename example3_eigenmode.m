%% example3_eigenmode
% Derivative-based optimization of an L3 cavity mode.

function [fval, x, f_vis] = example3_eigenmode(case_name, varargin)

        %
        % Parse inputs.
        %

    options = my_parse_options(struct(  'delta', 0.5 * [0 ones(1, 5)], ... 
                                        'width', 0.1 * ones(1, 6), ...
                                        'sim_only', false, ...
                                        'flatten', false), ...
                                varargin, mfilename);


        %
        % Set up the optimization problem.
        %


    switch case_name
        case 'L3'
            [fun, x] = case3_L3('grad_f', 'flatten', options.flatten);
            [f_vis] = case3_L3('get_fields', 'flatten', options.flatten);
        otherwise
            error('Invalid case_name.');
    end


    function vis_progress(hist)
        fprintf('fval: %e\n', hist(end));
        figure(2);
        plot(hist, '.-');
        xlabel('optimization iterations');
        ylabel('fval');
        title('structure optimization progress');
    end

    if ~options.sim_only
        [x, fval, hist] = my_grad_descent(fun, x,   'init_step', 0.1, ...
                                                    'max_delta', 0.1, ...
                                                    'vis_progress', @vis_progress);
    end

end
