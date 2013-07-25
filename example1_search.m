%% example1_search
% Derivative-free optimization of a nanophotonic grating coupler.


function [fval, x0, f_vis] = example1_search(case_name, varargin)

        %
        % Parse inputs.
        %

    options = my_parse_options(struct(  'iters', 200, ...
                                        'flatten', false), ...
                                varargin, mfilename);

        
        %
        % Set up the optimization problem.
        %

    switch case_name
        case 'grating'
            [f, x0] = case1_grating('fval', 'flatten', options.flatten);
            [f_vis] = case1_grating('get_fields', 'flatten', options.flatten);
        otherwise
            error('Invalid case_name.');
    end

    search_options = optimset(  'Display', 'iter', ...
                                'TolX', 0.05, ...
                                'Tolfun', 1e-16, ...
                                'MaxFunEvals', 10 * options.iters, ...
                                'MaxIter', options.iters, ...
                                'FunValCheck', 'on', ...
                                'PlotFcns', {[]}, ...
                                'OutputFcn', {[]});


        %
        % Perform the optimization.
        %

    [x, fval] = fminsearch(f, x0, search_options);
end


