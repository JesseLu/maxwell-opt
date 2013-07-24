function [optim_step] = my_linesearch(f, f0, x0, dx, varargin) 
% Backtracking line search.

    options = my_parse_options(struct(  'alpha', 0.9, ...

                                        'solver_args', {{}}, ...
                                        'fitness', [], ...
                                        'check_gradients', false), ...
                                varargin, mfilename);