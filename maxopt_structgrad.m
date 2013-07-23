%% maxopt_structgrad
% Structural gradient.

%%% Syntax

%%% Description
%
%
% Note that we assume mu = 1 here.
%

function [df_dp] = maxopt_structgrad(grid, eps, E, Egrad, make_eps, p0, varargin)

        %
        % Check and validate inputs.
        %

    my_validate_grid(grid, mfilename);

    my_validate_field(eps, grid.shape, 'eps', mfilename);

    my_validate_field(E, grid.shape, 'E', mfilename);
    
    my_validate_field(Egrad, grid.shape, 'Egrad', mfilename);

    validateattributes(make_eps, {'function_handle'}, {}, 'make_eps', mfilename);

    % Optional arguments.
    options = my_parse_options(struct(  'delta_p', 1e-6), ...
                                varargin, mfilename);

    validateattributes(options.delta_p, {'numeric'}, ...
                        {'nonnan', 'finite', 'real'}, ...
                        'delta_p', mfilename);

        %
        % Initiate adjoint solve.
        %

    fprintf('Initiating adjoint field solve...\n');
    for k = 1 : 3
        Egrad{k} = -conj(Egrad{k}); % Set up to calculate -(A^-T Egrad^*).
    end
    cb = maxopt_solveAT(grid, eps, Egrad);


        %
        % Calculate dz/dp.
        %
    
    % Print progress.
    fprintf('Calculating structural gradient... ');
    progress_text = sprintf('(%d/%d)', 0, numel(p0));
    fprintf(progress_text);

    eps0 = make_eps(p0);
    for k = 1 : numel(p0)
        p = p0;
        p(k) = p(k) + options.delta_p;
        eps1 = make_eps(p);
        dz_dp(:,k) = sparse([   eps1{1}(:)-eps0{1}(:); ...
                                eps1{2}(:)-eps0{2}(:); ...
                                eps1{3}(:)-eps0{3}(:)]);

        % Print progress.
        fprintf(repmat('\b', 1, length(progress_text)));
        progress_text = sprintf('(%d/%d)', k, numel(p0));
        fprintf(progress_text); 
    end
    fprintf('\n');
    

        %
        % Put it all together.
        %

    % Wait for the solve field.
    fprintf('Waiting for adjoint solve to complete ...');
    while ~cb(); end
    [~, Ea] = cb(); % "Adjoint field".

    % Form the material derivative.
    df_dz = (conj([Ea{1}(:); Ea{2}(:); Ea{3}(:)]) .* ...
            (-grid.omega^2 * [E{1}(:); E{2}(:); E{3}(:)])).';

    % Form the structural (parameter) derivative.
    df_dp = real([df_dz * dz_dp].');
    

