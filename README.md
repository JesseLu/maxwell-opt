# maxwell-opt

Optimization using maxwell (maxwell-client)

## Example usage

    % User-supplied function which defines struxture
    function [grid, eps] = my_structure(params{:})
        ...
    end
    
    % User-supplied function which solves the simulation
    function [E, H] = my_sim(grid, eps)
        ...
    end
    
    % Build the figure of merit
    % Need functionality here to chec the field of merit.
    [f, dfdx] = mopt_fom();
    
    % Obtain the function for computing the structural derivative (optional).
    dfdz = mopt_zgrad(f, dfdx, @my_sim, @my_structure);
    
    % Feed into Matlab generic optimization functiona.
    
