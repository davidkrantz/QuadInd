function [x, w] = GaussLaguerre8()
    % GAUSSLAGUERRE8 8-point Gauss-Laguerre quadrature nodes and weights
    %
    %   [x, w] = GaussLaguerre8() returns the 8-point Gauss-Laguerre
    %   quadrature nodes and weights for integrals of the form:
    %
    %       integral_0^inf f(x) * exp(-x) dx ≈ sum(f(x) .* w)
    %
    % Outputs:
    %   x - Row vector of 8 quadrature nodes
    %   w - Row vector of 8 quadrature weights
    %
    % Note:
    %   These are hardcoded high-precision values. The 8-point rule
    %   is sufficient for the error-indicator evaluation integrals in this package.
    %
    % See also: quadind.indicator.UniformIndicatorBuilder
    
    x = [
        0.170279632305101
        0.903701776799380
        2.251086629866131
        4.266700170287659
        7.045905402393466
        10.758516010180996
        15.740678641278004
        22.863131736889265
    ]';
    
    w = [
        0.369188589341639
        0.418786780814343
        0.175794986637171
        0.033343492261216
        0.002794536235226
        0.000090765087734
        0.000000848574672
        0.000000001048001
    ]';
end
