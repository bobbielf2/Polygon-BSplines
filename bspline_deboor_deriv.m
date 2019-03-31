function [C, u] = bspline_deboor_deriv(n,t,P,u)
% Evaluate derivative of explicit B-spline at specified locations.
%
% Input arguments:
% n:
%    B-spline order (2 for linear, 3 for quadratic, etc.)
% t:
%    knot vector
% P:
%    control points, typically 2-by-m, 3-by-m, or 4-by-m (for weights)
% u (optional):
%    values where the B-spline is to be evaluated, or a positive
%    integer to set the number of points to automatically allocate
% Output arguments:
% C:
%    points of the B-spline curve

% 2019 Bowei Wu, based on Levente Hunyadi 2010 B-splines package

[dt,dP] = bspline_deriv(n,t,P);
if nargin < 4
    u = linspace(dt(n-1),dt(end-n+2),1000); % again 0 <= u <= 1 if normalized
end
C = bspline_deboor(n-1,dt,dP,u);
