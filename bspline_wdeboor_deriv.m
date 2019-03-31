function C = bspline_wdeboor_deriv(n,t,P,w,u)
% Evaluate derivative of explicit weighed B-spline at specified locations.
% (a.k.a., the derivative of NURBS)
%
% Input arguments:
% n:
%    B-spline order (2 for linear, 3 for quadratic, etc.)
% t:
%    knot vector
% P:
%    control points, typically 2-by-m, 3-by-m, or 4-by-m (for weights)
% w:
%    weight vector
% u:
%    values where the B-spline is to be evaluated, or a positive
%    integer to set the number of points to automatically allocate
% Output arguments:
% C:
%    points of the B-spline curve

% 2019 Bowei Wu, based on Levente Hunyadi 2010 B-splines package

% first term = { \sum_i w_i N'_{i,p}(xi) P_i } / W(xi) = = { \sum_i N_{i,p-1}(xi) (P .* w)'_i } / W(xi)
w = transpose(w(:));
Pw = bsxfun(@times, P, w);
[dt,dPw] = bspline_deriv(n,t,Pw);
C1 = bspline_deboor(n-1,dt,dPw,u); % not normalize yet!

% second term = { \sum_i w_i N_{i,p}(xi) P_i } * W'(xi) / W(xi)^2
[dt,dw] = bspline_deriv(n,t,w);
dW = bspline_deboor(n-1,dt,dw,u);
C2W = bspline_deboor(n,t,[Pw; w],u);
C2 = C2W(1:end-1,:);
W = C2W(end,:);

% apply correct weights
C2 = bsxfun(@times, C2, dW ./ W.^2);
C1 = bsxfun(@rdivide, C1, W);

% final result
C = C1 - C2;