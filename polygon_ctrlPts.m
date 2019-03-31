function [t, P, w] = polygon_ctrlPts(P, n)
% Input:
%   P   : polygon vertices, 2-m matrix
%   n   : order of the B-splines
% Output:
%   t   : knot vector
%   P   : control points
%   w   : weights

% Bowei Wu, 2019/3/31

% polynomial degree
p = n - 1; 

% uniform parameterization (put extra points on each edge)
P = polygon_refine(P);

% weight of the vertices
w = polygon_weight(P);

% Periodic extension of control points and weights 
%   according to polynomial degree (since polygon is a closed curve)

if nargin > 2
    assert(size(P,2) == size(w,2), 'polygon:ctrlPts:DimensionMismatch', ...
        'Number of weights must equal to the number of vertices.');
end
ind = mod( (1:p)-1 , size(P,2) ) + 1; % periodic index
P = [P, P(:,ind)]; % extend control points for periodicity
if nargout > 2
    w = [w, w(ind)]; % extend weights accordingly
end

% define knot vector
m = size(P,2);
t = 1:m+n; % knot vector (query points tq satisfy n <= tq <= m+1 )
t = (t - n)/(m+1-n); % normalized knots (query points become 0 <= tq <= 1)


function P_fine = polygon_refine(P)
% uniform parameterization (put extra points on each edge)

P = [P, P(:,1)];
len = sqrt(sum(diff(P,1,2).^2,1)); % lengths of polygon edges
nsubdiv = round(len / min(len) * 3); % how many subdivisions on each edge
P_fine = []; % record the refined polygon
for i = 1:size(P,2)-1
    P1_fine = linspace(P(1,i),P(1,i+1),nsubdiv(i)+1);
    P2_fine = linspace(P(2,i),P(2,i+1),nsubdiv(i)+1);
    P_fine = [P_fine, [P1_fine(1:end-1); P2_fine(1:end-1)]];
end


function w = polygon_weight(P,k)
% assign weights to polygon vertices based on exterior angles
if nargin < 2, k = 2; end
w = zeros(1,size(P,2));
P = [P(:,end), P, P(:,1)];
for i = 2:size(P,2)-1
    u = P(:,i) - P(:,i-1);
    v = P(:,i+1) - P(:,i);
    
    theta = real(acos(dot(u,v)/(norm(u)*norm(v))));
    w(i-1) = ((theta/pi + 1)/2)^k;
end
% w = ones(1,size(P,2));