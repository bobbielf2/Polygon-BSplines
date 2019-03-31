% define polygon vertices
P = [ 0.1993 0.4965 0.6671 0.7085 0.6809;
    0.8377 0.8436 0.7617 0.6126 0.212];
p = 8; % polynomial degree
n = p + 1;  % some call this the order of B-spline

% uniform parameterization (put extra points on each edge)
P_temp = [P, P(:,1)];
len = sqrt(sum(diff(P_temp,1,2).^2,1)); % lengths of polygon edges
nsubdiv = round(len / min(len) * 3); % how many subdivisions on each edge
P_fine = []; % record the refined polygon
for i = 1:size(P_temp,2)-1
    P1_fine = linspace(P_temp(1,i),P_temp(1,i+1),nsubdiv(i)+1);
    P2_fine = linspace(P_temp(2,i),P_temp(2,i+1),nsubdiv(i)+1);
    P_fine = [P_fine, [P1_fine(1:end-1); P2_fine(1:end-1)]];
end
P = P_fine;

% weight of the vertices
w = zeros(1,size(P,2));
P_temp = [P(:,end), P, P(:,1)];
for i = 2:size(P_temp,2)-1
    u = P_temp(:,i) - P_temp(:,i-1);
    v = P_temp(:,i+1) - P_temp(:,i);
    
    theta = real(acos(dot(u,v)/(norm(u)*norm(v))));
    w(i-1) = ((theta/pi + 1)/2)^2;
end
% w = ones(1,size(P,2));

% periodicity (since polygon is a closed curve)
ind = mod( (1:p)-1 , size(P,2) ) + 1; % periodic index
P = [P, P(:,ind)]; % extend control points for periodicity
w = [w, w(ind)]; % extend weights accordingly

% define knot vector
m = size(P,2);
t = 1:m+n; % knot vector (query points tq satisfy n <= tq <= m+1 )
t = (t - n)/(m+1-n); % normalized knots (query points become 0 <= tq <= 1)

% B-splines evaluation

% define query points
% must satisfy t(p+1) <= tq <= t(m+1)
tq = linspace(t(n),t(end-p),1000); % just tq = linspace(0,1,1000) if normalized

% Evaluate the B-splines, require package by Hunyadi (2010)
Y = bspline_deboor(n,t,P,tq); % unweighted B-spline
X = bspline_wdeboor(n,t,P,w,tq); % weighted B-spline, aka NURBS

% b-spline derivative (unweighted)
[dt,dP] = bspline_deriv(n,t,P);
dtq = linspace(dt(n-1),dt(end-n+2),1000); % again 0 <= dtq <= 1 if normalized
dY = bspline_deboor(n-1,dt,dP,dtq);

% b-spline derivative (weighted, aka NURBS derivative)
dX = bspline_wdeboor_deriv(n,t,P,w,tq);

% plot results
figure(1); clf
subplot(131);
hold on;
plot(P(1,:), P(2,:), 'k');
plot(Y(1,:), Y(2,:), 'b');
plot(X(1,:), X(2,:), 'r');
plot(X(1,1), X(2,1), 'og');
title(['num of continuous derivative = ',num2str(p-1)])
legend({'polygon','B-spline','NURBS','starting pt'},'location','sw')

subplot(132)
hold on;
plot(tq,X(1,:))
plot(tq,dX(1,:))
plot(tq(1), X(1,1), 'og');
% plot(dX(1,:), dX(2,:),'-o');
title(['x(t) and x''(t)'])
subplot(133)
hold on;
plot(tq,X(2,:))
plot(tq,dX(2,:))
plot(tq(1), X(2,1), 'og');
title(['y(t) and y''(t)'])