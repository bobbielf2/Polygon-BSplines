% define polygon vertices
P = [ 0.1993 0.4965 0.6671 0.7085 0.6809;
    0.8377 0.8436 0.7617 0.6126 0.212];
n = 9;  % order of B-spline (note: polynomial degree p = n - 1)

% polygon control points, knot vector, and weights
[t, P, w] = polygon_ctrlPts(P, n);

% B-splines evaluation

% define query points { must satisfy t(n) <= tq <= t(m+1) }
tq = linspace(t(n),t(end-n+1),1000); % just tq = linspace(0,1,1000) if normalized

% eval B-splines, require package by Hunyadi (2010)
Y = bspline_deboor(n,t,P,tq); % unweighted
X = bspline_wdeboor(n,t,P,w,tq); % weighted, aka NURBS

% eval B-splines derivative 
dY = bspline_deboor_deriv(n,t,P,tq); % unweighted
dX = bspline_wdeboor_deriv(n,t,P,w,tq); % weighted, aka NURBS derivative

% plot results
figure(1); clf
subplot(231);
hold on;
plot(P(1,:), P(2,:), 'k');
plot(Y(1,:), Y(2,:), 'b');
plot(X(1,:), X(2,:), 'r');
plot(X(1,1), X(2,1), 'og');
title(['B-splines order = ',num2str(n)])
legend({'polygon','B-spline','NURBS','starting pt'},'location','sw


subplot(232)
hold on;
plot(tq,X(1,:))
plot(tq,dX(1,:))
plot(tq(1), X(1,1), 'og');
title('x(t) and x''(t), NURBS')
subplot(233)
hold on;
plot(tq,X(2,:))
plot(tq,dX(2,:))
plot(tq(1), X(2,1), 'og');
title('y(t) and y''(t), NURBS')

subplot(235)
hold on;
plot(tq,Y(1,:))
plot(tq,dY(1,:))
plot(tq(1), Y(1,1), 'og');
title('x(t) and x''(t), B-splines')
subplot(236)
hold on;
plot(tq,Y(2,:))
plot(tq,dY(2,:))
plot(tq(1), Y(2,1), 'og');
title('y(t) and y''(t), B-splines')