%% Get Best Fit Plane
clear all, clc, close all;
load('data_1.mat');
points_pl = planes{1};
fun=@(params)d3Dp(points_pl,params);
a=lsqnonlin(fun,[0 0 1 0]);

z = @(x,y)((-a(1)*x - a(2)*y - a(4))/a(3));

plot3(points_pl(:, 1), points_pl(:, 2), points_pl(:, 3), 'o')
hold on
fsurf(z)
axis square

%% Transformation Matrix
num_cyl = 10;
points_cyl = cylinders{1};
plot3(points_cyl(:, 1), points_cyl(:, 2), points_cyl(:, 3), 'o')

q = [0, 0, (-a(1)*0 - a(2)*0 - a(4))/a(3)];
r = [0, 100, (-a(1)*0 - a(2)*100 - a(4))/a(3)];
s = [100, 0, (-a(1)*100 - a(2)*0 - a(4))/a(3)];

N = cross((r-q), (s-q));
N = N/norm(N);

localz = N;

localx = r-q;   
unitx = localx/norm(localx);

localy = cross(localz, localx);  
unity = localy/norm(localy); 

T = [localx(:), localy(:), localz(:), q(:); 0 0 0 1];

%% Project and Determine Centroids
for j = 0:num_cyl-1
    for i = j*(length(points_cyl)/num_cyl)+1:(j+1)*(length(points_cyl)/num_cyl)
        v = points_cyl(i,:) - q;
        dist = dot(v, N);
        k = i - j*(length(points_cyl)/num_cyl);
        proj_p(k,:) = points_cyl(i,:) - dist*N;
    end
%     plot3(proj_p(:,1),proj_p(:,2),proj_p(:,3), 'o')

    C = [proj_p, ones(158,1)];
    Coor_2D = T \ C';
    Coor_2D = Coor_2D(1:2,:)';

    fun=@(params)(d2Dpp(Coor_2D,params(1:2))-params(3));
    opt_par=lsqnonlin(fun,[0 0 10]);
    % figure; axis equal;hold on
    % plot(Coor_2D(:,1),Coor_2D(:,2),'.');
    pc=[linspace(0,2*pi,100)' repmat(opt_par(3),100,1)];
    [pc(:,1),pc(:,2)]=pol2cart(pc(:,1),pc(:,2));
    pc=[pc(:,1)+opt_par(1),pc(:,2)+opt_par(2)];

    polyin = polyshape({pc(:,1)}, {pc(:,2)});
    [x,y] = centroid(polyin);
    % plot(polyin);
    % plot(x,y,'r*')

    Coor_3D = T*[x,y,0,1]';
    p(j+1,:) = Coor_3D(1:3)' + dist*N;
end

plot3(p(:,1),p(:,2), p(:,3), 'o')

%% Determine Deviation

fun=@(params)(d3Dsl(p,N,params(1:3)));
x0=fminimax(fun,[0 0 0]);
dev=2*max(fun(x0));
p0 = x0 + N*80;
p1 = [x0; p0];
plot3(p1(:,1),p1(:,2),p1(:,3), 'LineWidth',2)

caption = sprintf('Deviation = %d', dev);
title(caption, 'FontSize', 20);

