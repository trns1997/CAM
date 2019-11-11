%% Es 1: Least squares roundness

clear;close all;clc

load roundness
fun=@(params)(d2Dpp(points,params(1:2))-params(3));
opt_par=lsqnonlin(fun,[0 0 10]);
roundness=range(d2Dpp(points,opt_par(1:2)));
disp(['Roundness = ' num2str(ref_roundness)]);
figure; axis equal;hold on
plot(points(:,1),points(:,2),'.');
pc=[linspace(0,2*pi,100)' repmat(opt_par(3),100,1)];
[pc(:,1),pc(:,2)]=pol2cart(pc(:,1),pc(:,2));
pc=[pc(:,1)+opt_par(1),pc(:,2)+opt_par(2)];
plot(pc(:,1),pc(:,2),'g');

%% Es 2: MZ straightness

clear;close all;clc

load straightness
fun=@(params)d2Dsl(points,params);
opt_par=fminimax(fun,[0 1 0]);
straightness=2*max(fun(opt_par));
figure; axis equal;hold on
plot(points(:,1),points(:,2),'.');
pc=linspace(-5,5,100)';
ordinate=-opt_par(1)/opt_par(2)*pc-opt_par(3);
plot(pc,ordinate,'r');
[~,maximum]=max(fun(opt_par));
d=-(opt_par(1)*points(maximum,1)+opt_par(2)*points(maximum,2));
ordinate=-opt_par(1)/opt_par(2)*pc-d;
plot(pc,ordinate,'g');
ordinate=-opt_par(1)/opt_par(2)*pc-(2*opt_par(3)-d);
plot(pc,ordinate,'g');
disp(['Straightness = ' num2str(ref_straightness)]);

%% Es 3: MC cylinder

clear;close all;clc

load cylindricity
fun=@(params)(d3Dsl(points,params(1:3),params(4:6)));
cylinder=fminimax(fun,[0 0 1 0 0 0]);
radius=max(fun(cylinder));
disp(['Reference cylinder parameters = ' num2str(ref_cylinder)]);


%% Es 4: planes perpendicularity

clear;close all;clc

load perpendicularity
fun=@(params)d3Dp(datum,params);
a=lsqnonlin(fun,[0 0 1 0]);
fun=@(params)d3Dp(feature,params);
b=fminimax(fun,[0 1 0 0],[],[],[a(1:3) 0],0);
perpendicularity=2*max(fun(b));
figure; hold on; axis equal
a=a/norm(a(1:3));
b=b/norm(b(1:3));
plot3(datum(:,1),datum(:,2),datum(:,3),'.r');
plan=[-6 -6
    -6 6
    6 6
    6 -6];
plan(:,3)=(-plan(:,1)*a(1)-plan(:,2)*a(2)-a(4))/a(3);
patch('Faces',[1 2 3 4],'Vertices',plan,'FaceColor','r','facealpha',0.4);
plot3(feature(:,1),feature(:,2),feature(:,3),'.b');
plan=[-6 0 -6
    -6 0 6
    6 0 6
    6 0 -6];
% plan(:,2)=(-plan(:,1)*b(1)-plan(:,3)*b(3)-b(4))/b(2);
% patch('Faces',[1 2 3 4],'Vertices',plan,'FaceColor','b','facealpha',0.4);
delta=perpendicularity/2;
plan(:,2)=(-plan(:,1)*b(1)-plan(:,3)*b(3)-b(4)+delta)/b(2);
patch('Faces',[1 2 3 4],'Vertices',plan,'FaceColor','b','facealpha',0.2);
plan(:,2)=(-plan(:,1)*b(1)-plan(:,3)*b(3)-b(4)-delta)/b(2);
patch('Faces',[1 2 3 4],'Vertices',plan,'FaceColor','b','facealpha',0.2);

%% Es 5: location tolerance

clear;close all;clc

load location
fun=@(params)d3Dp(datum1,params);
a=lsqnonlin(fun,[0 0 1 0])';a=a/norm(a(1:3));
fun=@(params)d3Dp(datum2,params);
b=lsqnonlin(fun,[0 1 0 0])';b=b/norm(b(1:3));
fun=@(params)d3Dp(datum3,params);
c=lsqnonlin(fun,[1 0 0 0])';c=c/norm(c(1:3));
if a(3)<0
    a=-a;
end
if b(2)<0
    b=-b;
end
if c(1)<0
    c=-c;
end
b1=b(1:3)-a(1:3)*(a(1:3)'*b(1:3));b1=b1/norm(b1);
c1=cross(b1,a(1:3));
rotmat=[c1,b1,a(1:3)];
x0=[a(1:3)';b(1:3)';c(1:3)']\[-a(4);-b(4);-c(4)];
frt=(feature-repmat(x0',size(feature)./size(x0')))*rotmat;
location=2*max(d2Dpp(frt(:,1:2),[100 68]));






