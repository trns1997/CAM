function [d]=d3Dsl(p,a,x0)

% function [d]=d3Dsl(p,a,x0)
% this function calculates the distance of a set of n points from a 3D
% straight line
% input: p, n x 3 matrix, each row containing the [x y z] coordinates of a
%        point
%        a, vector [a b c] defining the direction of the line
%        x0, vector [x0 y0 z0] of a point beliìonging to the line
% output; d, n x 1 vector of the distances

% determine the number of points

n=size(p,1);

% standardize the a vector

a=a/norm(a);

% calculate cross products

crosses=cross(repmat(a(:)',n,1),p-repmat(x0(:)',n,1));

% calculate the norms

d=zeros(n,1);
for h=1:n
    d(h)=norm(crosses(h,:));
end