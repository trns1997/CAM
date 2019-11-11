function [d]=d3Dp(p,a)

% function [d]=d3Dp(p,a)
% this function calculates the distance of a set of n points from a 3D
% plane
% input: p, n x 3 matrix, each row containing the [x y z] coordinates of a
%        point
%        a, vector [a b c d] defining the plane parameters
% output; d, n x 1 vector of the distances

% standardize the a vector

a=a(:)';
a=a/norm(a(1:3));

% calculate distances

d=abs(p*a(1:3)'+a(4));

