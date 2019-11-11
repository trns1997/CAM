function [d]=d2Dpp(p,x0)

% function [d]=d2Dpp(p,x0)
% this function calculates the distance of a set of n points from a 2D
% point
% input: p, n x 2 matrix, each row containing the [x y] coordinates of a
%        point
%        x0, vector [x0 y0] of the reference point
% output; d, n x 1 vector of the distances

% determine the distance from the reference point

x0=x0(:)';
d=sqrt(sum((p-repmat(x0,size(p)./size(x0))).^2,2));