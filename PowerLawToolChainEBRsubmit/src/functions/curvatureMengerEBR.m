function [curvature] = curvatureMengerEBR(tripletXY)
% Adapted from Maria Botero's menger_curvature.m
% https://github.com/ucam-ceb-como/HRTEMFringeMapping
% A detailed description can be found in the papers: 
% Botero et al., 2016 doi:10.1016/j.carbon.2015.09.077
% Botero et al., 2019 doi:10.1016/j.carbon.2018.09.063
%
%% demo code for Exp Brain Res review paper of Fraser et al., 2024
% Created May 2024
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk
%
%% input is tripletXY is 3 euclidean points
% of 2 rows, x 3 columns
%% output is curvature at mid poit

% Curvature of three points; is the curvature of the unique circle that 
% passes through these three points. 
% The Menger curvature of such circle is 4 times it's area -  divided by the 
% product of the length of the three sides   
% curvature, c = 1/R = (4*A) / (x-y)(y-z)(z-x)
% where x y and z are non colinear euclidean points

tripleX = tripletXY(1,:);
tripleY = tripletXY(2,:);

% calculate the length of the Triangle Sides
sideLength(1) = sqrt((tripleX(2)-tripleX(1))^2 + (tripleY(2)-tripleY(1))^2);
sideLength(2) = sqrt((tripleX(3)-tripleX(1))^2 + (tripleY(3)-tripleY(1))^2);
sideLength(3) = sqrt((tripleX(3)-tripleX(2))^2 + (tripleY(3)-tripleY(2))^2);

% calculate the area of the Triangle
A = 1/2*((tripleX(1)*(tripleY(2)-tripleY(3))+ tripleX(2)*(tripleY(3)-tripleY(1)) + tripleX(3)*(tripleY(1)-tripleY(2))));

% calculate Curvature
curvature = ( 4 * A ) / ( sideLength(1) * sideLength(2) * sideLength(3) );




