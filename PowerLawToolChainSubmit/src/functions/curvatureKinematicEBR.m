function [curvature] = curvatureKinematicEBR(velocityX, velocityY,accelerationX, accelerationY)
% Curvature formula from Viviani & Stucchi, 1992 DOI 10.1037//0096-1523.18.3.603
% adapted from Python implementation of Matic & Gomez-Marin, 2019 DOI : 10.1016/j.jneumeth.2019.108398
% https://github.com/adam-matic/KinematicCognition-Analysis
% curvature in terms of only acceleration and velocity
%% demo code for Exp Brain Res review paper of Fraser et al., 2024
% Created May 2024
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk
%
%% inputs
% velocity and acceleration in x and y
%% outputs 
% curvature


denominator = ( ( velocityX.^2 + velocityY.^2 ) .^0.5 ) .^3;

numerator = abs( (velocityX .* accelerationY) - (velocityY .* accelerationX) );

curvature = numerator ./ denominator;



