% Demo Code for One-Third Power Law Principled Protocol
% Power Law Tool Chain EBR v0.9.0
%
% Run 'PowerLawDemo_ExpBrainRes.m' in src 
to generate synthetic ellipse trajectories
% obeying the 1/3 power law, or linear velocites as per
% Maoz et al., 2005 DOI: 10.13140/2.1.4401.8884
% Schaal & Sternad, 2001 DOI: 10.1007/s002210000505
%
% white / pink noise of a range of Standard Deviations may be added
%
% Velocity Gain Factor and Beta Power Law exponents are then calculated via
% choices of differentiation, filtering, curvature calculation, and regression
% as outlined in Exp Brain Res review paper of Fraser et al., 2024
%
% Created May 2024
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk
%
% requires functions subfolder
% - curvatureKinematicEBR.m
% - differentiateKinematicsEBR.m
% - curvatureMengerEBR.m
% - regressDataEBR.m
%
% require req subfolder with
% https://uk.mathworks.com/matlabcentral/fileexchange/34874-interparc
% - interparc/
%             -interparc.m
%






