% Example force calculation script
%
% This script calculates the force for a spherical particle in a plane
% wave beam.  This is a good comparison with the FDTD force estimates.
%
% Part of ot-cpp-fdtd, Copyright 2020 Isaac Lenton
% See LICENSE for details about using/distributing this file.


% Ensure OTT is on Matlab's path

%% Simulation properties

index_medium = 1.0;
index_sphere = 1.59;

wavelength0 = 1064e-9;
wavelength_medium = wavelength0 ./ index_medium;
wavelength_sphere = wavelength0 ./ index_sphere;

%% Setup T-matrix for sphere

radius = 0.3*wavelength_medium;
Tmatrix = ott.TmatrixMie(radius, 'wavelength0', wavelength0, ...
    'index_medium', index_medium, 'index_particle', index_sphere);

%% Setup the beam

beam = ott.BscPlane(0, 0, 'Nmax', Tmatrix.Nmax, ...
  'wavelength0', wavelength0, 'index_medium', index_medium);

%% Calculate forces

force = ott.forcetorque(beam, Tmatrix);
disp(['OTT Force: ' num2str(force.')]);

%% Load FDTD output files, convert forces and display values

% TODO

