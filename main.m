% Main file that load the snapshot matrices, create the state space
% matrices for the system and finally perform the MPC control

% Values required for the current example
StartPoint = 201;
Velocity = 60;

% Obtain state-space matrices
[DMD_Matrices] = aDMDc(Velocity,StartPoint);

% Perform the MPC

%PLACEHOLDER