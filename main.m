clear all; close all; clc

% Main file that load the snapshot matrices, create the state space
% matrices for the system and finally perform the MPC control
% This is a very simplified example were only one ROM, at a specific
% velocity, is used. The controller will try to minimize the root bending
% moment due to the a vertical gust. At the steady state condition it was
% registered the value of the bending mode amplitude.
% This reference is then used as objective of the optimisation.
% The snapshots were generated using a mean angle of attack of 5°, no
% side-slip, no rotation rates.
% It is important to note that the DMD matrices are already in discrete
% time. Thus, they are only valid at the particular time step size used in
% the simulation.

%% Values required for the current example
Velocity = 60;              % Velocity used in the example
dt = 0.006;                 % time step, [s]
simTime = 3.0;              % simulation time, [s]
rollActuation = 0;          % roll actuation input
omega = [0;0;0];            % aircraft rotation rates inputs, [rad]
gustTime = 0.5;             % gust length, [s]
deltaAlpha = 1.5;           % gust induced change in angle of attack, [°]
refState = 0.7203;          % reference bending moment amplitude

MPC_on = true;              % run simulation with or without MPC, can be used,
                            % for example, to compute the difference
                            % between actuated and non-actuated wing.

%% Obtain state-space matrices
[ROMs] = aDMDc(Velocity);

%% Initialise the state of the reduced order model: V = 60m/s
ROMs.State = zeros(size(ROMs.Atil,1),1);         % initial state of aDMDc
ROMs.Input = zeros(6,1);                         % initial input of aDMDc

%% run MPC case
[inputMorphing, bendingOut, alpha] = runMPC(ROMs,dt,simTime,gustTime,deltaAlpha,...
    refState,Velocity,rollActuation,omega,MPC_on);

%% plot gust induced angle of attack
figure
plot(0:dt:simTime,alpha*180/pi,'b')
grid on
title('angle of attack')

%% plot actuation
figure
plot(0:dt:simTime,inputMorphing,'b')
grid on
title('morphing actuation')

%% plot bending mode
figure
plot(0:dt:simTime,bendingOut/bendingOut(1)-1,'b')
grid on
title('first bending mode')