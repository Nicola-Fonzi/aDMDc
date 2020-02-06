function [inputMorphing, bendingOut, alpha] = runMPC(ROMs,dt,simTime,gustTime,deltaAlpha,refState,...
                                                        Velocity,rollActuation,omega,MPC_on)

% Function to run MPC test case
%
% This function runs a simulation of a morphing wing at gust encounter, 
% reducing the first bending mode amplitude using MPC
%
% In this example, the bending mode is assessed using only the aDMDc model,
% whereas in the published manuscript, the full order model is also used. The 
% full order model is made available on request.
%
% MPC code based on the work from Eurika Kaiser
% Sparse Identification of Nonlinear Dynamics for Model Predictive Control in the Low-Data Limit
% https://github.com/eurika-kaiser/SINDY-MPC
%
% Inputs
%
% ROMs:             aDMDc matrices
% dt:               time step, [s]
% simTime:          simulation time, [s]
% Velocity:         flight speed, [m/s]
% rollActuation:    roll actuation input
% omega:            aircraft rotation rates inputs
% gustTime:         gust length, [s]
% deltaAlpha:       gust induced change in angle of attack, [°]
% refState:         reference bending moment amplitude
%
% Outputs
%
% inputMorphing:    morphing wing actuattion input for load alleviation
% bendingOut:       resulting bending mode amplitude
% alpha:            gust induced angle of attack
%

%% define 1-cos gust
Nts = round(simTime/dt);    % number of time steps
Vgust = ones(1,round(gustTime/dt)) + 1/2*(1-cos(2*pi*(0:round(gustTime/dt)-1)/round(gustTime/dt))) - 1;
Vz = zeros(1,Nts+1);
Vz(1:length(Vgust)) = Vgust;
alpha = (atan2(Velocity*tand(deltaAlpha)*Vz,Velocity*ones(1,Nts+1))); % gust induced angle of attack


%% Parameters MPC
options = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
    'MaxIterations',100);

dtC = 0.018;        % controller update delta time

N  = 10;            % Prediction horizon (number of iterations), we assume it
                    % equal to the control horizon      
  
Q = 10;             % State weights
R = 1;              % du weights
Ru = 0;             % u weights

LB = -1*ones(N,1);	% Lower bound of control input
UB = 1*ones(N,1);	% Upper bound of control input
LBdu = -10*dtC;     % Lower bound of control input rate
UBdu = 10*dtC;      % Upper bound of control input rate

LBo = nan;          % Lower bound of output
UBo = nan;          % Upper bound of output

% Reference state, which shall be achieved
xref = refState*ones(Nts+1,1);

% Initialize variables
uopt0    = 0;
uopt     = uopt0.*ones(N,1);
uHistory = zeros(1,Nts+1); uHistory(1)   = uopt(1);
tHistory = zeros(1,Nts+1); tHistory(1)   = 0;

%% run simulation
for i = 1:(Nts+1)

    %% run MPC
    if i == 1 || ~MPC_on
        inputMorphing(i) = 0;   % actuation velocity
        inputMorphingD(i) = 0;  % actuation
    else
        if mod(i*dt,dtC) == 0
            % MPC
            COSTFUN = @(u) ObjectiveFCN_modelsParam(u,N,xref(i),uHistory(:,i),...
                diag(Q),R,Ru,ROMs,inputMorphing(i-1),alpha(i),alpha(i-1));
            CONSFUN = @(u) ConstraintFCN_modelsParam(u,uHistory(:,i),N,LBo,UBo,LBdu,UBdu);
            uopt  = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,CONSFUN,options);
            % smooth morphing actuation inputs
            inputMorphingD(i) = (uopt(1)-inputMorphing(i-1))*dt/dtC;    % actuation velocity
            inputMorphing(i) = inputMorphing(i-1) + inputMorphingD(i);  % actuation
        else
            inputMorphingD(i) = inputMorphingD(i-1);                    % actuation velocity
            inputMorphing(i) = inputMorphing(i-1) + inputMorphingD(i);  % actuation
        end
    end
    uHistory(:,i+1) = uopt(1);
    tHistory(:,i+1) = i*dt;     
    
    
    %% assess bending mode amplitude using the aDMDc model   
    Input_ROMs = [rollActuation; inputMorphing(i); omega; alpha(i)];
    ROMs.State = ROMs.Atil * ROMs.State + ROMs.Btil*ROMs.Input + ROMs.Ftil*Input_ROMs;
    ROMs.Input = Input_ROMs;

    ROMs.Gamma = ROMs.Uhat*ROMs.State + ROMs.Xmean;  % This is the state in physical coordinates
    bendingOut(i) = ROMs.Gamma(610);                 % The position in the vector 
                                                     % is dependent on the
                                                     % position of the
                                                     % bending mode in the
                                                     % base

end

