function [c, ceq] = ConstraintFCN_modelsParam(u,uold,N,LBo,UBo,LBdu,UBdu)
%% Constraint function of nonlinear MPC for F8 system
%
% Inputs:
%   u:      optimization variable, from time k to time k+N-1 
%   x:      current state at time k
%   Ts:     controller sample time
%   N:      prediction horizon length
%   uold:   latest applied control input
%   LBo:    Lower bound of output x
%   UBo:    Upper bound of output x
%   LBdu:   Lower bound for input difference uk - uk-1
%   UBdu:   Upper bound for input difference uk - uk-1
%   p:      Parameters for model
%   select_model: Selects model future-state prediction
%
% Output:
%   c:      inequality constraints applied across prediction horizon
%   ceq:    equality constraints  (empty)
%

%% Inequality constraints calculation
c = zeros(N,1);
% Apply N population size constraints across prediction horizon, from time
% k+1 to k+N

% max actuator speed range 
duk = abs(u(1)-uold);
for ct=1:N

    c(ct) = duk-UBdu;
    
    if ct<N
        duk = abs(u(ct+1)-u(ct));
    end
end

%% No equality constraints
ceq = [];
