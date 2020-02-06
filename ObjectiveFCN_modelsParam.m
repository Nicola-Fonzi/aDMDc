function J = ObjectiveFCN_modelsParam(u,N,xref,u0,Q,R,Ru,ROMs,u_old,alphaPlane,alphaPlaneLastIteration)
%% Cost function of nonlinear MPC for HIV system
%
% Inputs:
%   u:      optimization variable, from time k to time k+N-1
%   x:      current state at time k
%   Ts:     controller sample time
%   N:      prediction horizon
%   xref:   state references, constant from time k+1 to k+N
%   u0:     previous controller output at time k-1
%
% Output:
%   J:      objective function cost
%


%% Integrate system
for ii = 1:N
    if ii > 1
        u_old = u(ii-1);
        alphaPlaneIold = alphaPlane;
    else
        alphaPlaneIold = alphaPlaneLastIteration;
    end
    
    u_state = [0,u(ii),0,0,0,alphaPlane]';
    u_old_state = [0,u_old,0,0,0,alphaPlaneIold]';
    
    ROMs.State = ROMs.Atil * ROMs.State + ROMs.Btil*u_old_state + ROMs.Ftil*u_state;
    ROMs.Gamma = ROMs.Uhat*ROMs.State+ ROMs.Xmean;
    
    xk(ii) = ROMs.Gamma(610); % first bending mode
end


%% Cost Calculation
% Set initial plant states, controller output and cost.
uk = u(1);
J = 0;
% Loop through each prediction step.
for ct=1:N
    % Obtain plant state at next prediction step.
    xk1 = xk(:,ct);
    
    % accumulate state tracking cost from x(k+1) to x(k+N).
    J = J + (xk1-xref')'*Q*(xk1-xref');
    % accumulate MV rate of change cost from u(k) to u(k+N-1).
    if ct==1
        J = J + (uk-u0)'*R*(uk-u0) + uk'*Ru*uk;
    else
        J = J + (uk-u(ct-1))'*R*(uk-u(ct-1)) + uk'*Ru*uk;
    end
    % Update uk for the next prediction step.
    if ct<N
        uk = u(ct+1);
    end
end
