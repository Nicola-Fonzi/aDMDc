function [DMD_Matrices] = aDMDc(Velocity)
 % function [DMD_Matrices] = aDMDc(Velocity)
 % 
 %
 % The function computes the linearisation State-Space matrices of a nonlinear
 % system based on observations of the system's evolution in time.
 %
 % It is applied to algebraic differential problems.
 %
 % It requires a snapshot matrix for:
 %  1) The evolution of the state
 %  2) The inputs sequence as a function of time
 %
 % The function is applied to LPV systems, where the scheduling parameter
 % may be of any kind. In the context of the code, this is called Velocity.
 %
 % The snapshot matrices must be placed in a folder called "V$PARAMETER$",
 % where the user should substitute $PARAMETER$ with a unique integer number
 % identifying a run with a fixed value of the scheduling parameter. The
 % snapshot matrix containing the evolution of the state will be called:
 %  "Snapshot.mat".
 % While the snapshot matrix containing the set of inputs provided to the
 % system will be called:
 %  "Snapshot_input.mat".
 %
 % The present script will be placed in the directory containing all the
 % folders defining different runs, with different scheduling parameter.
 %
 % The result is saved in a structure and also saved in the same folder 
 % where all the snapshots are saved. The struct is called DMD_Matrices.mat 
 % and contains Atil, Btil, Ftil and Uhat.
 % x^{n+1} = Atil*x^{n} + Btil*u^{n} + Ftil*u^{n+1}

 
 
% Load the snapshot matrix.
% The matrix is evolving in time column-wise
SnapshotLoad = load(strcat(pwd,filesep,'V',num2str(Velocity),filesep,'Snapshot.mat'),'Snapshot');
Snapshot = SnapshotLoad.Snapshot;
clear SnapshotLoad

X = Snapshot(:,1:end-1);
Xp = Snapshot(:,2:end);

% center matrices
Xmean = mean(X,2);
X = X - Xmean;
Xp = Xp - Xmean;

% Same as before for the input snapshots
SnapshotLoad = load(strcat(pwd,filesep,'V',num2str(Velocity),filesep,'Snapshot_input.mat'));
Snapshot_input = SnapshotLoad.Snapshot_input;
clear SnapshotLoad
UPS = Snapshot_input(:,1:end-1);
UPSp = Snapshot_input(:,2:end);

% We can now stack the matrices together and perform the SVD. This will be
% truncated based on the automatic truncation explained in D. L. Donoho 
% and M. Gavish, "The Optimal Hard Threshold for Singular Values is
% 4/sqrt(3)". The obatined SVD will be used to approximate the state-space
% matrices of the system. Note: they would have the same dimension as the
% state, in order to obtain a ROM, we project those matrices on the PODs
% obained with an SVD of the state only --> 2nd SVD in this script.
XModified = [X;UPS;UPSp];
[U,Sig,V] = svd(XModified,'econ');
figure
semilogy(Sig,'bo')
title('Singular values in descreasing order of importance of the SVD of the full Snapshot matrix')
thresh = optimal_SVHT_coef(size(Sig,2)/size(Sig,1),Sig);
rtil = length(find(diag(Sig)>thresh));
Util = U(:,1:rtil);
Sigtil = Sig(1:rtil,1:rtil);
Vtil = V(:,1:rtil);

[U,Sig,V] = svd(Xp,'econ');
figure
semilogy(Sig,'bo')
title('Singular values in descreasing order of importance of the SVD of the matrix Xprime')
thresh = optimal_SVHT_coef(size(Sig,2)/size(Sig,1),Sig);
r = length(find(diag(Sig)>thresh));
Uhat = U(:,1:r);
Sighat = Sig(1:r,1:r);
Vhat = V(:,1:r);

% Here we extract the matrices.
n = size(X,1);
q = size(UPS,1);
U_1 = Util(1:n,:);
U_2 = Util(n+1:n+q,:);
U_3 = Util(n+q+1:n+q+q,:);
Atil = Uhat'*(Xp)*Vtil*inv(Sigtil)*U_1'*Uhat;
Btil = Uhat'*(Xp)*Vtil*inv(Sigtil)*U_2';
Ftil = Uhat'*(Xp)*Vtil*inv(Sigtil)*U_3';

% Here we extract the eigenvalues for plotting purpose.
[DMtil,EIGtil] = eig(Atil);

figure
DiagonalEig=diag(EIGtil);
plot(real(DiagonalEig),imag(DiagonalEig),'gd')
title('Eigenvalues of the A matrix of the system')
hold on, grid on
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
for i=1:length(diag(EIGtil))
    if norm([real(DiagonalEig(i)),imag(DiagonalEig(i))])>1
        disp('Unstable system at Velocity:')
        disp(paramFSI.inputCreate.Velocity)
        plot(real(DiagonalEig(i)),imag(DiagonalEig(i)),'ro')
    end
    break
end

% Construct the struct for immediate usage
DMD_Matrices.Atil = Atil;
DMD_Matrices.Btil = Btil;
DMD_Matrices.Ftil = Ftil;
DMD_Matrices.Uhat = Uhat;
DMD_Matrices.Xmean = Xmean;

% We just save everything is a structure.
save(strcat(pwd,filesep,'V',num2str(Velocity),filesep,'DMD_Matrices.mat'), 'Atil' , 'Btil' , 'Ftil' , 'Uhat' , 'Xmean');
end