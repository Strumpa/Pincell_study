% Question 1, Pre-Doctoral Exam : 
% Method of Collision Probabilities solution to 2D Pincell neutron flux.
% Author : R. Guasch, combining and adapting scripts written by A. Hébert.
% available at https://https://moodle.polymtl.ca/course/view.php?id=1233

side = sqrt(4.9) ; %cm
surfaces = [side, side, side, side] ;
beta = 1 ; % albedo for isotropic boundary conditions.

nangle = 14 ; %number of angles
ngauss = 4 ; %number of gauss quadrature points
Vol_i = [0.4, 0.7, 0.4, 1.3, 2.1] ; %2D volumes : cm^2
sig_tot = [0.2, 0.0, 0.5, 0.0, 0.3]; % total macroscopic cross sections : cm^-1
sig_scattering = [0.05, 0.0, 0.05, 0.0, 0.05] ; % scattering macroscopic cross sections : cm^-1
nu_sig_f = [1.4, 0.0, 0.0, 0.0, 0.0] ; % neutron production cross section = nu (avg number of neutrons per fission) times Sigma_f, fission macroscopic xs.

nsurf = 4 ;
nvol = 5 ;

S0=diag(sig_scattering);
Qfiss=diag(nu_sig_f);
% Computing combined volumes in order to get radii more easily.
Vol_combined = zeros(1,size(Vol_i,2)) ;
Vol_combined(1) = Vol_i(1) ;
for i=2:5
   Vol_combined(i) = Vol_combined(i-1)+Vol_i(i) ;
end

% In a CARCEL, there are 1 less radii than volumes as the last volume is
% inscribed in the square boundary, but r>rmax = radii(-1) (last element in radii array).
radii = zeros(1,size(Vol_i,2)-1);
radii(1) = sqrt(Vol_i(1)/pi) ;
for i=2:size(Vol_i,2)-1
    radii(i) = sqrt(Vol_combined(i)/pi());
end
% 1) Generate tracking file :
tracks = sybt2d(side,radii,nangle,ngauss) ; % calling sybt2d to generate tracking file
% sybt2d.m retrieved from the ENE6101 moodle course page, presented in
% Appendix A.4 of "Applied Reactor Physics" - A. Hébert.

% 2) Compute the symmetric T matrix :
Tij = tij_2d(tracks,sig_tot) ; % calling tij_2d to compute the T matrix for a finite 2D unstructured geometry.
% This might need to be revisited as Hebert mentions that cij_f, di_f, ei_f
% were not written to deal with cross-sections approaching 0. but deal with
% sig = 0, so might not need to change them.

% build symmetric matrix T from Tij lower diagonal coefficients
% Store elements in t_ij matrix upper diagonal including diagonal
%t_ij = zeros(nsurf+nvol, nsurf+nvol) ;
%diag_indices = eye(nsurf+nvol) ;
%lower_diag_indices = tril(true(nsurf+nvol), -1) ;
%t_ij(diag_indices | lower_diag_indices) = Tij ;

% compute the upper diagonal symetric
%for j=2:(nsurf+nvol)
%    for i=1:(j-1)
%        t_ij(i,j) = t_ij(j,i) ;
%    end
%end
% 3) normalize with T sybrhl.m to get T_tilde

T_tilde=sybrhl(tracks,sig_tot,Tij) ;


t_ij = zeros(nsurf+nvol, nsurf+nvol) ;
diag_indices = eye(nsurf+nvol) ;
lower_diag_indices = tril(true(nsurf+nvol), -1) ;
t_ij(diag_indices | lower_diag_indices) = T_tilde ;
% compute the upper diagonal symetric
for j=2:(nsurf+nvol)
    for i=1:(j-1)
        t_ij(i,j) = t_ij(j,i) ;
    end
end
% recover P_SS, P_vS, pij and p_Sv from T_tilde entries --> these should be
% normalized probabilities

t_ss = t_ij(1:nsurf, 1:nsurf) ;
t_vS = t_ij(nsurf+1:nsurf+nvol, 1:nsurf) ;
t_Sv = t_ij(1:nsurf, nsurf+1:nsurf+nvol) ;
t_tilde_ij = t_ij(nsurf+1:nsurf+nvol, nsurf+1:nsurf+nvol) ;

P_ss = zeros(nsurf,nsurf) ;
for alpha=1:nsurf
    for beta=1:nsurf
        P_ss(alpha,beta) = t_ss(alpha,beta)*4/surfaces(alpha) ;
    end
end

p_Sv = zeros(nsurf, nvol) ;
for alpha=1:nsurf
    for j=1:nvol
        p_Sv(alpha,j) = t_Sv(alpha,j)*4/surfaces(alpha) ;
    end
end

P_vS = zeros(nvol,nsurf) ;
for i=1:nvol
    for alpha=1:nsurf
        P_vS(i,alpha) = t_vS(i,alpha)/Vol_i(i) ;
    end
end

p_ij = zeros(nvol,nvol) ;
for i=1:nvol
    for j=1:nvol
        p_ij(i,j) = t_tilde_ij(i,j)/Vol_i(i);
    end
end

% 3) bis check normalization
sumI = zeros(1, nvol) ;
sumAlpha = zeros(1,nsurf) ;
for i=1:nvol
    for beta=1:nsurf
        sumI(i) = sumI(i)+P_vS(i,beta) ;
    end
    for j=1:nvol
        sumI(i) = sumI(i)+p_ij(i,j)*sig_tot(j) ;
        disp(sig_tot(j)) ;
    end
end

for alpha=1:nsurf
    for beta=1:nsurf
        sumAlpha(alpha) = sumAlpha(alpha) + P_ss(alpha,beta) ;
    end
    for j=1:nvol
        sumAlpha(alpha) = sumAlpha(alpha) + p_Sv(alpha,j)*sig_tot(j) ;
    end
end

% 4) Compute the closed reduced collision probability matrix :
% Use eq. 3.350 and 3.351

% 5) Compute the scattering reduced probability matrix W

%W = (eye(nvol)-Pvv_tilde*S0)*Pvv_tilde*Qfiss ;

% 6) Compute 
%[iter,evect,eval] = al1eig(W,10^-8);
%Keff=eval;
%disp("Keff = ");
%disp(Keff);
