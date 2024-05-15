% Question 1, Pre-Doctoral Exam : 
% Method of Collision Probabilities solution to 2D Pincell neutron flux.
% Author : R. Guasch, combining and adapting scripts written by A. Hébert.
% available at https://https://moodle.polymtl.ca/course/view.php?id=1233

side = sqrt(4.9) ; %cm
surfaces = [side, side, side, side] ;
beta = 1 ; % albedo for isotropic boundary conditions.

nangle = 14 ; %number of angles
ngauss = 4 ; %number of gauss quadrature points
Vol = [0.4, 0.7, 0.4, 1.3, 2.1] ; %2D volumes : cm^2
sig_tot = [0.2, 0.0, 0.5, 0.0, 0.3]; % total macroscopic cross sections : cm^-1
sig_scattering = [0.05, 0.0, 0.05, 0.0, 0.05] ; % scattering macroscopic cross sections : cm^-1
nu_sig_f = [1.4, 0.0, 0.0, 0.0, 0.0] ; % neutron production cross section = nu (avg number of neutrons per fission) times Sigma_f, fission macroscopic xs.

nsurf = 4 ;
nvol = 5 ;

S0=diag(sig_scattering);
Qfiss=diag(nu_sig_f);
% Computing combined volumes in order to get radii more easily.
Vol_combined = zeros(1,size(Vol,2)) ;
Vol_combined(1) = Vol(1) ;
for i=2:5
   Vol_combined(i) = Vol_combined(i-1)+Vol(i) ;
end

% In a CARCEL, there are 1 less radii than volumes as the last volume is
% inscribed in the square boundary, but r>rmax = radii(-1) (last element in radii array).
radii = zeros(1,size(Vol,2)-1);
radii(1) = sqrt(Vol(1)/pi) ;
for i=2:size(Vol,2)-1
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

% 3) extract collision, escape and transmission probabilities using 
% Eq 3.339

% Initialize matrices
% Transmission components
t_ab = zeros(nsurf, nsurf);  % S_a*(P_SaSb) / 4, eq. 3.339
P_SS = zeros(nsurf, nsurf); % = t_ab * 4/S_a, here S_a = side for all a since it's a square pincell.

% Entrance components
%t_ib = zeros(nvol, nsurf); % Vi*P_iSb, eq 3.339
%P_vS = zeros(nvol, nsurf); % = t_ib/Vi
t_bi = zeros(nsurf, nvol) ;
p_Sv = zeros(nsurf, nvol) ;

% collision components
t_ij = zeros(nvol, nvol) ; % Vi*p_ij, eq 3.339 ! Reduced p_ij here !
p_ij = zeros(nvol, nvol); % p_ij = t_ij/Vi


% Compute some lengths and indices to facilitate extraction from T_ij
% vector
n_SS = (nsurf*nsurf - (nsurf-1)*nsurf/2) ; % number of t_ab terms calculated by tij_2d
n_IJ = (nvol*nvol - (nvol-1)*nvol/2) ; % number of t_ij terms calculated by tij_2d
n_IJ_start = size(Tij,2) - n_IJ + 1 ;

% Store elements in P_SS upper diagonal including diagonal
diag_indices_surf = eye(nsurf);
upper_diag_indices_surf = triu(true(nsurf), 1);
t_ab(diag_indices_surf | upper_diag_indices_surf) = Tij(1:n_SS) ;

% extracting upper diagonal terms from Tij (normalized) vector
P_SS(diag_indices_surf | upper_diag_indices_surf) = Tij(1:n_SS)*4/side ; 
for i=2:nsurf
    for j=1:(i-1)
        P_SS(i,j) = P_SS(j,i) * surfaces(i) / surfaces(j) ; % calculating full transmission probbility matrix.
    end
end
pss = sybpss(tracks, sig_tot) ;

% Store elements in P_vS
t_bi (:) = Tij(n_SS+1:n_IJ_start-1);
%P_vS = t_ib/Vi
for alpha=1:nsurf
    p_Sv(alpha,:) = t_bi(alpha,:)*4/surfaces(alpha) ;
end

% Store elements in P_ij upper diagonal including diagonal
diag_indices_vol = eye(nvol);
upper_diag_indices_vol = triu(true(nvol), 1);
t_ij(diag_indices_vol | upper_diag_indices_vol) = Tij(n_IJ_start:end);
disp(t_ij)
%p_ij = t_ij/Vi
for i=1:nvol
    p_ij(i,:) = t_ij(i,:) / Vol(i) ;
end
% use Vi*pij = Vj*pji !
for i=2:nvol
    for j = 1:(i-1)
        p_ij(i,j) = p_ij(j,i)*(Vol(i)/Vol(j));
    end
end

% Get P_vS from p_Sv :
P_vS = zeros(nvol,nsurf) ;
for i=1:nvol
    for j=1:nsurf
        P_vS(i,j) = p_Sv(j,i)*surfaces(j)/(4*Vol(i)) ;
    end
end

% 4) Normalize probabilities using conservation equations :
sums_for_alpha = zeros(1,nsurf) ;
sums_for_i = zeros(1,nvol); % create arrays to store normalization coefficients
% 4.1) compute normalization factors
for i=1:nvol
    for beta=1:nsurf
        sums_for_i(i) = sums_for_i(i) + P_vS(i,beta) ;
    end
    for j=1:nvol
        sums_for_i(i) = sums_for_i(i) + p_ij(i,j)*sig_tot(j) ;
    end
end

for alpha=1:nsurf
    for beta=1:nsurf
        sums_for_alpha(alpha) = sums_for_alpha(alpha) + P_SS(alpha,beta) ;
    end
    for j=1:nvol
        sums_for_alpha(alpha) = sums_for_alpha(alpha) + p_Sv(alpha,j)*sig_tot(j) ;
    end
end

% 4.2) now renormalize by deviding by the total sum for each alpha or i.
for alpha=1:nsurf
    P_SS(alpha,:) = P_SS(alpha,:)/sums_for_alpha(alpha) ;
    p_Sv(alpha,:) = p_Sv(alpha,:)/sums_for_alpha(alpha) ;
end
for i=1:nvol
    p_ij(i,:) =  p_ij(i,:)/sums_for_i(i) ;
    P_vS(i,:) =  P_vS(i,:)/sums_for_i(i) ;
end

% 4.3)  check normalization
sumI = zeros(1, nvol) ;
sumAlpha = zeros(1,nsurf) ;
for i=1:nvol
    for beta=1:nsurf
        sumI(i) = sumI(i)+P_vS(i,beta) ;
    end
    for j=1:nvol
        sumI(i) = sumI(i)+p_ij(i,j)*sig_tot(j) ;
    end
end
for alpha=1:nsurf
    for beta=1:nsurf
        sumAlpha(alpha) = sumAlpha(alpha) + P_SS(alpha,beta) ;
    end
    for j=1:nvol
        sumAlpha(alpha) = sumAlpha(alpha) + p_Sv(alpha,j)*sig_tot(j) ;
    end
end

% 5) Compute the closed reduced collision probability matrix :
% Use eq. 3.350 and 3.351

A = eye(nsurf)*beta ;
I = eye(nsurf) ;
PSS_tilde = A*inv(I-P_SS*A) ; % eq 3.355
% PSS_tilde = A/(I-P_SS*A) for optimized inverse

P_Sv = p_Sv.*sig_tot ;

Pvv_tilde = p_ij + P_vS*PSS_tilde*P_Sv ; % eq 3.354

% 6) Compute the scattering reduced probability matrix W
%W=(eye(size(pwb,1))-pwb*S0)^-1*pwb*Qfiss;

W = (eye(nvol)-Pvv_tilde*S0)*Pvv_tilde*Qfiss ;

% 7) Compute 
[iter,evect,eval] = al1eig(W,10^-8);
Keff=eval;
disp("Keff = ");
disp(Keff);


