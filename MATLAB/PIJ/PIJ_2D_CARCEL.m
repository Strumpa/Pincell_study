% Question 1, Pre-Doctoral Exam : 
% Method of Collision Probabilities solution to 2D Pincell neutron flux.
% Author : R. Guasch, combining and adapting scripts written by A. Hébert.
% available at https://https://moodle.polymtl.ca/course/view.php?id=1233

cote = sqrt(4.9) ; %cm
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
tracks = sybt2d(cote,radii,nangle,ngauss) ; % calling sybt2d to generate tracking file
% sybt2d.m retrieved from the ENE6101 moodle course page, presented in
% Appendix A.4 of "Applied Reactor Physics" - A. Hébert.

% 2) Compute the symmetric T matrix :
Tij = tij_2d(tracks,sig_tot) ; % calling tij_2d to compute the T matrix for a finite 2D unstructured geometry.
disp(size(Tij,2))
% This might need to be revisited as Hebert mentions that cij_f, di_f, ei_f
% were not written to deal with cross-sections approaching 0. but deal with
% sig = 0, so might not need to change them.

% 3) Normalize Tij using the Villarino-Stamm'ler method :
T_tilde = sybrhl(tracks,sig_tot,Tij) ;
disp(size(T_tilde,2))

% 3) bis, extract collision, escape and transmission probabilities using 
% Eq 3.339

% Define the array of length 45

% Initialize matrices
P_SS = zeros(nsurf, nsurf);
P_vS = zeros(nvol, nsurf);
P_ij = zeros(nvol, nvol);

n_SS = (nsurf*nsurf - (nsurf-1)*nsurf/2) ;
n_IJ = (nvol*nvol - (nvol-1)*nvol/2) ;
n_IJ_start = size(T_tilde,2) - n_IJ + 1 ;
n_VS_start = n_IJ + 1 ;
disp("nSS=")
disp(n_SS)
disp("nVS start=")
disp(n_VS_start)
disp("nIJ=")
disp(n_IJ)
disp("nIJ start=")
disp(n_IJ_start)
% Store elements in P_SS upper diagonal including diagonal
diag_indices_surf = eye(nsurf);
upper_diag_indices_surf = triu(true(nsurf), 1);
P_SS(diag_indices_surf | upper_diag_indices_surf) = T_tilde(1:n_SS);

% Store elements in P_vS
P_vS(:) = T_tilde(11:30);

% Store elements in P_ij upper diagonal including diagonal
diag_indices_vol = eye(nvol);
upper_diag_indices_vol = triu(true(nvol), 1);
P_ij(diag_indices_vol | upper_diag_indices_vol) = T_tilde(n_IJ_start:end);


% 4) Compute the closed reduced collision probability matrix :
% Use eq. 3.350 and 3.351

%PSS_tilde = A*(I-PSS*A)^-1
%Pvv_tilde = Pvv + PvS*PSS_tilde*PSv

% 5) Compute the scattering reduced probability matrix W
%W=(eye(size(pwb,1))-pwb*S0)^-1*pwb*Qfiss;

% 6) Compute 
%[iter,evect,eval] = al1eig(W,10^-8);
%Keff1=eval;
%disp("Keff1 = ");
%disp(Keff1);
tij=zeros(1,(5+4)*(5+4+1)/2)


