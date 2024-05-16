% Question 1, Pre-Doctoral Exam : 
% Method of Collision Probabilities solution to 2D Pincell neutron flux.
% Author : R. Guasch, combining and adapting scripts written by A. Hébert.
% available at https://https://moodle.polymtl.ca/course/view.php?id=1233

side = sqrt(4.9) ; %cm
albe = 1.0 ; % albedo for isotropic boundary conditions.

nangle = 14 ; %number of angles
ngauss = 4 ; %number of gauss quadrature points
Vol = [0.4, 0.7, 0.4, 1.3, 2.1] ; %2D volumes : cm^2
sig_tot = [0.2, 0.0, 0.5, 0.0, 0.3]; % total macroscopic cross sections : cm^-1
sig_scattering = [0.05, 0.0, 0.05, 0.0, 0.05] ; % scattering macroscopic cross sections : cm^-1
nu_sig_f = [1.4, 0.0, 0.0, 0.0, 0.0] ; % neutron production cross section = nu (avg number of neutrons per fission) times Sigma_f, fission macroscopic xs.


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
% 1) Generate tracking :
tracks = sybt2d(side,radii,nangle,ngauss) ; % calling sybt2d to generate tracking file
% sybt2d.m retrieved from the ENE6101 moodle course page, presented in
% Appendix A.4 of "Applied Reactor Physics" - A. Hébert.
nsurf = tracks(1) ;
nvol = tracks(2) ;
surfaces = tracks(6:5+nsurf) ;
volumes = tracks(6+nsurf:5+nsurf+nvol) ;

% 2) Compute the symmetric T matrix :
Tij = tij_2d(tracks,sig_tot) ; % calling tij_2d to compute the Tij compressed vector
% tij_2d.m recovered from "Applied Reactor Physics" moodle page / presented
% in "Applied Reactor Physics" Chapter 3.8.6.


% 3) Normalize using the sybrhl.m script, implemeting the Stamm'ler
% normalization algorithm
Tij=sybrhl(tracks,sig_tot,Tij);

% 4) extract collision, escape and transmission probabilities using 
% indpos function which associates Tij entry with upper diag
% entries of full T matrix
indpos=@(i,j) max(i,j).*(max(i,j)-1)./2+min(i,j) ;
p_matrix = zeros(nsurf+nvol,nsurf+nvol) ;

for i=1:nsurf
    p_matrix(i,1:nvol+nsurf) = Tij(indpos(i,1:nvol+nsurf)).*(4.0/surfaces(i)) ;
end
for i=nsurf+1:nsurf+nvol
    p_matrix(i,1:nvol+nsurf) = Tij(indpos(i,1:nvol+nsurf))/volumes(i-nsurf) ;
end


% decompose into sub-blocks : easier for me to think about it like this! 

p_ij = p_matrix(nsurf+1:nsurf+nvol, nsurf+1:nsurf+nvol) ;
p_vS = p_matrix(nsurf+1:nsurf+nvol, 1:nsurf) ;
p_Sv = p_matrix(1:nsurf, nsurf+1:nsurf+nvol) ;
p_SS = p_matrix(1:nsurf, 1:nsurf) ;



%pss = sybpss(tracks, sig_tot) ; % compare with pss given by sybpss 
% sybpss recovered from "Applied Reactor Physics" Appendix A. 
% They're the same --> confident that so far, so good.


% 5.1) Use equation 3.355 to compute P_SS_tilde
%A = eye(nsurf,nsurf)*albe ;
p_SS_tilde = albe.*inv(eye(nsurf,nsurf)-albe.*p_SS) ;

% 5.2) use equation 3.354 to compute P_vv_tilde, the closed CP matrix.

P_vv_tilde = p_ij + p_vS*p_SS_tilde*p_Sv ;

% 6) use equation 3.261 to compute W, the scattering-reduced CP matrix.

W = inv(eye(nvol) - P_vv_tilde*S0)*P_vv_tilde ;

system = W*Qfiss ;
% 7) Compute 
[iter,evect,eval] = al1eig(system,10^-8);
Keff=eval;
disp("Keff = ");
disp(Keff);
disp("convergence reached in iterations =") ;
disp(iter) ;
