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
% 1) Generate tracking file :
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

% 3) extract collision, escape and transmission probabilities using 
% indpos function which associates Tij entry with upper diag
% entries of full T matrix

indpos=@(i,j) max(i,j).*(max(i,j)-1)./2+min(i,j) ;
T_matrix = zeros(nsurf+nvol,nsurf+nvol) ;

for i=1:nsurf+nvol
    for j=1:nsurf+nvol
        if j>=i
            T_matrix(i,j) = Tij(indpos(i,j)) ;
        end
    end
end

for indi=2:(nsurf+nvol)
    for indj=1:(indi-1)
        T_matrix(indi,indj) = T_matrix(indj,indi) ;
    end
end

%check_symmetric = issymmetric(T_matrix) ;

% decompose into sub-blocks. 

t_SS = T_matrix(1:nsurf, 1:nsurf) ;
t_Sv = T_matrix(1:nsurf, nsurf+1:nsurf+nvol) ;
t_vS = T_matrix(nsurf+1:nsurf+nvol, 1:nsurf) ;
t_ij = T_matrix(nsurf+1:nsurf+nvol, nsurf+1:nsurf+nvol) ;


% recover probabilities / reduced probabilities using eq 3.339
P_SS = zeros(nsurf,nsurf) ;
for alpha=1:nsurf
    for beta=1:nsurf
        P_SS(alpha,beta) = t_SS(alpha,beta)*4/surfaces(alpha) ;
    end
end

%pss = sybpss(tracks, sig_tot) ; % compare with pss given by sybpss 
% sybpss recovered from "Applied Reactor Physics" Appendix A. 
% ---> They're the same.

p_ij = zeros(nvol,nvol) ;

for i=1:nvol
    p_ij(i,:) = t_ij(i,:)/volumes(i) ;
end
disp(p_ij) ;

P_vS = zeros(nvol,nsurf);
for i=1:nvol
    P_vS(i,:) = t_vS(i,:)/volumes(i) ;
end

p_Sv = zeros(nsurf,nvol) ;
for alpha=1:nsurf
    p_Sv(alpha,:) = t_Sv(alpha,:)*4/surfaces(alpha) ;
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

PSS_tilde = albe.*inv(eye(nsurf,nsurf) -albe*P_SS) ; % eq 3.355

Pvv_tilde = p_ij + P_vS*PSS_tilde*p_Sv ; % eq 3.354

% 6) Compute the scattering reduced probability matrix W

W = (eye(nvol,nvol)-Pvv_tilde*S0)\Pvv_tilde ;

% 7) Compute 
[iter,evect,eval] = al1eig(W*Qfiss,10^-8);
Keff=eval;
disp("Keff = ");
disp(Keff);


