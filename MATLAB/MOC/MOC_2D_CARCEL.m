% Question 2, Pre-Doctoral Exam : 
% Method of Characteristics solution to 2D Pincell neutron flux.
% Author : R. Guasch, combining and adapting scripts written by A. Hebert.
% available at https://https://moodle.polymtl.ca/course/view.php?id=1233
side = sqrt(4.9) ; %cm
nangle = 14 ; %number of angles
ngauss = 4 ; %number of gauss quadrature points
nmu = 4 ; %polar angle quadrature order 2, 3, 4 to test to compare with DRAGON5

errtol = 1e-8 ;
maxit = 1000000 ;
beta = 1 ;

Vol_i = [0.4, 0.7, 0.4, 1.3, 2.1] ; %2D volumes : cm^2
sig_tot = [0.2, 0.0, 0.5, 0.0, 0.3]; % total macroscopic cross sections : cm^-1
sig_scattering = [0.05, 0.0, 0.05, 0.0, 0.05] ; % scattering macroscopic cross sections : cm^-1
nu_sig_f = [1.4, 0.0, 0.0, 0.0, 0.0] ; % neutron production cross section = nu (avg number of neutrons per fission) times Sigma_f, fission macroscopic xs.
q = [0.0, 0.0, 0.0, 0.0, 1.4, 0.0, 0.0, 0.0, 0.0]' ;
nsurf = 4 ;
nvol = 5 ;
phi = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]' ;

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
% 2) compute Pii : self-collision factors
pii = mcgpii(tracks, sig_tot, nmu) ;

% 3) iterate over Phi until convergence
[Phi, error, iter,Keff ] = free(phi,q,"mcgsis",errtol,maxit,tracks,sig_tot,sig_scattering,pii,nmu,beta) ;
