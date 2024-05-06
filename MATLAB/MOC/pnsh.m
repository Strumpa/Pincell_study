function f=pnsh(l,m,zmu,eta,xi)
% return the real spherical harmonics corresponding to a set of
% direction cosines.
% function f=pnsh(l,m,zmu,eta,xi)
% (c) 2009 Alain Hebert, Ecole Polytechnique de Montreal
   test=zmu*zmu+eta*eta+xi*xi;
   if(abs(test-1.0) > 1.0e-5)
      error('pnsh: invalid direction cosines.')
   end
   if((l == 0) && (m == 0))
      f=1.0;
   elseif((l == 1) && (m == -1))
      f=xi;
   elseif((l == 1) && (m == 0))
      f=zmu;
   elseif((l == 1) && (m == 1))
      f=eta;
   elseif((l == 2) && (m == -2))
      f=sqrt(3.0)*eta*xi;
   elseif((l == 2) && (m == -1))
      f=sqrt(3.0)*zmu*xi;
   elseif((l == 2) && (m == 0))
      f=0.5*(3.0*zmu*zmu-1.0);
   elseif((l == 2) && (m == 1))
      f=sqrt(3.0)*zmu*eta;
   elseif((l == 2) && (m == 2))
      f=0.5*sqrt(3.0)*(eta*eta-xi*xi);
   elseif((l == 3) && (m == -3))
      f=sqrt(5./8.)*xi*(3.0*eta*eta-xi*xi);
   elseif((l == 3) && (m == -2))
      f=sqrt(15.0)*eta*xi*zmu;
   elseif((l == 3) && (m == -1))
      f=sqrt(3./8.)*xi*(5.0*zmu*zmu-1.0);
   elseif((l == 3) && (m == 0))
      f=0.5*zmu*(5.0*zmu*zmu-3.0);
   elseif((l == 3) && (m == 1))
      f=sqrt(3./8.)*eta*(5.0*zmu*zmu-1.0);
   elseif((l == 3) && (m == 2))
      f=sqrt(15.0/4.0)*zmu*(eta*eta-xi*xi);
   elseif((l == 3) && (m == 3))
      f=sqrt(5./8.)*eta*(eta*eta-3.0*xi*xi);
   elseif((l == 4) && (m == -4))
      f=0.5*sqrt(35.)*eta*xi*(eta*eta-xi*xi);
   elseif((l == 4) && (m == -3))
      f=0.5*sqrt(0.5*35.)*zmu*xi*(3.*eta*eta-xi*xi);
   elseif((l == 4) && (m == -2))
      f=sqrt(5.)*(21.*zmu*zmu-3.)*eta*xi/6.;
   elseif((l == 4) && (m == -1))
      f=0.5*sqrt(2.5)*zmu*xi*(7.*zmu*zmu-3.);
   elseif((l == 4) && (m == 0))
      f=(35.*zmu^4-30.*zmu*zmu+3.)/8.;
   elseif((l == 4) && (m == 1))
      f=0.5*sqrt(2.5)*zmu*eta*(7.*zmu*zmu-3.);
   elseif((l == 4) && (m == 2))
      f=sqrt(5.)*(21.*zmu*zmu-3.)*(eta*eta-xi*xi)/12.;
   elseif((l == 4) && (m == 3))
      f=0.5*sqrt(0.5*35.)*zmu*eta*(eta*eta-3.*xi*xi);
   elseif((l == 4) && (m == 4))
      f=sqrt(35.)*(eta^4-6.*(eta*xi)^2+xi^4)/8.;
   elseif((l == 5) && (m == -5))
      f=21.*xi*(5.*eta^4-10.*(eta*xi)^2+xi^4)/(8.*sqrt(14.));
   elseif((l == 5) && (m == -4))
      f=0.5*105.*zmu*eta*xi*(eta*eta-xi*xi)/sqrt(35.);
   elseif((l == 5) && (m == -3))
      f=35.*(9*zmu*zmu-1.)*xi*(3.*eta*eta-xi*xi)/(8.*sqrt(70.));
   elseif((l == 5) && (m == -2))
      f=0.5*sqrt(105.)*zmu*(3.*zmu*zmu-1.)*eta*xi;
   elseif((l == 5) && (m == -1))
      f=sqrt(15.)*xi*(21.*zmu^4-14.*zmu*zmu+1.)/8.;
   elseif((l == 5) && (m == 0))
      f=zmu*(63.*zmu^4-70.*zmu*zmu+15.)/8.;
   elseif((l == 5) && (m == 1))
      f=sqrt(15.)*eta*(21.*zmu^4-14.*zmu*zmu+1.)/8.;
   elseif((l == 5) && (m == 2))
      f=0.25*sqrt(105.)*zmu*(3.*zmu*zmu-1.)*(eta*eta-xi*xi);
   elseif((l == 5) && (m == 3))
      f=35.*(9*zmu*zmu-1.)*eta*(eta*eta-3.*xi*xi)/(8.*sqrt(70.));
   elseif((l == 5) && (m == 4))
      f=105.*zmu*(eta^4-6.*(eta*xi)^2+xi^4)/(8.*sqrt(35.));
   elseif((l == 5) && (m == 5))
      f=21.*eta*(eta^4-10.*(eta*xi)^2+5.*xi^4)/(8.*sqrt(14.));
   else
      error('pnsh: legendre order not available.')
   end
