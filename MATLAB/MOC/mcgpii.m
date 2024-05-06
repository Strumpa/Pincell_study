function pii=mcgpii_dd(track,sigt,lmax,nmu,ielem)
% compute the pii components for source isolation with the MOC.
% P0/P1 linear/quadratic diamond differencing scheme.
% function pii=mcgpii_dd(track,sigt,lmax,nmu,ielem)
% (c) 2010 Alain Hebert, Ecole Polytechnique de Montreal
  if ielem>2, error('ielem must be equal to 1 or 2.'), end
  nsurf=track(1) ; nreg=track(2) ; nangle=track(3) ; nbtr=track(4) ;
  nfunl=(lmax+1)*(lmax+2)/2 ;
%----
%   scattering anisotropy treatment
%----
  eta_2d=track(5+nsurf+nreg+1:5+nsurf+nreg+nangle) ;
  xi_2d=track(6+nsurf+nreg+nangle:5+nsurf+nreg+2*nangle) ;
  [zmu,wzmu]=lmcd(nmu) ;
  trharp=zeros(nfunl,nmu*nangle) ; trharm=trharp ;
  for iangle=1:nangle
    for imu=1:nmu
      eta_3d=eta_2d(iangle)/zmu(imu) ;
      xi_3d=xi_2d(iangle)/zmu(imu) ;
      mu_3d=sqrt(1.-1./zmu(imu)^2) ;
      iof=0 ;
      for il=0:lmax
        for im=-il:il
          if mod(il+im,2) == 1; continue; end
          iof=iof+1;
          trharp(iof,(iangle-1)*nmu+imu)=pnsh(il,im,mu_3d,eta_3d,xi_3d) ;
          trharm(iof,(iangle-1)*nmu+imu)=pnsh(il,im,-mu_3d,-eta_3d,-xi_3d) ;
        end
      end
    end
  end
%----
%   pii calculation
%----
  k=5+track(1)+track(2)+2*track(3) ;
  volnum=zeros(nreg,1) ; pii=zeros(nfunl,nfunl,nreg) ;
  if ielem==1
    for iline=1:nbtr
      iangle=track(k+1) ; weitf=track(k+4) ; km=track(k+5) ; kgar=k+5 ; k=k+5+km ;
      nom=track(kgar+1:kgar+km) ; htf=track(k+1:k+km) ; h=zeros(1,km) ;
      for imu=1:nmu
        ww=weitf*wzmu(imu) ; h(:)=htf(:).*zmu(imu) ;
        angmod=trharp(:,(iangle-1)*nmu+imu)*trharp(:,(iangle-1)*nmu+imu)' ;
        for i=1:km
          nomi=nom(i) ;
          volnum(nomi)=volnum(nomi)+2*h(i)*ww ;
          denom=h(i)*sigt(nomi)+2 ;
          pii(:,:,nomi)=pii(:,:,nomi)+(ww*h(i)^2/denom)*angmod ;
        end
        angmod=trharm(:,(iangle-1)*nmu+imu)*trharm(:,(iangle-1)*nmu+imu)' ;
        for i=km:-1:1
          nomi=nom(i) ;
          denom=h(i)*sigt(nomi)+2 ;
          pii(:,:,nomi)=pii(:,:,nomi)+(ww*h(i)^2/denom)*angmod ;
        end
      end
      k=k+km ;
    end
  elseif ielem==2
    for iline=1:nbtr
      iangle=track(k+1) ; weitf=track(k+4) ; km=track(k+5) ; kgar=k+5 ; k=k+5+km ;
      nom=track(kgar+1:kgar+km) ; htf=track(k+1:k+km) ; h=zeros(1,km) ;
      for imu=1:nmu
        ww=weitf*wzmu(imu) ; h(:)=htf(:).*zmu(imu) ;
        angmod=trharp(:,(iangle-1)*nmu+imu)*trharp(:,(iangle-1)*nmu+imu)' ;
        for i=1:km
          nomi=nom(i) ;
          volnum(nomi)=volnum(nomi)+2*h(i)*ww ;
          q11=h(i)*sigt(nomi) ; q12=2*sqrt(3) ;
          q21=q12 ; q22=-h(i)*sigt(nomi)-6 ;
          denom=q11*q22-q12*q21 ;
          pii(:,:,nomi)=pii(:,:,nomi)+(ww*h(i)^2/denom)*q22*angmod ;
        end
        angmod=trharm(:,(iangle-1)*nmu+imu)*trharm(:,(iangle-1)*nmu+imu)' ;
        for i=km:-1:1
          nomi=nom(i) ;
          q11=h(i)*sigt(nomi) ; q12=2*sqrt(3) ;
          q21=q12 ; q22=-h(i)*sigt(nomi)-6 ;
          denom=q11*q22-q12*q21 ;
          pii(:,:,nomi)=pii(:,:,nomi)+(ww*h(i)^2/denom)*q22*angmod ;
        end
      end
      k=k+km ;
    end
  end
  for ii1=1:nreg
    pii(:,:,ii1)=pii(:,:,ii1)/volnum(ii1) ;
  end