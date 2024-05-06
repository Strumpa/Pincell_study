function [A B]=mcgaca(track,sigt,nmu,beta)
% matrix assembly for the second-order ACA operator.
% function [A B]=mcgaca(track,sigt,nmu,beta)
% (c) 2008 Alain Hebert, Ecole Polytechnique de Montreal
  nsurf=track(1) ; nreg=track(2) ; nbtr=track(4) ;
  [zmu,wzmu]=lmcd(nmu) ;
%----
%  matrix assembly
%----
  k=5+track(1)+track(2)+2*track(3) ;
  A=sparse(nreg+nsurf,nreg+nsurf) ; B=A ;
  for iline=1:nbtr
    isurf=track(k+2) ; jsurf=track(k+3) ; weitf=track(k+4) ; km=track(k+5) ;
    kgar=k+5 ; k=k+5+km ;
    nom=track(kgar+1:kgar+km) ; htf=track(k+1:k+km) ; h=zeros(1,km) ;
    for imu=1:nmu
      ww=weitf*wzmu(imu) ; h(:)=htf(:).*zmu(imu) ;
      dinv=zeros(1,km) ; b=zeros(1,km) ;
      for i=1:km
        taud=h(i)*sigt(nom(i)) ;
        if taud ~= 0.
          alpha=2.0/(1.0-exp(-taud))-2.0/taud-1.0 ;
          dinv(i)=taud/(2.0+taud*alpha) ;
          b(i)=0.5*h(i)*(dinv(i)-alpha) ;
        end
      end
      ind2p=nsurf+nom(1) ;
      A(isurf,isurf)=A(isurf,isurf)+ww*(1.0+beta+dinv(1)*(1.0-beta))/2. ;
      A(isurf,ind2p)=A(isurf,ind2p)-ww*(1.0-b(1)*sigt(nom(1))) ;
      B(isurf,isurf)=B(isurf,isurf)-ww*(1.0-dinv(1))/2. ;
      B(isurf,ind2p)=B(isurf,ind2p)+ww*b(1) ;
      ind2m=nsurf+nom(km) ;
      A(jsurf,jsurf)=A(jsurf,jsurf)+ww*(1.0+beta+dinv(km)*(1.0-beta))/2. ;
      A(jsurf,ind2m)=A(jsurf,ind2m)-ww*(1.0-b(km)*sigt(nom(km))) ;
      B(jsurf,jsurf)=B(jsurf,jsurf)-ww*(1.0-dinv(km))/2. ;
      B(jsurf,ind2m)=B(jsurf,ind2m)+ww*b(km) ;
      if km == 1
        nomi=nom(1) ; ind1=nsurf+nomi ;
        A(ind1,isurf)=A(ind1,isurf)+0.5*ww*(1.0-beta) ;
        A(ind1,ind1)=A(ind1,ind1)+ww*h(1)*sigt(nomi) ;
        A(ind1,jsurf)=A(ind1,jsurf)+0.5*ww*(1.0-beta) ;
        B(ind1,isurf)=B(ind1,isurf)+0.5*ww ;
        B(ind1,ind1)=B(ind1,ind1)+ww*h(1) ;
        B(ind1,jsurf)=B(ind1,jsurf)+0.5*ww ;
      else
        for i=1:km
          nomi=nom(i) ; ind1=nsurf+nomi ;
          if i == 1
            noml=nom(i+1) ; dhp=dinv(i+1)+dinv(i) ; ind2p=nsurf+noml ;
            A(ind1,isurf)=A(ind1,isurf)+0.5*ww*(1.0-beta) ;
            A(ind1,ind1)=A(ind1,ind1)+ww*((1.0-b(i)*sigt(nomi))/dhp+h(i)*sigt(nomi)) ;
            A(ind1,ind2p)=A(ind1,ind2p)-ww*(1.0-b(i+1)*sigt(noml))/dhp ;
            B(ind1,isurf)=B(ind1,isurf)+0.5*ww ;
            B(ind1,ind1)=B(ind1,ind1)-ww*(b(i)/dhp-h(i)) ;
            B(ind1,ind2p)=B(ind1,ind2p)+ww*b(i+1)/dhp ;
          elseif i == km
            nomj=nom(i-1) ; dhm=dinv(i)+dinv(i-1) ; ind2m=nsurf+nomj ;
            A(ind1,ind2m)=A(ind1,ind2m)-ww*(1.0-b(i-1)*sigt(nomj))/dhm ;
            A(ind1,ind1)=A(ind1,ind1)+ww*((1.0-b(i)*sigt(nomi))/dhm+h(i)*sigt(nomi)) ;
            A(ind1,jsurf)=A(ind1,jsurf)+0.5*ww*(1.0-beta) ;
            B(ind1,ind2m)=B(ind1,ind2m)+ww*b(i-1)/dhm ;
            B(ind1,ind1)=B(ind1,ind1)-ww*(b(i)/dhm-h(i)) ;
            B(ind1,jsurf)=B(ind1,jsurf)+0.5*ww ;
          else
            nomj=nom(i-1) ; dhm=dinv(i)+dinv(i-1) ; ind2m=nsurf+nomj ;
            noml=nom(i+1) ; dhp=dinv(i+1)+dinv(i) ; ind2p=nsurf+noml ;
            A(ind1,ind2m)=A(ind1,ind2m)-ww*(1.0-b(i-1)*sigt(nomj))/dhm ;
            A(ind1,ind1)=A(ind1,ind1)+ww*((1.0/dhm+1.0/dhp)*(1.0-b(i)*sigt(nomi))+h(i)*sigt(nomi)) ;
            A(ind1,ind2p)=A(ind1,ind2p)-ww*(1.0-b(i+1)*sigt(noml))/dhp ;
            B(ind1,ind2m)=B(ind1,ind2m)+ww*b(i-1)/dhm ;
            B(ind1,ind1)=B(ind1,ind1)-ww*((1.0/dhm+1.0/dhp)*b(i)-h(i)) ;
            B(ind1,ind2p)=B(ind1,ind2p)+ww*b(i+1)/dhp ;
          end
        end
      end
    end
    k=k+km ;
  end