function phi=mcgsis(phi,q,track,sigt,sigw,pii,nmu,beta)
% single MOC scattering iteration with source isolation
% function phi=mcgsis(phi,q,track,sigt,sigw,pii,nmu,beta)
% (c) 2009 Alain Hebert, Ecole Polytechnique de Montreal
  nsurf=track(1) ; nreg=track(2) ; nbtr=track(4) ;
  s=q+phi.*[beta.*ones(1,nsurf) sigw]' ;
  disp([beta.*ones(1,nsurf) sigw]') ;       
  [zmu,wzmu]=lmcd(nmu) ;
  %----
  %  flux calculation
  %----
 k=5+track(1)+track(2)+2*track(3) ;
volsur=zeros(nsurf+nreg,1) ; phi=zeros(nsurf+nreg,1) ;
for iline=1:nbtr
    isurf=track(k+2) ; jsurf=track(k+3) ; weitf=track(k+4) ; km=track(k+5) ;
    kgar=k+5 ; k=k+5+km ;
    nom=track(kgar+1:kgar+km) ; htf=track(k+1:k+km) ; h=zeros(1,km) ;
    for imu=1:nmu
        ww=weitf*wzmu(imu) ; h(:)=htf(:).*zmu(imu) ;
        b=mcgsces(nom,h,sigt) ;
        rp=s(isurf) ; rm=s(jsurf) ;
        for i=1:km
            nomi=nom(i) ;
            volsur(nsurf+nomi)=volsur(nsurf+nomi)+2*h(i)*ww ;
            phi(nsurf+nomi)=phi(nsurf+nomi)+b(1,i)*rp*ww ;
            rp=rp+b(1,i)*(s(nsurf+nomi)-sigt(nomi)*rp) ;
        end
        for i=km:-1:1
            nomi=nom(i) ;
            phi(nsurf+nomi)=phi(nsurf+nomi)+b(1,i)*rm*ww ;
            rm=rm+b(1,i)*(s(nsurf+nomi)-sigt(nomi)*rm) ;
        end
        phi(jsurf)=phi(jsurf)+rp*ww ; phi(isurf)=phi(isurf)+rm*ww ;
        volsur(jsurf)=volsur(jsurf)+ww ; volsur(isurf)=volsur(isurf)+ww ;
    end
    k=k+km ;
end
phi=phi./volsur ;
unk=phi(nsurf+1:nsurf+nreg)+pii.*q(nsurf+1:nsurf+nreg) ;
phi(nsurf+1:nsurf+nreg)=unk./(ones(nreg,1)-pii.*sigw') ; 