function pii=mcgpii(track,sigt,nmu)
% compute the pii components for source isolation with the MOC.
% function pii=mcgpii(track,sigt,nmu)
% (c) 2008 Alain Hebert, Ecole Polytechnique de Montreal
nreg=track(2) ; nbtr=track(4) ;
[zmu,wzmu]=lmcd(nmu) ;
k=5+track(1)+track(2)+2*track(3) ;
volnum=zeros(nreg,1) ; pii=zeros(nreg,1) ;
for iline=1:nbtr
    weitf=track(k+4) ; km=track(k+5) ; kgar=k+5 ; k=k+5+km ;
    nom=track(kgar+1:kgar+km) ; htf=track(k+1:k+km) ; h=zeros(1,km) ;
    for imu=1:nmu
            ww=weitf*wzmu(imu) ; h(:)=htf(:).*zmu(imu) ;
            b=mcgsces(nom,h,sigt) ;
        for i=1:km
            nomi=nom(i) ;
            volnum(nomi)=volnum(nomi)+h(i)*ww ;
            pii(nomi)=pii(nomi)+b(2,i)*ww ;
        end
    end
    k=k+km ;
end
pii=pii./volnum ;