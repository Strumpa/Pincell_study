function pss=sybpss(track,sigt)
% integration of transmission probabilities. The tracks are computed by sybt2d.
% function pss=sybpss(track,sigt)
% (c) 2008 Alain Hebert, Ecole Polytechnique de Montreal
%--
% define anonymous function indpos
%--
indpos=@(i,j) max(i,j).*(max(i,j)-1)./2+min(i,j) ;
%--
% pss integration
%--
nsurf=track(1) ; surfa=track(6:5+nsurf) ;
k=5+track(1)+track(2)+2*track(3) ; tij=zeros(1,nsurf*(nsurf+1)/2) ;
for itrk=1:track(4)
    isurf=track(k+2) ; jsurf=track(k+3) ; ind=indpos(isurf,jsurf) ;
    z1=track(k+4) ; km=track(k+5) ; kgar=k+5 ; k=k+5+km ;
    pop=sum(sigt(track(kgar+1:kgar+km)).*track(k+km:-1:k+1)) ;
    tij(ind)=tij(ind)+akin(3,pop)*z1 ; k=k+km ;
end
pss=zeros(nsurf,nsurf) ;
for i=1:nsurf
    pss(i,1:nsurf)=tij(indpos(i,1:nsurf)).*(4.*track(5)^2/surfa(i)) ;
end