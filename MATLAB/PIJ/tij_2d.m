function tij=tij_2d(track,sigt)
% integration of the collision, escape and transmission probabilities
% in unstructured finite 2D geometry.
% function tij=tij_2d(track,sigt)
% (c) 2009 Alain Hebert, Ecole Polytechnique de Montreal
  indpos=@(i,j) max(i,j).*(max(i,j)-1)./2+min(i,j) ;
  nsurf=track(1) ; nreg=track(2) ; k=5+track(1)+track(2)+2*track(3) ;
  tij=zeros(1,(nreg+nsurf)*(nreg+nsurf+1)/2) ;
  for itrk=1:track(4)
    isurf=track(k+2) ; jsurf=track(k+3) ; wei=track(k+4) ; km=track(k+5) ;
    kgar=k+5 ; k=k+5+km ; irs=isurf ; seg1=0. ; sig1=0. ;
    for ixi=1:km
      irt=irs ; irs=track(kgar+ixi) ; seg2=track(k+ixi) ; sig2=sigt(irs) ;
      irs=irs+nsurf ; iij=indpos(irs,irs) ;
      tij(iij)=tij(iij)+2.0*wei*di_f(sig2,seg2) ; tau0=0. ;
      for ixj=ixi:km
        jrs=track(kgar+ixj) ; seg3=track(k+ixj) ; sig3=sigt(jrs) ;
        jrs=jrs+nsurf ; iij=indpos(irt,jrs) ;
        if irt <= nsurf
          tij(iij)=tij(iij)+wei*ei_f(tau0,sig3,seg3) ;
        else
          wi3=cij_f(tau0,sig1,sig3,seg1,seg3) ;
          if jrs == irt, wi3=2.0*wi3; end
          tij(iij)=tij(iij)+wei*wi3 ;
        end
        tau0=tau0+seg3*sig3 ;
      end
      iij=indpos(irt,jsurf) ;
      if irt <= nsurf
        wi3=akin(3,tau0) ;
        if isurf == jsurf, wi3=2.0*wi3; end
        tij(iij)=tij(iij)+wei*wi3 ;
      else
        tij(iij)=tij(iij)+wei*ei_f(tau0,sig1,seg1) ;
      end
      seg1=seg2 ; sig1=sig2 ;
    end
    iij=indpos(irs,jsurf) ; tij(iij)=tij(iij)+wei*ei_f(0.0,sig1,seg1) ;
    k=k+km ;
  end
  tij(:)=tij(:).*track(5)^2 ;