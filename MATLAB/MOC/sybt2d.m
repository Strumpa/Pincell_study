function track=sybt2d(a,rad,nangle,ngauss)
% produce a Gauss-Jacobi tracking in 2D square pincell geometry
% function track=sybt2d(a,rad,nangle,ngauss)
% (c) 2009 Alain Hebert, Ecole Polytechnique de Montreal
  nreg=1+size(rad,2) ; radd = [0. rad] ; na2=2*nangle ;
  if ngauss == 1
    alp= .6666666667 ; pwr= .5 ; zx= 0. ; wx= 2. ;
  elseif ngauss == 2
    alp=[ .3550510257,.8449489743 ] ; pwr=[ .1819586183,.3180413817 ] ;
    zx=[ -.577350259,.577350259 ] ; wx=[ 1.,1. ] ;
  elseif ngauss == 3
    alp=[ .2123405382,.5905331356,.9114120405 ] ;
    pwr=[ .0698269799,.2292411064,.2009319137 ] ;
    zx=[ -.774596691,0.,.774596691 ];
    wx=[ .555555556,.888888889,.555555556 ] ;
  elseif ngauss == 4
    alp=[ .1397598643,.4164095676,.7231569864,.9428958039 ] ;
    pwr=[ .0311809710,.1298475476,.2034645680,.1355069134 ] ;
    zx=[ -.861136317,-.339981049,.339981049,.861136317 ] ;
    wx=[ .347854853,.652145147,.652145147,.347854853 ] ;
  elseif ngauss == 5
    alp=[ .0985350858,.3045357266,.5620251898,.8019865821,.9601901429 ] ;
    pwr=[ .0157479145,.0739088701,.1463869871,.1671746381,.0967815902 ] ;
    zx=[ -.906179845,-.538469315,0.,.538469315,.906179845 ] ;
    wx=[ .236926883,.478628665,.568888843,.478628665,.236926883 ] ;
  elseif ngauss == 6
    alp=[ .0730543287,.2307661380,.4413284812,.6630153097,.8519214003,.9706835728 ] ;
    pwr=[ .0087383018,.0439551656,.0986611509,.1407925538,.1355424972,.0723103307 ] ;
    zx=[ -.932469487,-.661209404,-.238619193,.238619193,.661209404,.932469487 ] ;
    wx=[ .171324492,.360761583,.467913926,.467913926,.360761583,.171324492 ] ;
  else
    error('invalid number of Gauss-Jacobi points.')
  end
  if 2.0*radd(nreg) > a
    error('a radius is greater than half a side.')
  end
  track_w=zeros(1,9+nreg+4*na2*(2+ngauss*nreg*(5+2*(2*nreg-1)))) ;
  kstart=9+nreg+8*na2 ; track_w(1:3)=[4, nreg, 4*na2] ; zn1=0. ;
  ao2=a/2. ; wa=2./real(na2) ; track_w(6:9)=a ; vol=a*a ;
  for jjj=nreg:-1:1
    r2=pi*radd(jjj)^2 ; track_w(9+jjj)=vol-r2 ; vol=r2 ;
  end
%----
%  track generation
%----
  k=kstart ;
  for ia=1:na2
    za=(2.0*real(ia)-1.)/real(na2)-1. ; phi=0.25*pi*(za+1.) ;
    zn1=zn1+sin(phi)*wa ; si=sin(phi) ; co=cos(phi) ; ta=si/co ;
    track_w(9+nreg+ia)=si ; track_w(9+nreg+4*na2+ia)=co ;
    track_w(9+nreg+na2+ia)=co ; track_w(9+nreg+5*na2+ia)=-si ;
    if phi <= 0.25*pi
      track_w(9+nreg+2*na2+ia)=co ; track_w(9+nreg+6*na2+ia)=-si ;
      track_w(9+nreg+3*na2+ia)=si ; track_w(9+nreg+7*na2+ia)=co ;
      jsu=4 ; x1=0. ; xlim=a ; dlim=ao2*co+(ao2-xlim)*si ;
    else
      track_w(9+nreg+2*na2+ia)=co ; track_w(9+nreg+6*na2+ia)=si ;
      track_w(9+nreg+3*na2+ia)=si ; track_w(9+nreg+7*na2+ia)=-co ;
      jsu=2 ; x1=a/ta ; xlim=0.5*(a+a/ta) ; dlim=0. ;
    end
    for k0=nreg:-1:1
      km=nreg-k0+1 ; x2=min(xlim,xlim-(radd(k0)-dlim)/si) ;
      if ((x1 < xlim) && (phi <= 0.25*pi)) || ((x1 < x2) && (phi > 0.25*pi))
        l3=k ; vap=zeros(1,nreg) ;
        for ix=1:ngauss
          if k0 == nreg
            s=0.5*(x2-x1)*si*wx(ix) ;
            x=x1+0.5*(x2-x1)*(1.0+zx(ix)) ;
          else
%           Flurig change of variable.
            s=2.*(x2-x1)*si*pwr(ix) ;
            x=x1+(x2-x1)*alp(ix)^2 ;
          end
          track_w(k+1:k+5)=[ia, 1, jsu, s*wa/4., 2*km-1] ;
          track_w(k+6:k+2*km+4)=abs(km-1:-1:1-km)+1+nreg-km ;
          k=k+5+(2*km-1) ;
          c=ao2*si-(ao2-x)*co ; d=(ao2*co+(ao2-x)*si)^2 ; sumtrk=0. ;
          for kk=nreg:-1:k0+1
            corde=sqrt(rad(kk-1)^2-d) ; del=c-corde ; sumtrk=sumtrk+del ;
            track_w(k+nreg-kk+1)=del ; vap(kk)=vap(kk)+del*s ;
            c=corde ;
          end
          if km ~= 1
            del=2.0*corde ; track_w(k+km)=del ; vap(k0)=vap(k0)+del*s ;
            sumtrk=sumtrk+del+sum(track_w(k+km-1:-1:k+2)) ;
            track_w(k+km+1:k+2*km-2)=track_w(k+km-1:-1:k+2) ;
            vap(k0+1:k0+km-2)=vap(k0+1:k0+km-2)+track_w(k+km-1:-1:k+2).*s ;
          end
          k=k+2*km-1 ;
          if phi <= 0.25*pi
            del=x/co-sumtrk ;
          else
            del=a/si-sumtrk ;
          end
          track_w(k)=del ; vap(nreg)=vap(nreg)+del*s ;
        end
%----
%  volume normalization
%----
        if k0 < nreg
          dlim1=ao2*co+(ao2-x2)*si ; dlim2=ao2*co+(ao2-x1)*si ;
          vw1=0. ; sumvap=0. ;
          for i=k0:nreg-1
            sumvap=sumvap+vap(i) ; rw=rad(i) ;
            vex1=rw*rw*acos(dlim1/rw)-dlim1*sqrt(rw*rw-dlim1*dlim1) ;
            if rw > dlim2
              vex1=vex1-(rw*rw*acos(dlim2/rw)-dlim2*sqrt(rw*rw-dlim2*dlim2)) ;
            end
            vap(i)=(vex1-vw1)/vap(i) ; vw1=vex1 ;
          end
          vex1=0.5*(a*si-(a-x1-x2)*co)*(x2-x1)*si ;
          if phi <= 0.25*pi
            vex2=0.5*ta*(x2*x2-x1*x1)-vex1 ;
          else
            vex2=(x2-x1)*a-vex1 ;
          end
          vex1=(vex1-0.5*vw1)/(vex1-0.5*sumvap) ;
          vex2=(vex2-0.5*vw1)/(vex2-0.5*sumvap) ;
          for ix=1:ngauss
            l3=l3+5 ; km=(track_w(l3)+1)/2 ; l3=l3+2*km-1 ;
            track_w(l3+km)=track_w(l3+km)*vap(k0) ;
            fact=vap(k0+1:k0+km-2) ;
            track_w(l3+km-1:-1:l3+2)=track_w(l3+km-1:-1:l3+2).*fact ;
            track_w(l3+km+1:l3+2*km-2)=track_w(l3+km+1:l3+2*km-2).*fact ;
            track_w(l3+1)=track_w(l3+1)*vex1 ; 
            track_w(l3+2*km-1)=track_w(l3+2*km-1)*vex2 ;
            l3=l3+2*km-1 ;
          end
        end
        track_w(4)=track_w(4)+ngauss ; x1=x2 ;
      end
    end
  end
  track_w(5)=1./sqrt(0.25*pi*zn1) ; kend=k ;
%----
%  apply symmetries
%----
  track_w(kend+1:2*kend-kstart)=track_w(kstart+1:kend) ;
  for itrk=1:track_w(4)
    track_w(k+1)=na2+track_w(k+1) ; nseg=track_w(k+5) ;
    if track_w(k+3) == 2
      track_w(k+2)=4 ; track_w(k+3)=3 ;
    elseif track_w(k+3) == 4 ;
      track_w(k+3)=3 ;
    end
    track_w(k+6+nseg:k+5+2*nseg)=track_w(k+5+2*nseg:-1:k+6+nseg) ;
    k=k+5+2*track_w(k+5) ;
  end
  track_w(4)=2*track_w(4) ; kend=k ;
  track_w(kend+1:2*kend-kstart)=track_w(kstart+1:kend) ;
  for itrk=1:track_w(4)
    track_w(k+1)=2*na2+track_w(k+1) ;
    if track_w(k+3) == 4
      track_w(k+2)=4 ; track_w(k+3)=2 ;
    elseif track_w(k+3) == 2
      track_w(k+2)=3 ; track_w(k+3)=4 ;
    elseif track_w(k+2) == 1 && track_w(k+3) == 3
      track_w(k+2)=3 ; track_w(k+3)=2 ;
    elseif track_w(k+2) == 4 && track_w(k+3) == 3
      track_w(k+2)=1 ; track_w(k+3)=2 ;
    end
    k=k+5+2*track_w(k+5) ;
  end
  track_w(4)=2*track_w(4) ;
%
  if k > 9+nreg+4*na2*(2+ngauss*nreg*(5+2*(2*nreg-1)))
    error('tracking overflow.')
  end
  track=track_w(1:k) ;