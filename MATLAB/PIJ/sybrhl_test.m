function pij=sybrhl_test(track,sigt,pij)
% Stamm'ler normalisation of collision, escape and transmission
% probabilities in unstructured finite 2D geometry.
% function pij=sybrhl(track,sigt,pij)
% (c) 2008 Alain Hebert, Ecole Polytechnique de Montreal
%----
%  define anonymous function indpos
%----
  indpos=@(i,j) max(i,j).*(max(i,j)-1)./2+min(i,j) ;
%
  epscon=1.0e-6 ; nitmax=20 ; cptlb=3 ; cptac=3 ;
  nsurf=track(1) ; nreg=track(2) ; nun=nsurf+nreg ; weig=zeros(nun,3) ;
  weig(:,2:3)=0.5 ; g=track(6:5+nsurf+nreg) ; g(1:nsurf)=g(1:nsurf)./4. ;
  g(nsurf+1:nun)=g(nsurf+1:nun).*sigt ; pij_w=pij ;
  disp(g) ;
  disp(weig) ;
  %----
  %  reconstruct real probabilities
  %----
  iij=nsurf*(nsurf+1)/2 ;
  for ir=1:nreg
    disp("iij is :") ;
    disp(iij) ;
    disp("iij+1 :") ;
    disp(iij+1) ;
    disp("iij+nsurf is :") ;
    disp(iij+nsurf);
    disp("values between iij+1 and iij+nsurf") ;
    disp(pij_w(iij+1:iij+nsurf)) ;
    disp(pij_w(iij+1:iij+nsurf).*sigt(ir) );
    disp("ir is :") ;
    disp(ir) ;
    pij_w(iij+1:iij+nsurf)=pij_w(iij+1:iij+nsurf).*sigt(ir) ;
    %iij = iij+nsurf ;
    iij=iij+nsurf+ir ;
  end
  iij=nsurf*(nsurf+1)/2 ;
  for jr=1:nreg
    iij=iij+nsurf ;
    disp("iij is :") ;
    disp(iij) ;
    disp("iij+1 :") ;
    disp(iij+1);
    disp(iij+jr) ;
    disp("values between iij+1 and iij+jr") ;
    disp(pij_w(iij+1:iij+jr)) ;
    disp(pij_w(iij+1:iij+jr).*(sigt(1:jr)*sigt(jr)));
    disp(jr) ;
    disp(sigt(1:jr)*sigt(jr)) ;
    disp(sigt(1:jr)) ;
    disp(sigt(jr)) ;
    pij_w(iij+1:iij+jr)=pij_w(iij+1:iij+jr).*(sigt(1:jr)*sigt(jr)) ;
    iij=iij+jr ;
  end
  %----
  %  main iteration loop
  %----
  nit=0 ; totcon=9999. ;
  while totcon >= epscon
    nit=nit+1 ;
    if nit > nitmax
      error('weights not converged.')
    end
    for ir=1:nun
      disp("ir is :") ;
      disp(ir) ;
      disp(indpos(ir,ir)) ;
      wfspad=g(ir)+pij_w(indpos(ir,ir))*weig(ir,3) ;
      wfspad=wfspad-sum(weig(1:nun,3).*pij_w(indpos(ir,1:nun))') ;
      wfsp=pij_w(indpos(ir,ir))+sum(pij_w(indpos(ir,1:nun))) ;
      if wfsp ~= 0
        weig(ir,3)=wfspad/wfsp ;
      end
    end
    %----
    %  acceleration by residual minimization
    %----
    zmu=1.0 ;
    if mod(nit-1,cptac+cptlb) >= cptac
      nom=sum((weig(:,2)-weig(:,1)).*(weig(:,3)-2.*weig(:,2)+weig(:,1))) ;
      denom=sum((weig(:,3)-2.*weig(:,2)+weig(:,1)).^2) ;
      if denom ~= 0.0
        zmu=-nom/denom ;
      end
      if (zmu > 5.0) || (zmu < 0.0)
        zmu=1.0 ;
      end
      weig(1:nun,3)=weig(1:nun,2)+zmu*(weig(1:nun,3)-weig(1:nun,2)) ;
      weig(1:nun,2)=weig(1:nun,1)+zmu*(weig(1:nun,2)-weig(1:nun,1)) ;
    end
    %----
    %  calculations of square distance between 2 iterations and updating
    %  of the solution
    %----
    totcon=max(abs(weig(:,3)-weig(:,2))./weig(:,3)) ;
    weig(1:nun,1)=weig(1:nun,2) ; weig(1:nun,2)=weig(1:nun,3) ;
  end
  %----
  %  renormalize "pij" symmetric matrix
  %----
  iprb=0 ;
  for ir=1:nun
    disp("in renormalize at ir =") ;
    disp(ir) ;
    disp("weights are :") ;
    disp((weig(ir,1)+weig(1:ir,1)')) ;
    pij(iprb+1:iprb+ir)=pij(iprb+1:iprb+ir).*(weig(ir,1)+weig(1:ir,1)') ;
    iprb=iprb+ir ;
  end