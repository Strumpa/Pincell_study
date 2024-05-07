function b=mcgsces(nom,h,xst)
% calculate coefficients of a track for the characteristics integration.
% step-characteristic scheme with exact exponential calls.
% function b=mcgsces(nom,h,xst)
% (c) 2008 Alain Hebert, Ecole Polytechnique de Montreal
km=length(h) ; b=zeros(2,km) ;
for i=1:km
    nomi=nom(i) ;
    if xst(nomi) == 0.
        b(1,i)=h(i) ;
        b(2,i)=0.5*h(i)^2 ;
    else
        b(1,i)=(1-exp(-h(i)*xst(nomi)))/xst(nomi) ;
        b(2,i)=(h(i)-b(1,i))/xst(nomi) ;
    end
end
