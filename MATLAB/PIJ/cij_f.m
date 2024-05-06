% function f=cij_f(tau0,sig,segment)
% (c) 2008 Alain Hebert, Ecole Polytechnique de Montreal
function f=cij_f(tau0,sigi,sigj,segmenti,segmentj)
    if sigi ~= 0 && sigj ~= 0
        f=(akin(3,tau0)-akin(3,tau0+sigi*segmenti)-akin(3,tau0+sigj*segmentj)+ ...
            akin(3,tau0+sigi*segmenti+sigj*segmentj))/(sigi*sigj) ;
    elseif sigi == 0 && sigj ~= 0
        f=(akin(2,tau0)-akin(2,tau0+sigj*segmentj))*segmenti/sigj ;
    elseif sigi ~= 0 && sigj == 0
        f=(akin(2,tau0)-akin(2,tau0+sigi*segmenti))*segmentj/sigi ;
    else
        f=akin(1,tau0)*segmenti*segmentj ;
    end
