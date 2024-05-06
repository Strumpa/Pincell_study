% function f=di_f(tau0,sig,segment)
% (c) 2008 Alain Hebert, Ecole Polytechnique de Montreal
function f=di_f(sig,segment)
    if sig ~= 0
        f=segment/sig-(akin(3,0)-akin(3,sig*segment))/sig^2 ;
    else
        f=pi*segment^2/4 ;
    end
    