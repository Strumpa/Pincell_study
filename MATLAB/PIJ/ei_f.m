function f=ei_f(tau0,sig,segment)
% function f=ei_f(tau0,sig,segment)
% (c) 2008 Alain Hebert, Ecole Polytechnique de Montreal
  if sig ~= 0
    f=(akin(3,tau0)-akin(3,tau0+sig*segment))/sig ;
  else
    f=segment*akin(2,tau0) ;
  end