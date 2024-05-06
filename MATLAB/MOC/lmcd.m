function [zmu,wzmu]=lmcd(nmu)
% set the Leonard-McDaniel quadrature base points and weights
% function [zmu,wzmu]=lmcd(nmu)
% (c) 2008 Alain Hebert, Ecole Polytechnique de Montreal
   if nmu == 2
      zmu=[ 1.15511584  3.65419436 ] ; wzmu=[ 0.744970262  0.03816792 ] ;
   elseif nmu == 3
      zmu=[ 1.12178171  2.52822733  10.0188332 ] ;
      wzmu=[ 0.707642138  0.0745818987  0.00175868778 ] ;
   elseif nmu == 4
      zmu=[ 1.10098934  2.15439725  6.02394629  23.7072144 ] ;
      wzmu=[ 0.677962959  0.10134203  0.00533553911  0.000130677523 ] ;
   else
      error('invalid number of base points.')
   end