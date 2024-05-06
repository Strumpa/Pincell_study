function [iter,evect,eval] = al1eig(A,eps)
% find the fundamental eigenvalue and corresponding eigenvector of
% equation (a-eval)*evect=0 using the power method.
% function [iter,evect,eval] = al1eig(A,eps)
% (c) 2020 Alain Hebert, Polytechnique Montreal
  n=size(A,1);
  evect=ones(n,1);
  iter=0;
  while true
    iter=iter+1;
    if iter > 500, error('al1eig: unable to converge.'), end
    gar=evect;
    evect=A*evect;
    eval=norm(evect);
    evect=evect/eval;
    err1=max(abs(evect)); err2=max(abs(gar-evect));
    if err2 <= err1*eps, break, end
  end