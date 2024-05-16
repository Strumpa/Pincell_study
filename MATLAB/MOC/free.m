function [x, error, iter, Keff] = free(x, b, atv, errtol, maxit, varargin)
% free-iteration linear equation solver
% function [x, error, total_iters] = free(x0, b, atv, errtol, maxit, varargin)
% input parameters:
% x      initial estimate
% b      right hand side
% atv    character name of a matrix-vector product routine returning x+M(b-Ax)
%        when x is input. The format for atv is "function x = atv(x,b,p1,p2,...)"
%        where p1 and p2 are optional parameters
% errtol relative residual reduction factor
% maxit  maximum number of iterations
% varargin optional parameters (p1,p2,...) for atv
% output parameters:
% x      solution of the linear system
% error  vector of residual norms for the history of the iteration
% iter   number of iterations
% (c) 2007 Alain Hebert, Polytechnique Montreal
  errtol=errtol*norm(b);
  rho=Inf; error=[ ];
  iter=0;
  eval_ini = 1 ;
  Fission_integral = b'*x ;
  while((rho > errtol) && (iter < maxit))
    iter=iter+1;
    r=feval(atv,x,b,varargin{:})-x;
    rho=norm(r); error=[error;rho];
    x=x+r;
    new_Fission_integral = b'*x  ;
    new_eval = eval_ini*(Fission_integral/new_Fission_integral) ;
    eval_ini = new_eval ;
    Fission_integral = new_Fission_integral ;
  end
  Keff = 1/eval_ini ;
