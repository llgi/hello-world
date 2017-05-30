function [err,res,x,y,z] = fmagcalibfull(X, S, r)


    x0 = X(1);
    y0 = X(2);
    z0 = X(3);
    sx = X(4);
    sy = X(5);
    sz = X(6);
   
    x  = S(:,1);
    y  = S(:,2);
    z  = S(:,3);
  
    x = sx*(x - x0);
    y = sy*(y - y0);
    z = sz*(z - z0);
  
%   Error = difference between norm sphere radius and the sphere radius calculated 
    err = (r.^2 - (x.^2 + y.^2 + z.^2));
    res=sum(err);
  end
  
