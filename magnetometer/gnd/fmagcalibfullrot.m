function [err,res,x,y,z,scale,offset,a] = fmagcalibfullrot(X, S, r)

%   Initial values
    x0    = X(1);
    y0    = X(2);
    z0    = X(3);
    sx    = X(4);
    sy    = X(5);
    sz    = X(6);
    phi   = X(7);
    theta = X(8);
    psi   = X(9);
 
%   Input measurement data 
    x  = S(:,1);
    y  = S(:,2);
    z  = S(:,3);
%    tx = S(:,4);
%    ty = S(:,5);
%    tz = S(:,6);

    x = sx*(x - x0);
    y = sy*(y - y0);
    z = sz*(z - z0);
 
%   Temperature anf offset correction
%    x = sx*(x*(1+sx*tx) - x0) 
%    y = sy*(y*(1+sy*ty) - y0) 
%    x = sz*(z*(1+sz*tz) - z0)
  
    cphi = cos(phi); 
    sphi = sin(phi); 
    cth  = cos(theta); 
    sth  = sin(theta); 
    cpsi = cos(psi); 
    spsi = sin(psi); 
  
    a = [cpsi*cth  -spsi*cphi+cpsi*sth*sphi  spsi*sphi+cpsi*cphi*sth; 
	       spsi*cth  cpsi*cphi+sphi*sth*spsi   -cpsi*sphi+sth*spsi*cphi; 
	       -sth      cth*sphi                  cth*cphi];

    A = a * [x y z]';
    A = pinv(a) * [sx*A(1,:);sy*A(2,:);sz*A(3,:)];
    x = A(1,:);
    y = A(2,:);
    z = A(3,:);
    
    scale = [sx;sy;sz];
    offset = [x0;y0;z0];
    
    err = r - sqrt(x.^2 + y.^2 + z.^2);
    res=sum(err);
  end
  