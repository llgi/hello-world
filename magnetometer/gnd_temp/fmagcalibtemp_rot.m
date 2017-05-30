% This function applies the temperature, offset and rotation correction to input data
function [err,res,xcorr,ycorr,zcorr,scale,offset,a] = fmagcalibtemp_rot(X, S, r)

%   Offset, scale, temperature and rotation parameters
    x0    = X(1);
    y0    = X(2);
    z0    = X(3);
    sx    = X(4);
    sy    = X(5);
    sz    = X(6);
    sxt   = X(7);
    syt   = X(8);
    szt   = X(9);
    phi   = X(10);
    theta = X(11);
    psi   = X(12);

%   Input data    
    x    = S(:,1);
    y    = S(:,2);
    z    = S(:,3);
    temp = S(:,4);
 
%   Offset , scale and temperature correction     
    xOff = sx*(x.*(1+sxt*temp) + x0);
    yOff = sy*(y.*(1+syt*temp) + y0);
    zOff = sz*(z.*(1+szt*temp) + z0);

    cphi = cos(phi); 
    sphi = sin(phi); 
    cth  = cos(theta); 
    sth  = sin(theta); 
    cpsi = cos(psi); 
    spsi = sin(psi); 

%   Rotation matrix 
    a = [cpsi*cth  -spsi*cphi+cpsi*sth*sphi  spsi*sphi+cpsi*cphi*sth; 
	 spsi*cth  cpsi*cphi+sphi*sth*spsi   -cpsi*sphi+sth*spsi*cphi; 
	 -sth      cth*sphi                  cth*cphi];

    A = a * [xOff yOff zOff]';
    A = pinv(a) * [sx*A(1,:);sy*A(2,:);sz*A(3,:)];

%   Data corrected
    xcorr = A(1,:);
    ycorr = A(2,:);
    zcorr = A(3,:);
    
    scale = [sx;sy;sz];
    offset = [x0;y0;z0];

%   Error    
    err = r - sqrt(xcorr.^2 + ycorr.^2 + zcorr.^2);
    res=sum(err);

  end
