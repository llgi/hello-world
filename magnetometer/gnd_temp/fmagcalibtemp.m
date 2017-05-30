% This function applies the temperature and offset correction to input data
function [err,res,xcorr,ycorr,zcorr] = fmagcalibtemp(X, S, r)

%   Offset, scale and temperature corrections
    x0  = X(1);
    y0  = X(2);
    z0  = X(3);
    sx  = X(4);
    sy  = X(5);
    sz  = X(6);
    stx = X(7);
    sty = X(8);
    stz = X(9);

%   Input measurements
    x    = S(:,1);
    y    = S(:,2);
    z    = S(:,3);
    temp = S(:,4);

%   Offset, scale and temperature correction    
    xCorr = sx*(x.*(1+sxt*temp) + x0);
    yCorr = sy*(y.*(1+syt*temp) + y0);
    zCorr = sz*(z.*(1+szt*temp) + z0);

%   Error    
    err = (r.^2 - (xcorr.^2 + ycorr.^2 + zcorr.^2));
    res=sum(err);

  end
  
