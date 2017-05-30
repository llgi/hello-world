function tJ=juliandate(dateNum)
%this function converts a matlab dateNum to the julian date

    [Y M D H MN S]=datevec(dateNum);
    
    tJ = 367.*Y - floor((7*(Y+floor((M+9)/12)))/4) + floor(275*M/9) + D + 1721013.5 ...
         + H/24 + MN/(24*60) + S/(24*3600) - 0.5*sign(100*Y+M-190002.5) + 0.5;
     
%    fprintf('Julian date: %6.6f\n',tJ);
end