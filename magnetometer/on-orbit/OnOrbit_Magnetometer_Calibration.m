function OnOrbit_Magnetometer_Calibration( data_file, path_lib )
% this function performs an on-orbit magnetometer calibration based on the
% following inputs:
%   1. Three axis Magnetometer Data w/ timestamps
%       -units of milliGauss [mG]
%       -should be 1000+ points dispersed over the attitude sphere
%       -ideally collected over 24 hours or less (so one TLE is valid)
%       -ideally collected during a quiet period of geomagnetic activity (Kp < 5)
%           real-time: http://www.spaceweatherlive.com/en/auroral-activity/kp
%           historic: ftp://ftp.ngdc.noaa.gov/STP/GEOMAGNETIC_DATA/INDICES/KP_AP/
%           more info: http://www.swpc.noaa.gov/noaa-scales-explanation
%   2. Path of Util library. If it is not inserted, generic path will be used
%   3. TLE corresponding to the satellite position during measurement (as a parameter)
%
%   The TLE is used to generate the magnitude of the B-field over
%   time, which is used as the "truth" for a least-squares fit to the
%   collected magnetometer data

%  Check the number of inputs
   if( nargin < 1) 
     error('Not input data file introduced')  
   end

%% Setup
    if ( nargin == 2 ) 
      addpath( genpath( path_lib ) ); 
    else
      addpath( genpath( './../Utils' ) );
    end
    
%% User Inputs:
%   Start date
    auxDate    = input( 'Input start date in the correct format ("11 Oct 2015 15:50:00") [UTC]: ' );
    begDateNum = datenum(auxDate);   
%   End date
    auxDate    = input( 'Input end date in the correct format ("12 Oct 2015 15:50:00") [UTC]: ' );
    endDateNum = datenum( auxDate );
    
%   Check that the time interval is consistent
    if ( begDateNum > endDateNum )
      error( 'The end data is earlier than start date' )
    end  
    
%   TLE corresponding to the input date
    tle.line1 = input( 'Introduce the first line of the TLE: ');
    tle.line2 = input( 'Introduce the second line of the TLE:');
    
    % *** ground calibration constants also needed - see below *** %

%   On board or external magnetometer:
    magNum = input( 'Select on board magnetometer (1) or external (2): ' );

%   Magnetic field already calculated:
    calcBfield = input( 'Calculate the magnetic field? (YES=1 NO=0): ' ); 
   
    if ( calcBfield == 1 ) 
%     Save magnetic field
      saveBfield = input( 'Save the magnetic field calculated (YES=1 NO=0): ' );
    else
%     Do not save it, as it is already calculated and saved in a file
      saveBfield = 0;
    end

%% Load Data

    dat = csvread(data_file,2,0);

%   Select data to keep only the measurements that belong to the time interval    
    uT = dat(:,1);
    goodA = find(uT>=datenum2unix(begDateNum));
    goodB = find(uT<=datenum2unix(endDateNum));
    good = intersect(goodA,goodB);

%   If there are not any measurements: error message and stop function
    if ( isempty(good) )
      error('The input data not belong to the time interval selected.')
    end    
    
    timeDn = unix2datenum(uT(good));
    timeUt = uT(good);
    
    if magNum==1
    
%       board magnetometer
        mxGC = dat(good,4);
        myGC = dat(good,5);
        mzGC = dat(good,6);
        mValid = dat(good,7);
%       Temperature from temperature sensor next to internal gyro
        temp = dat(good,3);
        
        fileDir = './BoardMag/';
        fprintf('Board Magnetometer\n');
        
%       ground calibration parameters
        mOffsetXYZ = [-17.83 49.451 30.852];
        mScaleXYZ  = [0.95692 0.99022 1.0674];
        tScaleXYZ  = [0.000 0.000 0.000];
        rotMatVec  = [0.89801   -0.36724 -0.2423 ...
                      0.4391    0.71337  0.54616 ...
                      -0.027723 -0.59685 0.80187];

     elseif magNum==2
%      external magnetometer
        mxGC = dat(good,8);
        myGC = dat(good,9);
        mzGC = dat(good,10);
        mValid = dat(good,11);
        temp = dat(good,2);
        fileDir = './ExtMag/';
        fprintf('External Magnetometer\n');
        %ground calibration parameters
        mOffsetXYZ = [-22.314 21.691 12.578];
         mScaleXYZ = [1.1002 1.1180 1.0863];
         tScaleXYZ = [0.000 0.000 0.000];
         rotMatVec = [0.749700  -0.152340 0.644000 ...
                      -0.060309 0.953370  0.295720 ...
                      -0.659020 -0.260540 0.705560]

    else
        error('magNum must be 1 or 2');
    end
   
%   Array with the norma of the good measurents    
    normGC = sqrt(mxGC.^2+myGC.^2+mzGC.^2);
    
%   Make the directory for plots
    [~,~]=mkdir(fileDir);

%   Get rid of bad data
    badA=find(abs(normGC-nanmean(normGC))>nanstd(normGC)*5);
    badB=find(mValid~=1);
    bad=union(badA,badB);
    mxGC(bad)=NaN;
    myGC(bad)=NaN;
    mzGC(bad)=NaN;
    normGC(bad)=NaN;
    temp(bad)=NaN;
    
%%  Remove all NaNs from data
    better = find(~isnan(normGC));
    timeDn = timeDn(better);
    timeUt = timeUt(better);
    mxGC   = mxGC(better);
    myGC   = myGC(better);
    mzGC   = mzGC(better);
    normGC = normGC(better);
    temp   = temp(better);

%% Get B-field magnitude vector
   if calcBfield == 0
%       Load the already-calculated B-field
        disp('Loading B-field magnitudes...');
        load('normB.mat');
   else
        disp('Calculating B-field magnitudes...');
%       Units = [Tesla] Input must be 1 row x N column
        [B_ECI, S_ECI] = GenBRVfromTle(timeDn,tle);
        B_ECI = B_ECI.*1e4*1e3; 
        normB_mG = sqrt(B_ECI(1,:).^2+B_ECI(2,:).^2+B_ECI(3,:).^2);%[mG]    

%       Save the B field is the option selected
        if saveBfield == 1
            save('normB.mat','normB_mG');
        end
   end

%% Plot Ground Calibration
    PlotMagDataOnUnitSphere(1,[mxGC myGC mzGC],normGC','GroundCal',fileDir);

    PlotNormError(2,normGC'-normB_mG,'GroundCal',fileDir);
    fprintf('\n');
  
%% Translate back to raw data
    [xRaw,yRaw,zRaw]=Uncalibrate(mxGC,myGC,mzGC,temp,mOffsetXYZ,mScaleXYZ,tScaleXYZ,rotMatVec);
    magTemp = temp; %[deg C]
    normRaw = sqrt(xRaw.^2+yRaw.^2+zRaw.^2); %[mG]
    
%% Show Rawified Data
    %group data as the calibration wants
    S(:,1) = xRaw;     %[mG]
    S(:,2) = yRaw;     %[mG]
    S(:,3) = zRaw;     %[mG]
    S(:,4) = normRaw;  %[mG]
    S(:,5) = magTemp;  %[C]

    R(:,1) = normB_mG; %[mG]  
    
    PlotMagDataOnUnitSphere(2,S,normB_mG,'RawMagOut',fileDir);

    PlotNormError(21,normRaw-normB_mG,'RawMagOut',fileDir);
    
    fprintf('\n');
 
%% Perform Calibration
    options = optimset('Display', 'off');

    % do scale / offset calibration
    C0 = [mean(xRaw) mean(yRaw) mean(zRaw) 1 1 1];
    [C] = lsqnonlin(@(C)fmagcalib_orbit(C,S,R),C0,[],[],options);
    % correct the data
    [~,resc,xc,yc,zc]=fmagcalib_orbit(C,S,R);
    normCal = sqrt(xc.^2+yc.^2+zc.^2)';
    PlotNormError(3,normCal-normB_mG,'Scale + Offset',fileDir);
    fprintf('offset = [%1.6f %1.6f %1.6f] %%[mG]\nscale = [%1.6f %1.6f %1.6f] %%[non-dimensional]\n\n',...
            C(1),C(2),C(3),C(4),C(5),C(6));
    
    % do scale / offset / temperature calibration
    Ct0 = [mean(xRaw) mean(yRaw) mean(zRaw) 1 1 1 0 0 0];
    [Ct] = lsqnonlin(@(Ct)fmagcalibTemp_orbit(Ct,S,R),Ct0,[],[],options);
    % correct the data
    [~,resct,xct,yct,zct]=fmagcalibTemp_orbit(Ct,S,R);
    normCalt = sqrt(xct.^2+yct.^2+zct.^2)';
    PlotNormError(4,normCalt-normB_mG,'Scale + Offset + Temp',fileDir);
    fprintf('offset = [%1.6f %1.6f %1.6f] %%[mG]\nscale = [%1.6f %1.6f %1.6f] %%[non-dimensional]\n',...
            Ct(1),Ct(2),Ct(3),Ct(4),Ct(5),Ct(6));
    fprintf('tempScale = [%1.6e %1.6e %1.6e] %%[mG/degC]\n\n',Ct(7),Ct(8),Ct(9));
    
    % do scale / offset / non-orthogonal calibration
    Cn0 = [mean(xRaw) mean(yRaw) mean(zRaw) 1 1 1 0 0 0];
    [Cn] = lsqnonlin(@(Cn)fmagcalibNonOrthog_orbit(Cn,S,R),Cn0,[],[],options);
    % correct the data
    [~,resct,xcn,ycn,zcn]=fmagcalibNonOrthog_orbit(Cn,S,R);
    normCalt = sqrt(xcn.^2+ycn.^2+zcn.^2)';
    PlotNormError(5,normCalt-normB_mG,'Scale + Offset + Orthog',fileDir);
    fprintf('offset = [%1.6f %1.6f %1.6f] %%[mG]\nscale = [%1.6f %1.6f %1.6f] %%[non-dimensional]\n',...
            Cn(1),Cn(2),Cn(3),Cn(4),Cn(5),Cn(6));
    fprintf('rho=%1.2f*pi/180; phi=%1.2f*pi/180; lambda=%1.2f*pi/180; %%[radians]\n\n',Cn(7)*180/pi,Cn(8)*180/pi,Cn(9)*180/pi);

    % do scale / offset / temperature / non-orthogonal calibration
    Ctn0 = [mean(xRaw) mean(yRaw) mean(zRaw) 1 1 1 0 0 0 0 0 0];
    [Ctn] = lsqnonlin(@(Ctn)fmagcalibTempNonOrthog_orbit(Ctn,S,R),Ctn0,[],[],options);
    % correct the data
    [~,resctn,xctn,yctn,zctn]=fmagcalibTempNonOrthog_orbit(Ctn,S,R);
    normCaltn = sqrt(xctn.^2+yctn.^2+zctn.^2)';
    PlotNormError(6,normCaltn-normB_mG,'Scale + Offset + Temp + Orthog',fileDir);
    fprintf('offset = [%1.6f %1.6f %1.6f] %%[mG]\nscale = [%1.6f %1.6f %1.6f] %%[non-dimensional]\n',...
            Ctn(1),Ctn(2),Ctn(3),Ctn(4),Ctn(5),Ctn(6));
    fprintf('tempScale = [%1.6e %1.6e %1.6e] %%[mG/degC]\n',Ctn(7),Ctn(8),Ctn(9));
    fprintf('rho=%1.2f*pi/180; phi=%1.2f*pi/180; lambda=%1.2f*pi/180; %%[radians]\n\n',Ctn(10)*180/pi,Ctn(11)*180/pi,Ctn(12)*180/pi);

    % do scale / offset / rotational calibration
    Cf0 = [mean(xRaw) mean(yRaw) mean(zRaw) 1 1 1  0 0 0 ];
    LB = [-Inf -Inf -Inf 0.5  0.5  0.5 -pi -pi -pi];
    UB = [ Inf  Inf  Inf 2.0  2.0  2.0  pi  pi  pi];
    [Cf] = lsqnonlin(@(Cf)fmagcalibRot_orbit(Cf,S,R),Cf0,LB,UB,options);
    % correct the data
    [~,rescf,xcf,ycf,zcf,af]=fmagcalibRot_orbit(Cf,S,R);
    normCalf = sqrt(xcf.^2+ycf.^2+zcf.^2)';
    PlotNormError(7,normCalf-normB_mG,'Scale + Offset + Rotation',fileDir);
    fprintf('offset = [%1.6f %1.6f %1.6f] %%[mG]\nscale = [%1.6f %1.6f %1.6f] %%[non-dimensional]\n',...
            Cf(1),Cf(2),Cf(3),Cf(4),Cf(5),Cf(6));
    fprintf('phi=%1.2f*pi/180; theta=%1.2f*pi/180; psi=%1.2f*pi/180; %%[radians]\n',Cf(7)*180/pi,Cf(8)*180/pi,Cf(9)*180/pi);
    fprintf('Rotation Matrix:\n');
    fprintf('rotate = [%1.6f %1.6f %1.6f %1.6f %1.6f %1.6f %1.6f %1.6f %1.6f]\n',af(1,1),af(1,2),af(1,3),af(2,1),af(2,2),af(2,3),af(3,1),af(3,2),af(3,3));
    disp(af);

    % do scale / offset / temp / rotational calibration
    Ctf0 = [mean(xRaw) mean(yRaw) mean(zRaw) 1 1 1  0 0 0 0 0 0];
    LB = [-Inf -Inf -Inf 0.5  0.5  0.5 -Inf -Inf -Inf -pi -pi -pi];
    UB = [ Inf  Inf  Inf 2.0  2.0  2.0  Inf  Inf  Inf  pi  pi  pi];
    [Ctf] = lsqnonlin(@(Ctf)fmagcalibTempRot_orbit(Ctf,S,R),Ctf0,LB,UB,options);
    % correct the data
    [~,resctf,xctf,yctf,zctf,atf]=fmagcalibTempRot_orbit(Ctf,S,R);
    normCalf = sqrt(xctf.^2+yctf.^2+zctf.^2)';
    PlotNormError(8,normCalf-normB_mG,'Scale + Offset + Temp + Rotation',fileDir);
    fprintf('offset = [%1.6f %1.6f %1.6f] %%[mG]\nscale = [%1.6f %1.6f %1.6f] %%[non-dimensional]\n',...
            Ctf(1),Ctf(2),Ctf(3),Ctf(4),Ctf(5),Ctf(6));
    fprintf('tempScale = [%1.6f %1.6f %1.6f] %%[mG/degC]\n',Ctf(7),Ctf(8),Ctf(9));
    fprintf('phi=%1.2f*pi/180; theta=%1.2f*pi/180; psi=%1.2f*pi/180; %%[radians]\n',Ctf(10)*180/pi,Ctf(11)*180/pi,Ctf(12)*180/pi);
    fprintf('Rotation Matrix (mag frame to soft-magnetic frame):\n');
    fprintf('rotate = [%1.6f %1.6f %1.6f %1.6f %1.6f %1.6f %1.6f %1.6f %1.6f]\n',atf(1,1),atf(1,2),atf(1,3),atf(2,1),atf(2,2),atf(2,3),atf(3,1),atf(3,2),atf(3,3));
    disp(atf);
    
    PlotMagDataOnUnitSphere(9,[xctf yctf zctf],normB_mG','OnOrbitCal',fileDir);

    
end

function dateNumOut=unix2datenum(unixTimeSec)
%this function converts unix time to datenum.
    dateNumOut = unixTimeSec/86400+719529; %[matlab dateNum] 719529 = datenum(1970,1,1)
end

function unixTimeSec=datenum2unix(dateNumIn)
%Convert from datenum to unix time
    unixTimeSec=(dateNumIn-719529)*86400; %[matlab dateNum] 719529 = datenum(1970,1,1)
end

function [err,res,xCorr,yCorr,zCorr] = fmagcalib_orbit(X, S, r)
% Offset and scale correction
    x0 = X(1);
    y0 = X(2);
    z0 = X(3);
    sx = X(4);
    sy = X(5);
    sz = X(6);
    xRaw=S(:,1);
    yRaw=S(:,2);
    zRaw=S(:,3);
    
    xCorr = sx*(xRaw + x0);
    yCorr = sy*(yRaw + y0);
    zCorr = sz*(zRaw + z0);
    
    err = (r.^2 - (xCorr.^2 + yCorr.^2 + zCorr.^2));
    res=sum(err);
end

function [err,res,xCorr,yCorr,zCorr] = fmagcalibTemp_orbit(X, S, r)
% Offset, scale and temperature correction
    x0 = X(1);
    y0 = X(2);
    z0 = X(3);
    sx = X(4);
    sy = X(5);
    sz = X(6);
    sXt = X(7);
    sYt = X(8);
    sZt = X(9);
    xRaw=S(:,1);
    yRaw=S(:,2);
    zRaw=S(:,3);
    temp=S(:,5);
    
    xCorr = sx*(xRaw.*(1+sXt*temp) + x0);
    yCorr = sy*(yRaw.*(1+sYt*temp) + y0);
    zCorr = sz*(zRaw.*(1+sZt*temp) + z0);
    
    err = (r.^2 - (xCorr.^2 + yCorr.^2 + zCorr.^2));
    res=sum(err);
end

function [err,res,xCorr,yCorr,zCorr] = fmagcalibNonOrthog_orbit(X, S, r)
% Offset, scale and non-orthogonality correction
    % USES THE NEGATIVE OFFSET

    x0 = X(1);
    y0 = X(2);
    z0 = X(3);
    sx = X(4);
    sy = X(5);
    sz = X(6);
    rho = X(7);
    phi = X(8);
    lambda = X(9);

    xRaw=S(:,1); %[mG]
    yRaw=S(:,2); %[mG]
    zRaw=S(:,3); %[mG]
    
    xCorr = sx*(xRaw-x0);
    yCorr = 1/cos(rho)*(sy*(yRaw-y0)-xCorr*sin(rho));
    zCorr = 1/(cos(phi)*cos(lambda))*(sz*(zRaw-z0)-xCorr*sin(lambda)-yCorr*sin(phi)*cos(lambda));
    
    err = (r.^2 - (xCorr.^2 + yCorr.^2 + zCorr.^2));
    res=sum(err);
end

function [err,res,xCorr,yCorr,zCorr] = fmagcalibTempNonOrthog_orbit(X, S, r)
% Offset, scale, temperature and non-orthogonality correction
% USES NEGATIVE OFFSET
    x0 = X(1);
    y0 = X(2);
    z0 = X(3);
    sx = X(4);
    sy = X(5);
    sz = X(6);
    sXt = X(7);
    sYt = X(8);
    sZt = X(9);
    rho = X(10);
    phi = X(11);
    lambda = X(12);

    xRaw=S(:,1); %[mG]
    yRaw=S(:,2); %[mG]
    zRaw=S(:,3); %[mG]
    temp=S(:,5); %[degC]
    
    xCorr = sx*(xRaw.*(1+sXt*temp)-x0);
    yCorr = 1/cos(rho)*(sy*(yRaw.*(1+sYt*temp)-y0)-xCorr*sin(rho));
    zCorr = 1/(cos(phi)*cos(lambda))*(sz*(zRaw.*(1+sZt*temp)-z0)-xCorr*sin(lambda)-yCorr*sin(phi)*cos(lambda));
    
    err = (r.^2 - (xCorr.^2 + yCorr.^2 + zCorr.^2));
    res=sum(err);
end

function [err,res,xCorr,yCorr,zCorr,a] = fmagcalibRot_orbit(X, S, r)
% Offset, scale and rotation correction
    x0 = X(1);
    y0 = X(2);
    z0 = X(3);
    sx = X(4);
    sy = X(5);
    sz = X(6);
    phi    = X(7);
    theta  = X(8);
    psi    = X(9);

    x=S(:,1);
    y=S(:,2);
    z=S(:,3);
  
    xOff = x + x0;
    yOff = y + y0;
    zOff = z + z0;
  
    cphi = cos(phi); 
    sphi = sin(phi); 
    cth  = cos(theta); 
    sth  = sin(theta); 
    cpsi = cos(psi); 
    spsi = sin(psi); 
  
    a = [... 
	cpsi*cth  -spsi*cphi+cpsi*sth*sphi  spsi*sphi+cpsi*cphi*sth 
	spsi*cth  cpsi*cphi+sphi*sth*spsi   -cpsi*sphi+sth*spsi*cphi 
	-sth      cth*sphi                  cth*cphi ];
    
    A = a * [xOff yOff zOff]';
    A = a^-1 * [sx*A(1,:);sy*A(2,:);sz*A(3,:)];
    xCorr = A(1,:)';
    yCorr = A(2,:)';
    zCorr = A(3,:)';
    
    err = r.^2 - (xCorr.^2 + yCorr.^2 + zCorr.^2);
    res=sum(err);
end

function [err,res,xCorr,yCorr,zCorr,a] = fmagcalibTempRot_orbit(X, S, r)
% Offset, scale, temperature and rotation correction

    x0 = X(1);
    y0 = X(2);
    z0 = X(3);
    sx = X(4);
    sy = X(5);
    sz = X(6);
    sXt = X(7);
    sYt = X(8);
    sZt = X(9);
    phi    = X(10);
    theta  = X(11);
    psi    = X(12);

    xRaw=S(:,1);
    yRaw=S(:,2);
    zRaw=S(:,3);
    temp=S(:,5);
  
    xOff = (xRaw.*(1+sXt*temp) + x0);
    yOff = (yRaw.*(1+sYt*temp) + y0);
    zOff = (zRaw.*(1+sZt*temp) + z0);
  
    cphi = cos(phi); 
    sphi = sin(phi); 
    cth  = cos(theta); 
    sth  = sin(theta); 
    cpsi = cos(psi); 
    spsi = sin(psi); 
  
    a = [... 
	cpsi*cth  -spsi*cphi+cpsi*sth*sphi  spsi*sphi+cpsi*cphi*sth 
	spsi*cth  cpsi*cphi+sphi*sth*spsi   -cpsi*sphi+sth*spsi*cphi 
	-sth      cth*sphi                  cth*cphi ];

    A = a * [xOff yOff zOff]';
    A = a^-1 * [sx*A(1,:); sy*A(2,:); sz*A(3,:)];
    xCorr = A(1,:)';
    yCorr = A(2,:)';
    zCorr = A(3,:)';
    
    err = r.^2 - (xCorr.^2 + yCorr.^2 + zCorr.^2);
    res=sum(err);
end

function [err,res,xCorr,yCorr,zCorr,a] = fmagcalibTempRotEcl_orbit(X, S, r)
% Offset, scale, temperature, rotation and eclipse correction
    x0 = X(1);
    y0 = X(2);
    z0 = X(3);
    sx = X(4);
    sy = X(5);
    sz = X(6);
    sXt = X(7);
    sYt = X(8);
    sZt = X(9);
    phi    = X(10);
    theta  = X(11);
    psi    = X(12);
    oEx = X(13);
    oEy = X(14);
    oEz = X(15);

    xRaw=S(:,1);
    yRaw=S(:,2);
    zRaw=S(:,3);
    temp=S(:,5);
    ecl =S(:,6);
  
    xOff = (xRaw.*(1+sXt*temp) + x0 + oEx*ecl);
    yOff = (yRaw.*(1+sYt*temp) + y0 + oEy*ecl);
    zOff = (zRaw.*(1+sZt*temp) + z0 + oEz*ecl);
  
    cphi = cos(phi); 
    sphi = sin(phi); 
    cth  = cos(theta); 
    sth  = sin(theta); 
    cpsi = cos(psi); 
    spsi = sin(psi); 
  
    a = [... 
	cpsi*cth  -spsi*cphi+cpsi*sth*sphi  spsi*sphi+cpsi*cphi*sth 
	spsi*cth  cpsi*cphi+sphi*sth*spsi   -cpsi*sphi+sth*spsi*cphi 
	-sth      cth*sphi                  cth*cphi ];

    A = a * [xOff yOff zOff]';
    A = a^-1 * [sx*A(1,:); sy*A(2,:); sz*A(3,:)];
    xCorr = A(1,:)';
    yCorr = A(2,:)';
    zCorr = A(3,:)';
    
    err = r.^2 - (xCorr.^2 + yCorr.^2 + zCorr.^2);
    res=sum(err);
end

function [xRaw,yRaw,zRaw]=Uncalibrate(mx,my,mz,temp,mOffsetXYZ,mScaleXYZ,tScaleXYZ,rotMatVec)
 %this function converts from calibrated to raw sensor data
    %first, form the rotation matrix
    a = [rotMatVec(1) rotMatVec(2) rotMatVec(3);
         rotMatVec(4) rotMatVec(5) rotMatVec(6);
         rotMatVec(7) rotMatVec(8) rotMatVec(9)];

    %rotate into soft-magnet frame
    magInSoftFrame=a*[mx'; my'; mz'];
    
    %remove scale
    magInSoftFrameNoScale=[magInSoftFrame(1,:)/mScaleXYZ(1); magInSoftFrame(2,:)/mScaleXYZ(2); magInSoftFrame(3,:)/mScaleXYZ(3)];

    %rotate back to magnetometer frame
    magNoScale=(a^-1)*magInSoftFrameNoScale;
    
    %remove offset and temperature scale
    xRaw = (magNoScale(1,:)-mOffsetXYZ(1))./(1+tScaleXYZ(1)*temp');
    yRaw = (magNoScale(2,:)-mOffsetXYZ(2))./(1+tScaleXYZ(2)*temp');
    zRaw = (magNoScale(3,:)-mOffsetXYZ(3))./(1+tScaleXYZ(3)*temp');

end

function PlotMagDataOnSphere(figNum,S,radius)
 %Plot the information on a sphere or radius ~=1
    figure(figNum); clf;
    % draw data
    plot3( S(:,1), S(:,2), S(:,3), '.r' ); hold on;
    [X,Y,Z]=sphere(20);
    M=mesh(X*radius,Y*radius,Z*radius);
    set(M,'facecolor','none','edgecolor',0.5*[1 1 1]);
    grid;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    hold off;
    axis equal;
end

function PlotMagDataOnUnitSphere(figNum,S,normVec,name,fileDir)
 %Plot the information on a sphere or radius =1
    fN=figure(figNum); clf;

    % draw data
    plot3( S(:,1)./normVec', S(:,2)./normVec', S(:,3)./normVec', '.r' ); hold on;
    [X,Y,Z]=sphere(20);
    M=mesh(X,Y,Z);
    set(M,'facecolor','none','edgecolor',0.5*[1 1 1]);
    grid;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    hold off;
    axis equal;
    title(['Unit Sphere - ',name]);
    print(fN,'-dpng','-r250',[fileDir,'UnitSphere_',name,'.png']);
end

function PlotNormError(figNum,normErr,titleStr,fileDir)
    histBins = 100;
    lineWidth = 1;
    
    [~,~]=mkdir(fileDir);
    fN=figure(figNum); clf;
    subplot(2,1,1); hold on;
    plot(1:length(normErr),normErr,'k');
    plot([1 length(normErr)],nanmean(normErr,2)*[1 1],'g','LineWidth',lineWidth);
    plot([1 length(normErr)],nanstd(normErr,2)*[1 1],'c','LineWidth',lineWidth);
    xlabel('Sample'); ylabel('Error [mG]');
    title({titleStr,['std = ',num2str(nanstd(normErr),'%1.2f'),' mG | mean = ',num2str(nanmean(normErr,2),'%1.2f'),' mG']});
    grid on;
    subplot(2,1,2); hold on;
    hist(normErr,histBins);
    [N,X]=hist(normErr,histBins);
    plot(nanmean(normErr,2)*[1 1],[1 max(N)],'g','LineWidth',lineWidth);
    [~,modeInd] = max(N);
    plot(X(modeInd)*[1 1],[1 max(N)],'r--','LineWidth',lineWidth);
    plot(nanstd(normErr)*[1 1],[1 max(N)],'c','LineWidth',lineWidth);
    xlabel('Error [mG]'); ylabel('Number of Samples');
    legend('dat','mean','mode','std','Location','Best');
    subplot(2,1,1);
    plot([1 length(normErr)],X(modeInd)*[1 1],'r--','LineWidth',lineWidth);
    plot([1 length(normErr)],std(normErr)*[1 1],'c','LineWidth',lineWidth);
    
    percBelow1Sig = length(find(abs(normErr)<=nanstd(normErr)))/length(normErr)*100;
    fprintf([titleStr,': < 1sig: %2.2f | std: %2.2f mG | mean: %2.2f mG\n'],percBelow1Sig,nanstd(normErr),nanmean(normErr));
    print(fN,'-dpng','-r300',[fileDir,'NormError_',titleStr,'.png']);
end

function PlotAxesOfMag(figNum,timeDn,x,y,z,titleStr)
 
 normMag = sqrt(x.^2+y.^2+z.^2);
 figure(figNum)
 plot(timeDn,x,timeDn,y,timeDn,z,timeDn,normMag);
 title(titleStr);

end