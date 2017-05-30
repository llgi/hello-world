function [B_ECI, S_ECI, posECI_km, posECEF_km, velECI_kmps] = ...
  GenBRVfromTle(timeDateNum, tle)
%This function generates a vector of the B-field vs. time which the
%DualAttAndParamEst function uses to lookup H & Hdot given t

% INPUTS:
%   timeDateNum : an array containing all the time stamps where the magnetic 
%                 should be determined
%   tle : two line elements of the satellite. It should be closer as possible 
%         to the time interval to get better results.
  
%% Determine ECI position vector vs. time
    [posECEF_km posECI_km S_ECI GST velECI_kmps] = OrbitProp(timeDateNum, tle);

%% Use IGRF to Determine ECEF H-field vector vs. time
    B_ECEF = IGRF( posECEF_km ); %[Tesla]

%% Convert ECEF to ECI H-field vector
    B_ECI = zeros(3,size(B_ECEF,2));
    for i=1:size(posECEF_km,2)
        R_E2I = [cos(GST(i)) -sin(GST(i)) 0; sin(GST(i)) cos(GST(i)) 0; 0 0 1];
        B_ECI(:,i)=R_E2I*B_ECEF(:,i); %[Tesla]
    end
    
end

function [posECEF_km posECI_km posEarth2SunECI_unit GST velECI_kmps] = ...
  OrbitProp( time, tle )
%% INPUTS:
%   time
%   TLE
%   orbitElem - structure with
%     (1) per [km]
%     (2) ap [km]
%     (3) incl [deg]
%     (4) argPer [deg]
%     (5) RAAN [deg]
%     (6) epoch [dateNum]
%     (7) n [rad/s]
%     (8) ecc [deg]
%
%% OUTPUTS:
%   posECEF_km: [km] satellite position vector in ECEF frame, %LLGI
%   posECI_km: [km] satellite position vector in ECI frame, cartesian coordinates
%   posEarth2SunECI_unit: unit vector from earth to sun
%   GST: Greenwich Sideral Time 
%   velECI_kmps: [km/s] satellite velocity in ECI frame

    %% Constants
    % distance from earth to sun
      AU = 149597871; %[km]
    % average sun radius
      Rs_avg = 696342; %[km]
    % average earth radius 
      Re_avg = 6371; %[km]

    posECEF_km  = zeros(3,length(time));
    posECI_km   = zeros(3,length(time));
    velECI_kmps = zeros(3,length(time));
    
    %TLE epoch
    dateNumEpYear = datenum(['1 Jan 20',tle.line1(19:20)]);
    dateNumEpoch = dateNumEpYear + str2double(tle.line1(21:32));   
      
    for i=1:length(time)

      
      %epoch = get_tle_epoch(tle); %epoch in julian date
      t_min = (time(i)-dateNumEpoch)*24*60; % propagator uses time in minutes since tle epoch
      
      % Propagate tle
      [posECEF_km(:,i), posECI_km(:,i), ~, velECI_kmps(:,i)] = propagate_tle(tle,t_min);

    end

    %Convert position vector to ECEF frame
    %Determine Greenwich Sidereal Time Vector
    %(based on http://www.astro.umd.edu/~jph/GST_eqn.pdf)
    
    %GST at start of 2011
    GST_2011=6.6208844; %[hours]
    
    %day fractions from start of 2012 to each time
    dR2011 = time-datenum([2011 0 0 0 0 0]);
    
    %whole days from start of 2012 to each time
    dayR2011 = floor(dR2011);
    
    %hours from start of 2012 to each time
    hrR2011= mod(mod(dR2011,1)*24,24);
    
    GST = mod(GST_2011+0.0657098244*dayR2011+1.00273791*hrR2011,24)*pi/12; %[rad]
    
%% Compute sun direction in ECI frame

   %convert time vector to julian date
   tJ = juliandate(time);
   
   %calculate distance vector to sun 
   %  from http://www.dept.aoe.vt.edu/~cdhall/courses/aoe4140/attde.pdf
   T_UT1   = (tJ - 2451545)/36525;
   lambda_M_sun = 280.4606184+36000.77005361*T_UT1;
   M_sun   = 357.5277233+35999.05034*T_UT1;
   lambda  = lambda_M_sun + 1.914666471*sind(M_sun)+0.918994643*sind(2*M_sun);
   epsilon = 23.439291-0.0130042*T_UT1;
   
   %calculate unit vector from earth (or satellite) to sun
   posEarth2SunECI_unit = [cosd(lambda) ...
                           sind(lambda).*cosd(epsilon) ...
                           sind(lambda).*sind(epsilon)];
                           
   %calculate magnitude of earth to sun vector
   magS = AU*(1.000140612-0.016708617*cosd(M_sun)-0.000139589*cosd(2*M_sun));

%% Determine Umbra/Penumbra Eclipse Times
%  from http://celestrak.com/columns/v03n01/

    %position vector from satellite to sun
    posEarth2SunECI_km = [magS.*cosd(lambda) ...
                          magS.*sind(lambda).*cosd(epsilon) ...
                          magS.*sind(lambda).*sind(epsilon)];
    posSat2SunECI_km   = transpose(posEarth2SunECI_km) + posECI_km; %[km]

    %position vector from satellite to earth
    posSat2EarthECI_km = -posECI_km;
    
    %calculate distance magnitude of vectors for each timestep
    magE = sqrt( posSat2EarthECI_km(1,:).^2 + posSat2EarthECI_km(2,:).^2 + ...
                 posSat2EarthECI_km(3,:).^2 ); %[km]
                 
    %calculate semidiameters of earth & sun
    thetaE = asind(Re_avg./magE);
    thetaS = asind(Rs_avg./magS);
    
    %calculate the angle between the center of earth & the sun
    theta = acosd((posSat2EarthECI_km(1,:).*posSat2SunECI_km(1,:)+...
                   posSat2EarthECI_km(2,:).*posSat2SunECI_km(2,:)+...
                   posSat2EarthECI_km(3,:).*posSat2SunECI_km(3,:))./(magE.*magS));
    
    %calculate the inds which are in an umbral eclipse
    umbraInds = find(theta<thetaE-thetaS);
    
    %calculate the inds which are in penumbral eclipse
    penumbraInds = intersect(find(abs(thetaE-thetaS)<theta),find(theta<thetaE+thetaS));

%% Set sun vector during eclipse to zero
    posEarth2SunECI_unit(umbraInds, :)   = 0;
    posEarth2SunECI_unit(penumbraInds,:) = 0;
    
end

function [x_ecf_km x_eci_km jd v_eci_kmps] = propagate_tle( tle,t_min )
    %% INPUTS
    % tle: structure containing strings for each line of the tle
    %      tle.line1, tle.line2
    % t_min: time vector in minutes from epoch

    %% OUTPUTS
    % x_ecf_km: satellite position in ECF (km)
    % x_eci_km: satellite position in ECI (km)
    % jd
    % v_eci_kmps

    satrec = twoline2rvMOD(tle.line1,tle.line2);
    jd0    = satrec.jdsatepoch;

    x_ecf_km   = zeros(length(t_min),3);
    x_eci_km   = zeros(length(t_min),3);
    v_eci_kmps = zeros(length(t_min),3);

    for i = 1:length(t_min)

        jd(i) = jd0 + t_min(i)/60/24;

        [satrec, x_eci_km(i,:), v_eci_kmps(i,:)] = sgp4(satrec,t_min(i));
        x_ecf_km(i,:) = eci2ecf(x_eci_km(i,:),jd(i));
    end
    
    x_ecf_km   = x_ecf_km';
    x_eci_km   = x_eci_km';
    v_eci_kmps = v_eci_kmps';
end

function ShowSatSphereTrack_ECEF(posECEF_km)

    r=6371; %[km]
    
    figure(1);
    
    % draw data
    plot3(posECEF_km(1,:),posECEF_km(2,:),posECEF_km(3,:),'.r');
    hold on;
    [X,Y,Z]=sphere(30);
    surface(X*r,Y*r,Z*r,props);
    
    %shading flat;
    grid;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Sat Position, ECEF frame');
    hold off;
    axis equal;

end