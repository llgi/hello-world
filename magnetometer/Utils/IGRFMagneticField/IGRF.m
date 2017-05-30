function B_ECEF=IGRF(posECEF_km)
% INPUT
%   posECEF_km : [km] satellite position vector in ECEF frame, cartesian
%                coordinates
% OUTPUT
%   B_ECEF : Magnetic field in [T]

    %% Constants
    % Mean radius for IGRF (6371.2 km)
    R_mean = 6371.2;
    % permeability of free space
    mu0 = 4*pi*1e-7;

    %% IGRF stuff
    % Load IGRF coefficients
    [G,H] = LoadCoIGRF(2015);
    % max degree of geopotential
    nmax = 10;
    % max order of geopotential
    mmax = 10;
    % call function to compute schmidt coefficients
    Kschmidt = schmidt(nmax,mmax);
    % define output vector
    B_ECEF = zeros(3,size(posECEF_km,2));
    
    for i=1:size(posECEF_km,2)
        % call function to compute legendre polynomials
        [A,ctilde,stilde] = recursion(posECEF_km(:,i),nmax,mmax);
        % determine b-field B_E
        bepe = bfield(posECEF_km(:,i),nmax,mmax,Kschmidt,A,ctilde,stilde,G,H,R_mean); %[Tesla]
        B_ECEF(:,i) = bepe; %[Tesla] ECEF frame
    end

end