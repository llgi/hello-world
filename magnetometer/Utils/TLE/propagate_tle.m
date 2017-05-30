function [x_ecf_km x_eci_km jd v_eci_kmps] = propagate_tle(tle,t_min)

%% Inputs
% tle: structure containing strings for each line of the tle
%      tle.line1, tle.line2
% t_min: time vector in minutes from epoch

%% Outputs
% x_ecf_km: satellite position in ECF (km)
% x_eci_km: satellite position in ECI (km)

satrec = twoline2rvMOD(tle.line1,tle.line2);
jd0 = satrec.jdsatepoch;

x_ecf_km = zeros(length(t_min),3);
x_eci_km = zeros(length(t_min),3);

for i = 1:length(t_min)
%     [satrec, x_ecf_km(i,:)] = spg4_ecf(satrec,t_min(i));
    jd(i) = jd0 + t_min(i)/60/24;
%     x_eci_km(i,:) = ecf2eci(x_ecf_km(i,:),jd(i));
    [satrec, x_eci_km(i,:), v_eci_kmps(i,:)] = sgp4(satrec,t_min(i));
    x_ecf_km(i,:) = eci2ecf(x_eci_km(i,:),jd(i));
end
