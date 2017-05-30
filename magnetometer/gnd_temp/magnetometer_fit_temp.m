function [xcc,ycc,zcc,ncc] = magnetometer_fit_temp(file)
% Extract magnetometer calibration data
% magnetometer_fit(file)
% file: comma-seperated-file with three columns generated by 'mag loopraw'
% plus and extra column for temperature data
% returns [x,y,z,norm] in calibrated form

% Read the file and save all the inputs in "data" variable
data=dlmread(file);

% Plot the input data
figure(5)
plot(data);
grid
xlabel('Samples')
s=input('Input start sample: ');
e=input('Input end sample: ');

x    = data(s:e,1);
y    = data(s:e,2);
z    = data(s:e,3);
n    = data(s:e,4);
temp = data(s:e,5);

r = mean(n)
[X,Y,Z]=sphere(30);

%Array with all the input data
S=[x y z temp];

% Scale, offset and temperature calibration
C0 = [mean(x) mean(y) mean(z) 1 1 1 0 0 0];
[C] = lsqnonlin(@(C)fmagcalibtemp(C,S,r),C0);

% Correct the data
[~,resc,xc,yc,zc]=fmagcalibtemp(C,S,r);
nc = (xc.*xc + yc.*yc + zc.*zc).^0.5;

% Scale, offset, temperature and rotation correction (full calibration)
Cf0 = [mean(x) mean(y) mean(z) 1 1 1 0 0 0 0 0 0 ];
[Cf] = lsqnonlin(@(Cf)fmagcalibtemp_rot(Cf,S,r),Cf0);

% Correct the data
[~,rescc,xcc,ycc,zcc]=fmagcalibtemp_rot(Cf,S,r);
ncc = (xcc.*xcc + ycc.*ycc + zcc.*zcc).^0.5;

figure(1)
% draw data
plot3( x, y, z, '.r' );
hold on;
% full correction data
plot3( xcc, ycc, zcc, '.r' );
hold on;
% scale, offset and temperature correction data
plot3( xc, yc, zc, '.b' );
hold on;
M=mesh(X*r,Y*r,Z*r);
set(M,'facecolor','none')
grid
xlabel('x')
ylabel('y')
zlabel('z')
title('Raw magnetometer readings [mG]', 'Mag data (full calib)','Mag data (only offset)')
hold off
axis equal

% Calibration parameters
offset = Cf(1:3);
scale = Cf(4:6);
phi    = Cf(7);
theta  = Cf(8);
psi    = Cf(9);
cphi = cos(phi);
sphi = sin(phi);
cth  = cos(theta);
sth  = sin(theta);
cpsi = cos(psi);
spsi = sin(psi);

A = [cpsi*cth  -spsi*cphi+cpsi*sth*sphi  spsi*sphi+cpsi*cphi*sth;
     spsi*cth  cpsi*cphi+sphi*sth*spsi   -cpsi*sphi+sth*spsi*cphi;
     -sth      cth*sphi                  cth*cphi ];

% Print the calibration parameters
format short g
disp(sprintf('magOffset ='))
disp(-offset)
disp(sprintf('magScale='))
disp( scale )
disp(sprintf('magRotate='))
disp( A )
disp(sprintf('magRotateTranspose:'))
disp( A' )

% Write them in a file
fileID = fopen('calibVal.txt','w');
fprintf(fileID,'param set mag_offset %6.3f %6.3f %6.3f\n', -offset(1), -offset(2), -offset(3));
fprintf(fileID,'param set mag_scale %6.3f %6.3f %6.3f\n', scale(1), scale(2), scale(3));
fprintf(fileID,'param set mag_rotate %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n',...
    A(1,1), A(1,2), A(1,3), A(2,1), A(2,2), A(2,3), A(3,1), A(3,2), A(3,3));
fclose(fileID);


disp(sprintf('Pre calib std %f mG',std(nc)))
disp(sprintf('Post calib std %f mG',std(ncc)))


figure(3);
subplot(3,1,1)
plot(x,ncc,'r',x,n,'.')
title('x')
grid
subplot(3,1,2)
plot(y,ncc,'r',y,n,'.')
grid
title('y')
subplot(3,1,3)
plot(z,ncc,'r',z,n,'.')
grid
title('z')
legend('Calibrated Norm','Non-calib Norm')

figure(4)
subplot(3,2,1)
plot(x,y,'*')
axis equal
grid
xlabel('x')
ylabel('y')
title('Non calibrated [mG]')
subplot(3,2,2)
plot(xcc,ycc,'*')
axis equal
grid
xlabel('x')
ylabel('y')
title('Calibrated  [mG]')
subplot(3,2,3)
plot(x,z,'*')
axis equal
grid
xlabel('x')
ylabel('z')
subplot(3,2,4)
plot(xcc,zcc,'*')
axis equal
grid
xlabel('x')
ylabel('z')
subplot(3,2,5)
plot(y,z,'*')
axis equal
grid
xlabel('y')
ylabel('z')
subplot(3,2,6)
plot(ycc,zcc,'*')
axis equal
grid
xlabel('y')
ylabel('z')

figure(6)
plot(ncc)
hold on
plot(n,'r');
hold off
grid
legend('Calibrated Norm','Non-calib Norm')
end
