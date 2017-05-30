function V_out = Math_Qinv_Vec_Q_Mult(V_in, Q_in)
%#eml

V_out = [0;0;0];

R0 = -Q_in(2) * V_in(3) + Q_in(3) * V_in(2) + Q_in(4) * V_in(1);
R1 =  Q_in(1) * V_in(3) - Q_in(3) * V_in(1) + Q_in(4) * V_in(2);
R2 = -Q_in(1) * V_in(2) + Q_in(2) * V_in(1) + Q_in(4) * V_in(3);
R3  = Q_in(1) * V_in(1) + Q_in(2) * V_in(2) + Q_in(3) * V_in(3);
       
V_out(1) =   R0 * Q_in(4) + R1 * Q_in(3) - R2 * Q_in(2) + R3 * Q_in(1);
V_out(2) = - R0 * Q_in(3) + R1 * Q_in(4) + R2 * Q_in(1) + R3 * Q_in(2);
V_out(3) =   R0 * Q_in(2) - R1 * Q_in(1) + R2 * Q_in(4) + R3 * Q_in(3);


% Alternate representation
% DCM(1,1) = q1^2 - q2^2 - q3^2 - q4^2;
% DCM(1,2) = 2*q1*q2 + 2*q3*q4;
% DCM(1,3) = 2*q1*q3 - 2*q2*q4;
% 
% DCM(2,1) = 2*q1*q2 - 2*q3*q4;
% DCM(2,2) = -q1^2 + q2^2 - q3^2 + q4^2;
% DCM(2,3) = 2*q1*q4 + 2*q2*q3;
% 
% DCM(1,1) = 2*q1*q3 + 2*q2*q4;
% DCM(1,2) = -2*q1*q4 + 2*q2*q3;
% DCM(1,3) = -q1^2 - q2^2 + q3^2 + q4^2;

% DCM(1,1) = Q_in(1)^2 - Q_in(2)^2 - Q_in(3)^2 + Q_in(4)^2;
% DCM(1,2) = 2*Q_in(1)*Q_in(2) + 2*Q_in(3)*Q_in(4);
% DCM(1,3) = 2*Q_in(1)*Q_in(3) - 2*Q_in(2)*Q_in(4);
% 
% DCM(2,1) = 2*Q_in(1)*Q_in(2) - 2*Q_in(3)*Q_in(4);
% DCM(2,2) = -Q_in(1)^2 + Q_in(2)^2 - Q_in(3)^2 + Q_in(4)^2;
% DCM(2,3) = 2*Q_in(1)*Q_in(4) + 2*Q_in(2)*Q_in(3);
% 
% DCM(3,1) = 2*Q_in(1)*Q_in(3) + 2*Q_in(2)*Q_in(4);
% DCM(3,2) = -2*Q_in(1)*Q_in(4) + 2*Q_in(2)*Q_in(3);
% DCM(3,3) = -Q_in(1)^2 - Q_in(2)^2 + Q_in(3)^2 + Q_in(4)^2;
% 
% 
% 
% V_out(1) = DCM(1,:)*V_in;
% V_out(2) = DCM(2,:)*V_in;
% V_out(3) = DCM(3,:)*V_in;