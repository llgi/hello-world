function K = schmidt(nmax,mmax)

%+=====================================================================+
%
%     Programmers:  Carlos Roithmayr          			Feb 1997
%
%		    NASA Langley Research Center
%		    Spacecraft and Sensors Branch (CBC)
%		    757 864 6778
%		    c.m.roithmayr@larc.nasa.gov
%
%+---------------------------------------------------------------------+
%
%     Purpose:
%
%     Compute coefficients that relate Schmidt functions to associated
%     Legendre functions.
%
%+---------------------------------------------------------------------+
%
%     Argument definitions:
%
%     nmax              Maximum degree of contributing spherical harmonics
%
%     mmax              Maximum order of contributing spherical harmonics
%
%     K		            coefficients that relate Schmidt functions to
%	                    associated Legendre functions (Ref. [1]).
%
%+---------------------------------------------------------------------+
%
%     References:
%
%     1. Haymes, R. C., Introduction to Space Science, Wiley, New
%        York, 1971.
%
%     2. Roithmayr, C., "Contributions of Spherical Harmonics to
%        Magnetic and Gravitational Fields", EG2-96-02, NASA Johnson
%        Space Center, Jan. 23, 1996.
%
%+=====================================================================+

% The number 1 is added to degree and order since MATLAB can't have an array
% index of 0.


% Seed for recursion formulae
K(2,2) = 1;

% Recursion formulae

for n = 1:nmax
    i=n+1;

  for m = 0:n
     j=m+1;

    if m == 0
   	% Eq. (3), Ref. [2]
   	  K(i,j) = 1;

   	elseif ((m >= 1) & (n >= (m+1)))
    % Eq. (4), Ref. [2]
   	  K(i,j) = sqrt((n-m)/(n+m))*K(i-1,j);

   	elseif ((m >= 2) & (n >= m))
    % Eq. (5), Ref. [2]
   	  K(i,j) = K(i,j-1)/sqrt((n+m)*(n-m+1));
   	end

  end
end


