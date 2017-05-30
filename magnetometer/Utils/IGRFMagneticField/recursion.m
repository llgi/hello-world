function [A,ctilde,stilde] = recursion(repe,nmax,mmax)
%+---------------------------------------------------------------------+
%
%     Purpose:
%
%     Recursive calculations of derived Legendre polynomials and other
%     quantities needed for gravitational and magnetic fields.
%
%+---------------------------------------------------------------------+

% The number 1 is added to degree and order since MATLAB can't have an
% array index of 0.

    clear A;
    A=zeros(nmax+3,nmax+3);         % A(n,m) = 0, for m > n

    R_m = sqrt(repe'*repe);
    rhat = repe/R_m;

    u = rhat(3);                    % sin of latitude

    A(1,1)=1;                       % "derived" Legendre polynomials
    A(2,1)=u;
    A(2,2)=1;
    clear ctilde
    clear stilde
    ctilde(1) = 1; ctilde(2) = rhat(1);
    stilde(1) = 0; stilde(2) = rhat(2);

    for n = 2:nmax
        i=n+1;

        % Calculate derived Legendre polynomials and "tilde" letters
        % required for gravitational and magnetic fields.

        % Eq. (4a), Ref. [2]
        A(i,i) = prod(1:2:(2*n - 1));

        % Eq. (4b), Ref. [2]
        A(i,(i-1))= u*A(i,i);

        if n <= mmax
            %   p. 9,     Ref. [1]
            ctilde(i)  = ctilde(2) * ctilde(i-1) - stilde(2) * stilde(i-1);
            stilde(i)  = stilde(2) * ctilde(i-1) + ctilde(2) * stilde(i-1);
        end

        for m = 0:n
            j=m+1;
            if (m < (n-1)) && (m <= (mmax+1))
                %     Eq. I, Table 1, Ref. [2]
                A(i,j)=((2*n - 1)*u*A((i-1),j) - (n+m-1)*A((i-2),j))/(n-m);
            end
        end
    end
end