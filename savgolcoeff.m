function c = savgolcoeff(wl, wr, od, m)
%% SAVGOLCOEFF Calculates the coefficients for the Savitky-Golay filter
%
%  Calculates the coefficients for the Savitzky-Golay filter [1].
%
%  Arguments: wl - window: points used at the left
%             wr - window: points used at the right
%             od - order of derivative
%             m  - degree of smoothing polynomial
%
%  Returns:   c  - Savitzky-Golay coefficients
%
%  Reference:
%  [1]  A. Savitzky, M.J.E. Golay, Anal. Chem. 36, 1627 (1964)
%       doi: 10.1021/ac60214a047
%
%  Copyright 2020, Philipp Rahe
%
%  This is a script under the terms of the Creative Commons Attribution 
%  License (creativecommons.org/licenses/by/4.0), which permits 
%  unrestricted use, distribution, and reproduction in any medium, 
%  provided the original work is properly cited.
%   
%  14.08.2020, Philipp Rahe (prahe@uos.de)

  % check input arguments for consistency
  if ( (wl<0) || (wr<0) || (od>m) || (wl+wr < m) )
    error('invalid input arguments.');
  end
    
  % number of output elements in c
  nc = wl+wr+1;

  % initialise arrays
  A = zeros(m+1, m+1);
  b = zeros(m+1, 1);
  c = zeros(nc, 1);

  % equations for least-squares fit
  for ip=0:(2*m)
    lsfsum = 1.0;
    if(ip>0), lsfsum=0.0; end
    for k=1:wr, lsfsum = lsfsum + k^ip; end
    for k=1:wl, lsfsum = lsfsum + (-k)^ip; end
    mm = min([ip, 2*m-ip]);
    for im=-mm:2:mm
      A(1+(ip+im)/2, 1+(ip-im)/2) = lsfsum;
    end
  end
  
  % solve the normal equations
  b(od+1) = 1.0;
  b = A\b;
  
  % calculate the coefficients
  for k=-wl:wr
    lsfsum = b(1);
    f = 1.0;
    for mm=1:m
      f = f*k;
      lsfsum = lsfsum + b(1+mm)*f;
    end
    % change order
    c(k+wl+1) = lsfsum;
  end
  
  % if a derivatives of order k is requested, multiply c by (k!) .
  if(od > 1)
    c = c.* factorial(od);
  end
  
end

