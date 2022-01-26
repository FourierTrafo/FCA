function [xf, yf] = savgolfilter(x, y, F, ld, m)
%% SAVGOLFILTER Implementation of the Savitzky-Golay filter
% 
%  Implements the Savitzky-Golay filter as defined in [1].
%
%  Arguments: x  - vector with x values
%             y  - vector with y values
%             F  - full window size (symmetric window), must be odd
%             ld - order of derivative
%             m  - degree of smoothing polynomial
%
%  Returns:   xf - x axis of filtered data
%             yf - Savitzky-Golay filtered data
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
%  version 14.08.2020, Philipp Rahe (prahe@uos.de)
  
  % check inputs
  if( length(x) ~= length(y) )
    error('x and y have different length');
  end
  if( mod(F,2)<0 ) 
    error('F must be odd');
  end

  Nxy = length(x);
  dx = diff(x);
  % Assuming equal spacing for the last two values
  dx(end+1) = dx(end);
  
  HalfWin  = ((F+1)/2)-1;

  xf = x;
  yf = zeros(Nxy, 1);

  % start and end values: calculate SavGol coefficients for each n
  for n = 0:(HalfWin-1)
    % first points
    c = savgolcoeff(n, F-1-n, ld, m);
    yf(n+1) = dot(c,y(1:F));

    % last points
    c = savgolcoeff(F-1-n, n, ld, m);
    yf(Nxy-n) = dot(c,y((Nxy-F+1):Nxy));
  end

  % all other values: calculate coefficients once, then perform 
  % the convolution
  c = savgolcoeff(HalfWin, HalfWin, ld, m);
  for n = ((F+1)/2):(Nxy-(F+1)/2+1)
    yf(n) = dot(c,y(n - HalfWin:n + HalfWin));
  end
  
  if(ld > 0)
    % Calculate derivative
    yf = yf./dx.^ld;
  end

end

