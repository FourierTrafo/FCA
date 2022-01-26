function [initParam,driftCorr] = findInitialParam(data)
%% FINDINITIALPARAM Finds inital parameters for the FCA method
%
% Arguments:
%---------------------
%   data - data struct as documented in the paramount function 
%          applyFCA(data_path,save_path,SGfiltparam,nodeParam,options)
%---------------------------------------------------------------------                  
% Returns:
%---------------------
%   initParam - [1x4] double, contains the intial parameter needed 
%               by the fminsearch function. It contains the following
%               values: [Sz0, p1, p2, p3] with
%               Sz0 : inital Sz value (in m/V)
%               p1  : 3rd degree drift coefficent (in m/s^3)
%               p2  : 2nd degree drift coefficent (in m/s^2)
%               p3  : 1st order drift coefficient (in m/s)
% 
%               (The zero degree coefficient of the poly3 fit is always
%               zero due to the offset subtraction before application of
%               the fit function)
% 
%   driftCorr - Struct containing all data used in the drift correction.
%               It contains the following fields:
%
%               driftCorr.tcycle   - Number of cycle time t/tcycle 
%                                    (thus an index) 
%
%               driftCorr.driftVal - Drift value array relative to zero
%
%               driftCorr.fit      - fit (poly3) to the drift array
%
%               driftCorr.gof      - Goodness of polynomial fit
%
%               driftCorr.Sz0      - The average Sz value calculated 
%                                    from all (zpiezo)/(Voltage amplitude)
%                                    difference quotients. Is used as the 
%                                    initial Sz value for the initParam 
%                                    array.

% FCA method is based on
% [1] D. Heile, R. Olbrich, M. Reichling, P. Rahe
%     "Alignment method for the accurate and precise quantification of 
%     tip-surface forces"
%
% Copyright (C) 2020, Daniel Heile, Reinhard Olbrich, Michael Reichling,
% Philipp Rahe
%
% This is a script under the terms of the Creative Commons Attribution 
% License (creativecommons.org/licenses/by/4.0), which permits 
% unrestricted use, distribution, and reproduction in any medium, 
% provided the original work is properly cited.
%
% version 13.08.2020, Daniel Heile (dheile@uos.de)

%% Calculate inital guess for S

    % Iteratively calculate array of S values from the difference quotient 
    % of the z-piezo position and their corresponding voltage amplitude
    % setpoints.
    %
    diffS = zeros(1,size(data.Vamp,2)-1);                  
    for i = 1:size(data.Vamp,2)-1                                                                                 
        diffS(1,i) = abs(data.zp{1,i+1}(1)-data.zp{1,i}(1))/...
                        abs((data.Vamp(1,i+1)-data.Vamp(1,i))); % in m/V     
    end
    
    % Calculate the average S value from the diffS array. This will
    % be the initial S value for the fminsearch routine.
    % Note: The difference quotient value between the the last ramp-up and 
    % the first ramp-down point is infinite due to equal voltage amplitude 
    % setpoints and thus is ignored in this procedure.
    %
    Sz0 = mean([diffS(1,1:(size(data.Vamp,2)/2-1)) ...
                diffS(1,(size(data.Vamp,2)/2+1):end)]); 


%% Calculate drift value array

    % Initial guess of the drift values for every z-piezeo position
    % corresponding to a voltage amplitde setpoint at the time tcycle(i).
    % The first curve is considered drift free in this calculation, and
    % thus is utilised as the relative reference at zero. Furthermore, the
    % previously calculated initial S value is utilised to substract the
    % intrinsic z-piezo differences resulting from the different voltage
    % amplitude setpoints.
    % These drift values are considered as function of the number of cycle 
    % time tcycle. 
    %
    driftVal = zeros(size(data.Vamp,2),1);
    for i = 2:size(data.Vamp,2)   
        driftVal(i,1) = data.zp{1,i}(1)-...
                            (Sz0*(data.Vamp(1,i)-data.Vamp(1,1)))-...
                                data.zp{1,1}(1);  
    end


%% Calculate drift compensation using a 'poly3' fit
    
    % Fit a third degree-polynomial to the previously calculated array 
    % driftVal(tcycle).
    %
    tcycle = (0:size(driftVal,1)-1)';
    [f,gof] = fit(tcycle,driftVal,'poly3');
    %
    % Note that the constant f.p4 is equal to zero. This is due to the fact
    % that the drift values were calculated relative to the data.zp{1,1}(1)
    % point, which was substracted from the complete array. That starting
    % point is considered drift free, hence the fit has to be zero at that
    % position.

    
%% Store return parameters 

    % The parameters for the drift corrections are saved in the drift
    % correction struct driftCorr
    %
    driftCorr.tcycle = tcycle;     % index number of cycle time 
    driftCorr.driftVal = driftVal; % drift value array 
    driftCorr.fit = f;             % poly3 fit to the drift array
    driftCorr.gof = gof;           % goodnes of that fit
    driftCorr.Sz0 = Sz0;     % mean sensitivity factor relative to z
    
    % Build initParam array for the fminsearch function. The first value 
    % has to be S in m/V. After that the parameters of the drift fit are
    % inserted consecutively beginning with the third degree coefficient.   
    % The last coefficient representing the offset is neglected, due to the
    % fact that it has to be zero for a drift free starting point of
    % measurement
    %
    initParam=[Sz0 f.p1 f.p2 f.p3];        

end

