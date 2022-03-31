function [epsRMS,optOvl] = optOvlForceS(data,SGfiltparam,nodeParam,...
                                        S,p1,p2,p3)
% OPTOVLFORCES Optimises the overlap of force curves 
%   Each z-piezo axis in the given data-struct corresponds to a df(z)-curve
%   measured for a set voltage amplitude at a specific cycle time in the
%   measurement. This function shifts each z-piezo axis according to their
%   respective voltage setpoint by utilising the sensitivity factor S.
%   Additionally, a cycle time dependent, third degree polynomial composed
%   of the coefficents p1, p2 and p3 is substracted to compensate drift
%   effects. Subsequently a Sader-Jarvis Force deconvolution is carried out
%   based on the resulting new z-axes. This approach leads for matching
%   parameters S, p1, p2 and p3 to an ensemble of force curves in optOvl
%   overlapping each other in the scope of measurement accuracy. The epsRMS
%   value is a measure for the goodness of that overlap.
% 
%---------------------------------------------------  
%   Arguments:
%---------------------
%    data - data struct containing measurements for N different amplitude
%           setpoints with the fields:
%           
%    data.name : Name of measurement. Used for naming the output files.
%
%    data.Vamp : [1x2N double array], voltage amplitude setpoints (in V)   
% 
%    data.zp   : {1x2N cell}, absolute z-piezo positions for each 
%                data curve (in m)
% 
%    data.df   : {1x2N cell}, df data (in Hz) acquired 
%                with respect to data.zp 
%
%    data.TipParam : struct, containing the tip parameters for the 
%                measurement: 
% 
%        .TipParam.k0      : Spring constant (in N/m)
%        .TipParam.f0      : Eigenfrequency (in Hz)
%        .TipParam.Q       : Quality factor
%        .TipParam.SzGuess : Estimate for the Sz value, for example 
%                            measured via the gamma-method [2] (in m/V)  
%                                                   
%---------------------                                                
% 
%   SGfiltparam - Struct containing the parameters for the Savitzky-Golay
%                 filter [2], which is applied in the force deconvolution 
%                 using the Sader/Jarvis algorithm [3].
%                 Has to contain the following fields:
%
%   SGfiltparam.frame: frame length for SG filtering.
%                      sets the number of points for the polynomial fit
%                      (larger values correspond to stronger filtering)
%                                    
%   SGfiltparam.order: order of the fitted polynomial 
%                      (smaller values correspond to stronger filtering)
%
%---------------------
%
%   nodeParam - Struct containing all parameters for the node positioning
%               on the deconvoluted force curves. These nodes are
%               utilized by the FCA method for the evaluation of force
%               curve overlap goodness. This struct contains the fields
%                 
%   nodeParam.Find : [1xk positive integer array], indicies defining the 
%                    nodes on the attractive branch of the reference force
%                    curve counted from the smallest shared force
%                    index optOvl.IminFref.
%
%   nodeParam.avgRange : Integer, representing the single side range of 
%                        force averagig from the node position. The
%                        resulting averaging environment around the node is
%                        given by 2*avgRange+1.           
%
%---------------------
% 
%   S  - Sensitivity factor in m/V
%
%   p1 - Coefficient for the third degree of the cycle time in the      
%        drift compensation polynom
%
%   p2 - Coefficient for the second degree of the cycle time in       
%        the drift compensation polynom
%
%   p3 - Coefficient for the first degree of the cycle time in the      
%        drift compensation polynom
%
%------------------------------------------------------
%   Returns:
%---------------------
%   Note: j represents in the following the curve index and k the
%   respective node index. 
%   j ranges from 1 to 2N while k depends on the size of nodeParam.Find.
%       
%   epsRMS - Root mean square of all root mean squares of the 
%            distance slices for each Node through each curve. A
%            smaller value equals a better overlap of force curves.
%            This parameter is minimized by the fminsearch routine
%            during FCA optimisation
%
%---------------------
%
%   optOvl -  Struct containing the calculated force curves and other 
%             relevant variables used in this function. This struct 
%             contains the fields:
%
%   optOvl.resParam   : [SzResult p1 p2 p3] double, contains the 
%                        parameters for the resulting force curve overlap.
%                                   
%   optOvl.minI       : Index of the shortest force curve in the given 
%                        force curve ensemble.
%                         
%   optOvl.IminFref   : Index on the reference curve (j=1) to the smallest  
%                        absolute shared force of all force curves in the 
%                        ensemble. 
%
%
%   optOvl.eFslices   : [k x j] double, contains the average distances 
%                        along each slice of equal force across all force
%                        curves of the ensemble 1:j:2N.
%                 
%   optOvl.eFsliceRMS : [k x 1] double, summarises all average distance
%                        values in each slice of equal force into one Root
%                        Mean Square (RMS) value at position k.        
%
%   optOvl.epsRMS     : Root Mean Square of each Root Mean Square value of
%                        each slice. Saved in this struct as well, for 
%                        completeness.
%
%   optOvl.curve      : {1 x size(data.Vamp,2)} struct array, contains
%                        data for each calculated force curve in a 
%                        separate field. 
%                        Struct for each curve j contains the fields:                                  
%                 
%   optOvl.curve(j).Vamp : voltage amplitude setpoint (in V)
%
%   optOvl.curve(j).df   : df(zps)-curve (in Hz)
%
%   optOvl.curve(j).AS   : amplitude setpoint (in m), resulting from the 
%                           calculated sensitivity factor S and the voltage
%                           amplitude setpoint.              
%
%   optOvl.curve(j).ktsevencap : Cap-average of the force gradient 
%                                 (see also [4]). 
%
%   optOvl.curve(j).zps  : absolute z piezo axis, shifted (in m)
%
%   optOvl.curve(j).ztip : absolute z piezo axis shifted after 
%                           Sader-Jarvis force deconvolution (in m). 
%                           For the correct parameters S, p1, p2 and p3, 
%                           the ztip axis of each curve j is identical.
%                           (array is shorter than optOvl.curve(j).zps 
%                           due to the deconvolution)
%   
%   optOvl.curve(j).Ftip : Force (in N), calculated from ktsevencap and 
%                           zps. Filter applied according to SGfiltparam.
%                                     
%   optOvl.curve(j).minFtip  : Minimum force
%                               
%   optOvl.curve(j).IminFtip : Integer index minimum force value              
%
%   optOvl.curve(j).nodeI : [1 x k] integer, contains the k nodes 
%                            placed according to the given nodeParam               
%                            struct.               
%
%   optOvl.curve(j).nodeF : Force (in N) at the respective index nodeI 
%
%   optOvl.curve(j).nodeIenvi : [k x (2*nodeParam.avgEnvi+1)] integer
%                                index array, representing the averaging
%                                environment around each node k.
%                                                                                     
%   optOvl.curve(j).NtoRefN   : [1 x k] double, represents the distance
%                                (in m) between the average position 
%                                of node k on curve j and to the average 
%                                position of the same node on the 
%                                reference curve j=1.   
                                        
%   
% FCA method is based on
% [1] D. Heile, R. Olbrich, M. Reichling, P. Rahe
%     "Alignment method for the accurate and precise quantification of 
%     tip-surface forces"
%
% References: 
% [2] A. Savitzky, M.J.E. Golay, Anal. Chem. 36, 1627 (1964)
%     doi: 10.1021/ac60214a047
% [3] J.E. Sader and S.P. Jarvis, Appl. Phys. Lett. 84, 1801 (2004)
%     doi: 10.1063/1.1667267
% [4] H. Söngen, R. Bechstein, A. Kühnle, 
%     J. Phys. Cond. Matter 29, 274001 (2017)
%     doi: 10.1088/1361-648X/aa6f8b
%
% Copyright (C) 2020, Daniel Heile, Reinhard Olbrich, Michael Reichling,
% Philipp Rahe
%
% This is a script under the terms of the Creative Commons Attribution 
% License (creativecommons.org/licenses/by/4.0), which permits 
% unrestricted use, distribution, and reproduction in any medium, 
% provided the original work is properly cited.
%
% Last Version 13.08.2020, Daniel Heile (dheile@uos.de)


%% Initialise variables

    curve = struct();
    optOvl = struct();
    epsRMS = NaN; 


%% Z-Axis shift and subsequent force deconvolution 
    for j = 1:size(data.Vamp,2)           
        curve(j).Vamp = data.Vamp(1,j);
        curve(j).df=data.df{1,j};
        curve(j).AS = S*data.Vamp(1,j); %Calculate amplitude in m 

        % Calculate the Cap-average of the force gradient (see also [4])
        curve(j).ktsevencap = (data.TipParam.k0).*...
                              (1-(((data.TipParam.f0)+data.df{1,j})./...
                              (data.TipParam.f0)).^2);

        % Substract amplitude and drift correction from z piezo array. This
        % is the fundamental operation of the FCA method.
        curve(j).zps = data.zp{1,j}-curve(j).AS-...
                        (p1*((j-1)^3)+p2*((j-1)^2)+p3*(j-1));

        % Carry out the Sader-Jarvis force deconvolution [3]
        [curve(j).ztip,curve(j).Ftip] = ...
            Feven_deconv(curve(j).zps, curve(j).ktsevencap,...
                            curve(j).AS, SGfiltparam.frame, ...
                                SGfiltparam.order, false);

        % Find the index of minimim force in every curve. 
        [curve(j).minFtip,curve(j).IminFtip] = min(curve(j).Ftip); 
    end


%% Node placement

    % Find and save index of the force curve, containing the minimum lowest
    % force in the curve ensemble. optOvl.minI is indexing this 
    % (shortest) force curve.
    % Nodes will be placed with respect to the minimum force of the
    % shortest curve. This guarantees that all force curves are taken into
    % consideration during the optimisation process
    %
    [~, optOvl.minI]=max([curve(:).minFtip]);
    
    % Find that minimum force on the (drift free) reference curve. Searched
    % is only on the attractive (right) side of the minimum. This prevents
    % errors by doubles results, if the force curve extends beyond the
    % minimum into the repulsive regime.
    % Itmp is relative to the shorter reference curve considered in the
    % search above. Obtain the Index value relative to the whole reference
    % curve
    %
    [~,Itmp]=min(abs(curve(1).Ftip(curve(1).IminFtip:end)... 
                                -curve(optOvl.minI).minFtip));                      
    
    % This is the Index to the minimal force on the reference curve.
    % Relative to that the nodes are placed in the subsequent step.
    optOvl.IminFref=Itmp + curve(1).IminFtip;  
    clear Itmp
        
    % Place nodes at slices of equal force on the complete curve ensemble
    % relative to optOvl.IminFref utilizing the given node parameters
    %
    for k=1:size(nodeParam.Find,2) 
        % For each node:
        for j=1:size(data.Vamp,2)
            % Calculate the node points on every curve representing the
            % force slides. To prevent double placement, only the 
            % attractive (right) side of the curves is considered here.
            %
            [~,Itmp] = min(abs(curve(j).Ftip(curve(j).IminFtip:end)-... 
                curve(1).Ftip(optOvl.IminFref+nodeParam.Find(1,k))));
            % Now find the Index relative to the whole curve
            curve(j).nodeI(1,k) = Itmp + curve(j).IminFtip - 1;
            % and the associated force 
            curve(j).nodeF(1,k) = curve(j).Ftip(curve(j).nodeI(1,k));
            
            % Control if the node is placed on the right side of the
            % minimum
            %
            if curve(j).nodeI(1,k) < curve(j).IminFtip
                error(['Index : ' num2str(k) 'of curve nr.: ' num2str(j)...
                    'is not placed on the right side of the minimum']) 
            end            
        end
        % Control if the node placement on the reference curve equals the
        % given index Find(k)
        %
        if (optOvl.IminFref + nodeParam.Find(1,k)) ~= curve(1).nodeI(1,k)
          error(['Index ' num2str(k) ' is not placed on the expected '...
              'position on the reference curve'])
        end   
    end
    
    clear Itmp j k

    
%% Calculate epsRMS utilizing averaged node positions 


    for k = 1:size(nodeParam.Find,2) 
        for j = 1:size(data.Vamp,2) 
        % Calculate the index environment around each node for averaging
        %
        curve(j).nodeIenvi(k,:) = ...
            linspace(curve(j).nodeI(1,k) - nodeParam.avgRange,...
               curve(j).nodeI(1,k) + nodeParam.avgRange...
               ,2.*nodeParam.avgRange+1);
            
            % Control for possible index overflows respective to the
            % ztip-array size. For adequate data this should not be needed.
            for z = 1:size(curve(j).nodeIenvi,2)
                % Control for lower index overflow
                if curve(j).nodeIenvi(k,z) < 1
                 curve(j).nodeIenvi(k,z) = curve(j).nodeI(1,k);
                 disp([data.name ': Lower owerflow exception' newline ...
                        ' Correction of averaging index ' num2str(z) ...
                        ' around the node ' num2str(k) ...
                        ' at curve ' num2str(j)])
                end
                % Control for upper index overflow
                if curve(j).nodeIenvi(k,z) > size(curve(j).ztip) 
                 curve(j).nodeIenvi(k,z) = curve(j).nodeI(1,k);
                 disp([data.name ': Upper owerflow exception' newline ...
                        ' Correction of averaging index ' num2str(z) ...
                        ' around the node ' num2str(k) ...
                        ' at curve ' num2str(j)])
                end
            end
            
            % Calculate the distance to the reference node for every curve.
            % Here the average ztip position of each node (average of node
            % environment) is considered for the differences.
            % NtoRefN : Node to reference node distance (in m)
            curve(j).NtoRefN(k,1) = ...
                mean(curve(j).ztip(curve(j).nodeIenvi(k,:),1)) - ...
                    mean(curve(1).ztip(curve(1).nodeIenvi(k,:),1));
            % Thus each curve knows its own average distance to the
            % reference curve for each node.
            
            % Store the slices of average distances through each node of
            % equal forces in a general matrix.
            optOvl.eFslices(k,j) = curve(j).NtoRefN(k,1);            
        end
        % Calculate the Root Mean Square (RMS) of the distances for each
        % force slice k
        optOvl.eFsliceRMS(k,1) = ...
            sqrt(sum(optOvl.eFslices(k,:).^2)./size(data.Vamp,2)); 
    end
    
    % The Root Mean Square values of each distance slice for each node k
    % are condensed into another Root Mean Square value epsRMS
    epsRMS = sqrt((sum(optOvl.eFsliceRMS(:,1).^2))./size(nodeParam.Find,2));
    
 
%% Store variables in optOvl
    
    optOvl.resParam = [S p1 p2 p3];
    optOvl.epsRMS = epsRMS;
    optOvl.curve = curve;
end

