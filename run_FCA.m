%% run_FCA Run to start FCA evaluation utilising the here set parameters
% Script to run the Force Curve Alignment (FCA) method for the data-files
% from the given path.
%
% Implements the FCA method as described in 
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
% Version 13.08.2020, Daniel Heile (dheile@uos.de)

clear all
close all
clc 


%% Define initial Parameters

%-------------------------- Filter Parameters -----------------------------
%
% Set the parameters for the Savitzky-Golay (SG) filter [2], applied 
% during the Sader-Jarvis (SJ) force deconvolution [3].
%
% References: 
% [2] A. Savitzky, M.J.E. Golay, Anal. Chem. 36, 1627 (1964)
%     doi: 10.1021/ac60214a047
% [3] J.E. Sader and S.P. Jarvis, Appl. Phys. Lett. 84, 1801 (2004)
%     doi: 10.1063/1.1667267
%
SGfiltparam = struct('order',[],'frame',[]); 

SGfiltparam.frame = 11; % frame length for SG filtering 
                        % sets the number of points for the polynomial fit
                        % (larger values correspond to stronger filtering)

SGfiltparam.order = 2;  % order of the fitted polynomial 
                        % (smaller values correspond to stronger filtering)


%-------------------------- Optimisation Nodes ----------------------------
%
% Parameters for the node placement on the reference force curve. Based on
% these positions the nodes on the other force curves will be placed at
% equal forces as on the reference curve. 
% To address noise the forces around each node position are averaged. The
% range can be adjusted below.
nodeParam = struct('Find',[],'avgRange',[],'avgEnvi',[]);

nodeParam.Find = [15 20 25 30 35 40]; % Index-distance counted up from the 
                                      % smallest shared force in the  
                                      % dataset on the attractive branch of
                                      % the reference curve. At these 
                                      % indices the nodes will be placed.

nodeParam.avgRange = 3; % Single side index range for the force averaging   
                        % around each node. Has to be an integer.
                        
% For more details see the documentation in the function
% FCAData = applyFCA(data_path,save_path,SGfiltparam,nodeParam,options)                        
                                                    
%-------------------------- Fminsearch options ----------------------------
%
% Set options for MATLAB's fminsearch function  
% MaxIter : Number of maximal allowed iterations of the fminsearch options
% TolFun  : Tolerance of the function value epsRMS in m
% TolX    : Tolerance of the parameter variation
% MaxFunEvals : Number of maximal allowed function evaluations
% Display : 'iter' : Show each iteration step of the fminsearch function.
%           'off'  : iterations steps not shown during calculation
% (more detailed information from the MATLAB help optimset)
options = optimset('MaxIter', 1000, 'TolFun', 1e-16, 'TolX', 1e-16,...
                'MaxFunEvals',2400,'Display','iter','PlotFcns'...
                ,@optimplotfval);


%-------------------------- Input data ------------------------------------
%
% Absolute path to the input data file (.mat file)
% The input data file has to follow the structure defined in the 
% supplementary information of [1]
% Here: example_data.mat file is used
cdir = dir('example_data.mat'); % ADJUST TO YOUR DATA SET
data_path= [cdir.folder '\' cdir.name]; 


%-------------------------- Output directory ------------------------------
%
% Output directory for storing the complete FCA analysis. The function 
% FCAData = applyFCA(data_path,save_path,SGfiltparam,nodeParam,options)
% will save its results and a copy of the given dataset to a subfolder in 
% the given output directory.
output_dir = 'exampleEval'; % ADJUST TO YOUR DATA SET
save_path = [cdir.folder filesep output_dir];
mkdir(save_path);


%% -------------------------- run FCA method ------------------------------
%
% Executes FCA by calling applyFCA
disp(['FCA method started for ' data_path ' ...'])
[FCAData] = applyFCA(data_path,save_path,SGfiltparam,nodeParam,options);     


