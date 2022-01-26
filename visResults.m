function visResults(data,FCAData,savevz)
%% VISRESULTS Visualises the results of FCA and compares to input data
%
%---------------------
%   Arguments:
%--------------------- 
% 
%   data - data struct as documented in the paramount function 
%          applyFCA(data_path,save_path,SGfiltparam,nodeParam,options)  
% 
%---------------------
%   
%   FCAData - Struct containg the data corrected by the FCA method. It
%             contains the following fields. For further information see
%             the paramount function 
%             applyFCA(data_path,save_path,SGfiltparam,nodeParam,options)
% 
%---------------------
% 
%    savevz = Path to the directory, where the here created figure will be
%             saved. Is generated in the paramount function  
%

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
% Version 13.08.2020, Daniel Heile (dheile@uos.de)


%% Calculate the force curves in the non-FCA corrected form.
% Apply the force calculation function to the data using 
% the estimated sensitivity and zero drift coefficients. 
[~,dataForces] = optOvlForceS(data,FCAData.corrProg.fminIn.SGfiltParam,...
    FCAData.corrProg.fminIn.nodeParam,data.TipParam.SzGuess,0,0,0);

  
%% Visualise the results and compare with the FCA result
fig =  figure;
set(fig, 'Units', 'normalized', 'Position', [0.05, 0.05, 1, 0.85])
mcols = winter(size(data.Vamp,2));

% Plot raw force curves from data
ax1 = subplot(1,2,1);
hold on
for i=1:size(data.Vamp,2)
   plot(dataForces.curve(i).ztip.*1e9,dataForces.curve(i).Ftip*1e12,...
       'Color',mcols(i,:),'DisplayName',num2str(round(data.Vamp(1,i),3))) 
end
xlabel('z_{tip} in nm');
ylabel('F_{tip} in pN');
title('Force curves without corrections')
ax1.YLim=[ax1.YLim(1) 1];
leg1=legend('Location','SouthEast');
title(leg1,['S_{z,\gamma} = ' num2str(round(data.TipParam.SzGuess.*1e9)) 'nm/V'...
    newline '\epsilon_{RMS} = ' num2str(dataForces.epsRMS,2) 'm'...
    newline newline 'V_{A} in V']);

% Plot the force curves using the FCA-method
ax2=subplot(1,2,2);
hold on
for i=1:size(FCAData.Vamp,2)
   plot(FCAData.ztip{1,i}.*1e9,FCAData.Ftip{1,i}.*1e12,...
       'Color',mcols(i,:),'DisplayName',...
       num2str(round(data.Vamp(1,i),3))) 
end
xlabel('z_{tip} in nm');
ylabel('F_{tip} in pN');
title('FCA result')
ax2.YLim=[ax2.YLim(1) 1];
leg2=legend('Location','SouthEast');
title(leg2,['S_{z,FCA} = ' num2str(round(FCAData.TipParam.SzResult.*1e9)) 'nm/V'... 
    newline '\epsilon_{RMS} = ' ...
    num2str(FCAData.corrProg.Fcalc.optOvl.epsRMS,2) 'm'...
    newline newline 'V_{A} in V']);

set(findall(fig,'-property','FontSize'),'FontSize',16);

set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,'FCA_Result_visualisation','-dpdf','-r0')

movefile('FCA_Result_visualisation.pdf',savevz)


%% Print results in command window

disp(' ')
disp(['FCA Results for: ' data.name])
disp(' ')
disp(['Measured SzGuess was ' num2str(data.TipParam.SzGuess.*1e9) ' nm/V'])
disp(' ')
disp(['SzFCA is ' num2str(FCAData.TipParam.SzResult.*1e9) ' nm/V'])
disp(['epsRMS ' num2str(FCAData.corrProg.Fcalc.epsRMS) ' m'])


end

