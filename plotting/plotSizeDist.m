% Script to plot any number of size distribution variables
% Written by Dan Stechman
%
% Requires that size distribution data be processed through the 
% avgSizeDists script

close all;clearvars;

flight = '20150617';

probe = 'CIP';

outFileAppend = '';

plotAvg = 0; % Plot averaged size distribution data (as defined in sDistFile)
avgTime = 10; % Averaging time - used to determine which data file to pull in - only used if plotAvg is False

saveFigs    = 1;
noDisp      = 1;
% Ftype		= '-dpdf';
Ftype		= '-dpng';

if strcmp(Ftype,'-dpdf')
	Fres = '-painters'; % Use this for fully vectorized files
else
	Fres = '-r0'; % Use this for smaller files - saves figure with same resolution/size as displayed on screen
end

allSpirals = 1; % Plot every spiral for given flight - otherwise specified spirals given below are plotted
showMarkers = 0; % Plot markers at given (below) locations on line plots

savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';
dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

plotND				= 0;
plotMD				= 0;

plotNDtime			= 0;
plotMDtime			= 0;
plotNDtemp			= 0;
plotMDtemp			= 0;

plotNDtempBinned	= 0;
plotMDtempBinned	= 1;

plotNtTemp			= 0;
plotNtTempAll		= 0;
plotTWCtemp			= 0;
plotTWCtempAll		= 0;
plotDmmTemp			= 0;
plotDmmTempAll		= 0;

cntrLine = 'Median'; % Do we plot median or mean on spread plots?
% cntrLine = 'Mean';
plotNtTempSprd		= 0;
plotTWCtempSprd		= 0;
plotDmmTempSprd		= 0;

plotNtTempSprdF		= 0;
plotTWCTempSprdF	= 0;
plotDmmTempSprdF	= 0;


startT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'startT');
mcsStg = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mcsStg');
sprlZone = nc_varget([dataPath '/' flight '_PECANparams.nc'],'sprlZone');
mlTopTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlTopTime');
mlTopTemp = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlTopTemp');
mlBotTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlBotTime');
mlBotTemp = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlBotTemp');

sDistFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.' probe '.' num2str(avgTime) 'secAvg' outFileAppend '.mat'];

FLfile = [dataPath 'FlightLevelData/Processed/' flight '_FltLvl_Processed.mat'];

load(sDistFile,'-regexp', '^(?!outFileAppend$)\w');
load(FLfile, '-regexp', '^(?!flight$|dataPath$|lon$|lat$)\w');

sprlNames = fieldnames(time_secs_orig); % Variable used unimportant - just needs to be one of the structs

if ~allSpirals
	switch flight
		case '20150617'
			loopVctr = [3,4];
		case '20150706'
			loopVctr = [6,7]; % Spiral numbers to plot
			if showMarkers
				loopSDix = [28, 3]; %Index of 10 sec avg'd SD data where marker is desired
			end
			
		case '20150709'
			loopVctr = [2, 5, 9];
			if showMarkers
				loopSDix = [23, 2, 72]; %10 sec avg
			end
	end
else
	loopVctr = 1:length(sprlNames);
end

%% Sort spirals based on spiral zone and MCS evolution
tzSprls = find(sprlZone == 'T');
esrSprls = find(sprlZone == 'S');
raSprls = find(sprlZone == 'A');

formSprls = find(mcsStg == 'F');
matureSprls = find(mcsStg == 'M');
weakSprls = find(mcsStg == 'W');

%% Transpose select sDist variables
for i=1:length(sprlNames)
	conc_minR_avg.(sprlNames{i}) = conc_minR_avg.(sprlNames{i})';
	n_avg.(sprlNames{i}) = n_avg.(sprlNames{i})';
	
	conc_minR_orig.(sprlNames{i}) = conc_minR_orig.(sprlNames{i})';
	n_orig.(sprlNames{i}) = n_orig.(sprlNames{i})';
end


%% Create directories to save plots in if they don't already exist
if saveFigs
    saveDir = [savePath flight];
	if (exist(saveDir, 'dir') ~= 7)
        mkdir(saveDir)
	end
	if ~plotAvg
		if (plotND && exist([saveDir '/' probe '-ND_1s'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-ND_1s'])
		end
		if (plotMD && exist([saveDir '/' probe '-MD_1s'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-MD_1s'])
		end
		if (plotNDtime && exist([saveDir '/' probe '-ND-Time_1s'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-ND-Time_1s'])
		end
		if (plotMDtime && exist([saveDir '/' probe '-MD-Time_1s'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-MD-Time_1s'])
		end
		if (plotNDtemp && exist([saveDir '/' probe '-ND-Temp_1s'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-ND-Temp_1s'])
		end
		if (plotMDtemp && exist([saveDir '/' probe '-MD-Temp_1s'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-MD-Temp_1s'])
		end
		if (plotNDtempBinned && exist([saveDir '/' probe '-ND-TempBinned'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-ND-TempBinned'])
		end
		if (plotMDtempBinned && exist([saveDir '/' probe '-MD-TempBinned'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-MD-TempBinned'])
		end
		if ( (plotNtTemp || plotNtTempAll || plotNtTempSprd || plotNtTempSprdF) && exist([saveDir '/' probe '-Nt-Temp_1s'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-Nt-Temp_1s'])
		end
		if ( (plotTWCtemp || plotTWCtempAll || plotTWCtempSprd || plotTWCTempSprdF) && exist([saveDir '/' probe '-TWC-Temp_1s'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-TWC-Temp_1s'])
		end
		if ( (plotDmmTemp || plotDmmTempAll || plotDmmTempSprd || plotDmmTempSprdF) && exist([saveDir '/' probe '-Dmm-Temp_1s'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-Dmm-Temp_1s'])
		end
	else
		if (plotND && exist([saveDir '/' probe '-ND_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-ND_' num2str(avgTime) 's'])
		end
		if (plotMD && exist([saveDir '/' probe '-MD_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-MD_' num2str(avgTime) 's'])
		end
		if (plotNDtime && exist([saveDir '/' probe '-ND-Time_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-ND-Time_' num2str(avgTime) 's'])
		end
		if (plotMDtime && exist([saveDir '/' probe '-MD-Time_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-MD-Time_' num2str(avgTime) 's'])
		end
		if (plotNDtemp && exist([saveDir '/' probe '-ND-Temp_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-ND-Temp_' num2str(avgTime) 's'])
		end
		if (plotMDtemp && exist([saveDir '/' probe '-MD-Temp_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-MD-Temp_' num2str(avgTime) 's'])
		end
		if (plotNDtempBinned && exist([saveDir '/' probe '-ND-TempBinned'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-ND-TempBinned'])
		end
		if (plotMDtempBinned && exist([saveDir '/' probe '-MD-TempBinned'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-MD-TempBinned'])
		end
		if ( (plotNtTemp || plotNtTempAll || plotNtTempSprd || plotNtTempSprdF) && exist([saveDir '/' probe '-Nt-Temp_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-Nt-Temp_' num2str(avgTime) 's'])
		end
		if ( (plotTWCtemp || plotTWCtempAll || plotTWCtempSprd || plotTWCTempSprdF) && exist([saveDir '/' probe '-TWC-Temp_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-TWC-Temp_' num2str(avgTime) 's'])
		end
		if ( (plotDmmTemp || plotDmmTempAll || plotDmmTempSprd || plotDmmTempSprdF) && exist([saveDir '/' probe '-Dmm-Temp_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/' probe '-Dmm-Temp_' num2str(avgTime) 's'])
		end
	end
end

%% Specify flight and/or probe specific plotting parameters

% Default values for probes - can be changed flight-by-flight in the proceding switch
% If any limit vector below is set to be empty (except NDLogLim), auto limits are applied
if strcmp(probe,'CIP')
	diamLim = [0 2];
elseif strcmp(probe,'PIP')
	diamLim = [0 7];
end

tempRangeAll = [-18.5 22];

switch flight
	case '20150617'
		if strcmp(probe,'CIP')
			NDLim = [1e-4 10];
			NDtempBinLim = [5e-5 40];
			MDtempBinLim = [];
			NDLogLim = [-4, 2];
			MDLim = [1e-8 8e-6];
			MDLogLim = [-8, -4];
			NtRangeAll = [4e-6 100];
			TWCrangeAll = [2.5e-6 20];
			DmmRangeAll = [-0.1 1.8];
			NtRangeFill = [-4 2];
			TWCRangeFill = [-4 0.5];
		elseif strcmp(probe,'PIP')
			NDLim = [];
			NDLogLim = [-5, 1];
			MDLim = [1e-10 1e-5];
			NtRangeAll = [4.8e-5 10];
			TWCrangeAll = [2.3e-6 1];
			DmmRangeAll = [-0.1 6.5];
			NtRangeFill = [];
			TWCRangeFill = [];
		end
	case '20150620'
		if strcmp(probe,'CIP')
			NDLim = [1e-3 10];
			NDtempBinLim = [5e-5 30];
			MDtempBinLim = [];
			NDLogLim = [-4, 2];
			MDLim = [1e-7 8e-6];
			MDLogLim = [-8, -4];
			NtRangeAll = [4e-6 100];
			TWCrangeAll = [2.5e-6 20];
			DmmRangeAll = [-0.1 1.8];
			NtRangeFill = [-4 2];
			TWCRangeFill = [-4 0.5];
		elseif strcmp(probe,'PIP')
			NDLim = [];
			NDLogLim = [-5, 1];
			MDLim = [1e-10 1e-5];
			NtRangeAll = [4.8e-5 10];
			TWCrangeAll = [2.3e-6 1];
			DmmRangeAll = [-0.1 6.5];
			NtRangeFill = [];
			TWCRangeFill = [];
		end
	case '20150701'
		if strcmp(probe,'CIP')
			NDLim = [];
			NDtempBinLim = [9e-5 10];
			MDtempBinLim = [];
			NDLogLim = [-4, 2];
			MDLim = [];
			MDLogLim = [-8, -4];
			NtRangeAll = [4e-6 100];
			TWCrangeAll = [2.5e-6 20];
			DmmRangeAll = [-0.1 1.8];
			NtRangeFill = [-4 2];
			TWCRangeFill = [-4 0.5];
		elseif strcmp(probe,'PIP')
			NDLim = [];
			NDLogLim = [-5, 1];
			MDLim = [1e-10 1e-5];
			NtRangeAll = [4.8e-5 10];
			TWCrangeAll = [2.3e-6 1];
			DmmRangeAll = [-0.1 6.5];
			NtRangeFill = [];
			TWCRangeFill = [];
		end
	case '20150702'
		if strcmp(probe,'CIP')
			NDLim = [3e-5 0.2];
			NDtempBinLim = [5e-5 10];
			MDtempBinLim = [];
			NDLogLim = [-4, 2];
			MDLim = [2e-9 5e-7];
			MDLogLim = [-8, -4];
			NtRangeAll = [4e-6 100];
			TWCrangeAll = [2.5e-6 20];
			DmmRangeAll = [-0.1 1.8];
			NtRangeFill = [-4 2];
			TWCRangeFill = [-4 0.5];
		end
	case '20150706'
		if strcmp(probe,'CIP')
			NDLim = [3e-4 10];
			NDtempBinLim = [4e-5 20];
			MDtempBinLim = [];
			NDLogLim = [-4, 2];
			MDLim = [7e-8 4e-6];
			MDLogLim = [-8, -4];
			NtRangeAll = [4e-6 100];
			TWCrangeAll = [2.5e-6 20];
			DmmRangeAll = [-0.1 1.8];
			NtRangeFill = [-4 2];
			TWCRangeFill = [-4 0.5];
		elseif strcmp(probe,'PIP')
			NDLim = [1e-6 1];
			NDLogLim = [-5, 1];
			MDLim = [1e-10 1e-5];
			NtRangeAll = [4.8e-5 10];
			TWCrangeAll = [2.3e-6 1];
			DmmRangeAll = [-0.1 6.5];
			NtRangeFill = [];
			TWCRangeFill = [];
		end
	case '20150709'
		if strcmp(probe,'CIP')
			NDLim = [4e-4 0.5];
			NDtempBinLim = [6e-5 10];
			MDtempBinLim = [];
			NDLogLim = [-4, 2];
			MDLim = [4e-8 2e-6];
			MDLogLim = [-8, -4];
			NtRangeAll = [4e-6 100];
			TWCrangeAll = [2.5e-6 20];
			DmmRangeAll = [-0.1 1.8];
			NtRangeFill = [-4 2];
			TWCRangeFill = [-4 0.5];
		elseif strcmp(probe,'PIP')
			NDLim = [1e-6 1];
			NDLogLim = [-5, 1];
			MDLim = [1e-10 1e-5];
			NtRangeAll = [];
			TWCrangeAll = [];
			DmmRangeAll = [-0.1 6.5];
			NtRangeFill = [];
			TWCRangeFill = [];
		end
end

%% Plotting
if plotND
    for ix = loopVctr
		
		if plotAvg
			conc_minR = conc_minR_avg.(sprlNames{ix});
		else
			conc_minR = conc_minR_orig.(sprlNames{ix});
		end
		

		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		stairs(bin_min, nanmean(conc_minR,2), 'b', 'LineWidth', 2);
		
		if ~plotAvg
			title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - 1s Avg']);
		else
			title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - ' num2str(avgTime) 's Avg']);
		end
		
		
		xlabel('D [mm]');
		ylabel('N(D) [cm^{-4}]');
		set(gca,'Yscale','log');
		if ~isempty(diamLim)
			xlim(diamLim);
		end
		if ~isempty(NDLim)
			ylim(NDLim);
		end
		set(gca,'XMinorTick','on','YMinorTick','on');
		set(findall(gcf,'-property','FontSize'),'FontSize',28)
		grid


        if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			if ~plotAvg
				print([saveDir '/' probe '-ND_1s/' flight '_' probe '_ND_1s_S' num2str(ix) outFileAppend],Ftype,Fres)
			else
				print([saveDir '/' probe '-ND_' num2str(avgTime) 's/' flight '_' probe '_ND_' num2str(avgTime) 's_S' num2str(ix) outFileAppend],Ftype,Fres)
			end
        end
    end
end

if plotMD
    for ix = loopVctr
		
		if plotAvg
			mass_twc = mass_twc_avg.(sprlNames{ix});
		else
			mass_twc = mass_twc_orig.(sprlNames{ix});
		end
		
		meanMass = nanmean(mass_twc,1);
		
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		yyaxis left
		stairs(bin_min, meanMass, 'LineWidth', 2);
		
		if ~plotAvg
			title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - 1s Avg']);
		else
			title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - ' num2str(avgTime) 's Avg']);
		end
		
		xlabel('D [mm]');
		ylabel('M_{twc}(D) [g cm^{-4}]');
		set(gca,'Yscale','log');
		if ~isempty(diamLim)
			xlim(diamLim);
		end
		if ~isempty(MDLim)
			ylim(MDLim);
		end
		set(gca,'XMinorTick','on','YMinorTick','on');
		set(findall(gcf,'-property','FontSize'),'FontSize',28)
		grid;
		
		massCDF = zeros(size(meanMass));
		sumMeanMass = nansum(meanMass);
		for ii=1:length(meanMass)
			massCDF(ii) = nansum(meanMass(1:ii))/sumMeanMass;
		end
		
		yyaxis right
		plot(bin_min,massCDF*100,'LineWidth',2);
		ylabel('Mass_{twc}(D) CDF [%]');
		
		
        if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			if ~plotAvg
				print([saveDir '/' probe '-MD_1s/' flight '_' probe '_MD_1s_S' num2str(ix) outFileAppend],Ftype,Fres)
			else
				print([saveDir '/' probe '-MD_' num2str(avgTime) 's/' flight '_' probe '_MD_' num2str(avgTime) 's_S' num2str(ix) outFileAppend],Ftype,Fres)
			end
        end
    end
end

if plotNDtime
    for ix = loopVctr
        
		if plotAvg
			conc_minR = conc_minR_avg.(sprlNames{ix});
			time_secs = time_secs_avg.(sprlNames{ix});
		else
			conc_minR = conc_minR_orig.(sprlNames{ix});
			time_secs = time_secs_orig.(sprlNames{ix});
		end
		
        if saveFigs && noDisp
            figure('visible','off','Position', [10,10,1920,700]);
        else
            figure('Position', [10,10,1920,700]);
		end
			
		contourf(time_secs/24/3600,bin_mid,log10(conc_minR),NDLogLim(1):0.1:NDLogLim(2),'LineColor','none');
		
		if ~plotAvg
			title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - 1s Avg']);
		else
			title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - ' num2str(avgTime) 's Avg']);
		end
		
		
		ylabel('D [mm]');
		datetick('x','HH:MM:SS');
		set(gca,'XMinorTick','on','YMinorTick','on','XTickLabelRotation',45);
		colormap(jetmod); %Uses modified 'jet' colormap
		c=colorbar;
		ylabel(c,'log_{10}N(D) [cm^{-4}]');
		set(findall(gcf,'-property','FontSize'),'FontSize',28)
		set(gca, 'CLim', NDLogLim);

        
        if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			if ~plotAvg
				print([saveDir '/' probe '-ND-Time_1s/' flight '_' probe '_NDtime_1s_S' num2str(ix) outFileAppend],Ftype,Fres)
			else
				print([saveDir '/' probe '-ND-Time_' num2str(avgTime) 's/' flight '_' probe '_NDtime_' num2str(avgTime) 's_S' num2str(ix) outFileAppend],Ftype,Fres)
			end
        end
    end
end

if plotMDtime
    for ix = loopVctr
		if plotAvg
			mass_twc = mass_twc_avg.(sprlNames{ix});
			time_secs = time_secs_avg.(sprlNames{ix});
		else
			mass_twc = mass_twc_orig.(sprlNames{ix});
			time_secs = time_secs_orig.(sprlNames{ix});
		end
		
        if saveFigs && noDisp
            figure('visible','off','Position', [10,10,1920,700]);
        else
            figure('Position', [10,10,1920,700]);
		end
			
		contourf(time_secs/24/3600,bin_mid,log10(mass_twc'),MDLogLim(1):0.1:MDLogLim(2),'LineColor','none');
		
		if ~plotAvg
			title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - 1s Avg']);
		else
			title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - ' num2str(avgTime) 's Avg']);
		end
		
		
		ylabel('D [mm]');
		datetick('x','HH:MM:SS');
		set(gca,'XMinorTick','on','YMinorTick','on','XTickLabelRotation',45);
		colormap(jetmod); %Uses modified 'jet' colormap
		c=colorbar;
		ylabel(c,'log_{10}Mass_{twc}(D) [g cm^{-4}]');
		set(findall(gcf,'-property','FontSize'),'FontSize',28)
		set(gca, 'CLim', MDLogLim);

        
        if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			if ~plotAvg
				print([saveDir '/' probe '-MD-Time_1s/' flight '_' probe '_MDtime_1s_S' num2str(ix) outFileAppend],Ftype,Fres)
			else
				print([saveDir '/' probe '-MD-Time_' num2str(avgTime) 's/' flight '_' probe '_MDtime_' num2str(avgTime) 's_S' num2str(ix) outFileAppend],Ftype,Fres)
			end
        end
    end
end

if plotNDtemp
    for ix = loopVctr
		
		if plotAvg
			conc_minR = conc_minR_avg.(sprlNames{ix})';
			time_secs = time_secs_avg.(sprlNames{ix});	
		else
			conc_minR = conc_minR_orig.(sprlNames{ix})';
			time_secs = time_secs_orig.(sprlNames{ix});
		end
		
		tempCsprl = tempC_orig.(sprlNames{ix});
		time_fl = time_secsFL_orig.(sprlNames{ix});
		

        if saveFigs && noDisp
            figure('visible','off','Position', [10,10,1500,1000]);
        else
            figure('Position', [10,10,1500,1000]);
		end
		
		set(gcf,'defaultAxesColorOrder',[0 0 0; 0 0 0]);

        yyaxis left
		ax = gca;
        contourf(bin_mid,time_secs/24/3600,log10(conc_minR),NDLogLim(1):0.1:NDLogLim(2),'LineColor','none');
        xlabel('D [mm]');
        ylabel('Time');
		set(ax,'XMinorTick','on');
		colormap(jetmod);
        c=colorbar;
		set(c,'Location','southoutside');
		ylabel(c,'log_{10}N(D) [cm^{-4}]');
		set(ax, 'CLim', NDLogLim);
		
		% Plot ML top/bottom locations and annotate with temp
		hold on
		plot([0 1.9], [1 1]*mlTopTime(ix)/24/3600,'k--')
		topStr = sprintf('%.3f %cC',mlTopTemp(ix),char(176));
		tMT = text(0.1,mlTopTime(ix)/24/3600,topStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
		plot([0 1.9], [1 1]*mlBotTime(ix)/24/3600,'k--')
		botStr = sprintf('%.3f %cC',mlBotTemp(ix),char(176));
		tMB = text(0.1,mlBotTime(ix)/24/3600,botStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
		hold off

		% Swap time direction depending on spiral direction and define top/bottom spiral temps
		if tempCsprl(1) < tempCsprl(end)
			set(ax,'YDir','reverse');
			topT = floor(tempCsprl(1));
			botT = ceil(tempCsprl(end));
		else
			topT = floor(tempCsprl(end));
			botT = ceil(tempCsprl(1));
		end
		
		% Set the tick frequency for the time axis and adjust plot accordingly
		delta = 60; % 1 minute
		dtStrt = (floor(time_secs(1)/delta)*delta)/24/3600;
		dtEnd = time_secs(end)/24/3600;
		dtDelta = delta/24/3600; 
		set(ax,'YTick',dtStrt:dtDelta:dtEnd);
		datetick('y','HH:MM','keepticks','keeplimits');
		
		
        yyaxis right
		
		dummyX = zeros(size(time_secs));
		
		if abs(topT-botT) >= 30
			yTemp = botT:-4:topT;
		else
			yTemp = botT:-2:topT;
		end
		
		index = NaN(size(yTemp));
		yLbl = cell(size(yTemp));
		for iT = 1:length(yTemp)
			yLbl{iT} = num2str(yTemp(iT));
			[tmpMin,index(iT)] = min(abs(tempCsprl-yTemp(iT)));
			% fprintf('Deviation from desired %.1f: %.2f\n',yTemp(iT),tmpMin);
			if (index(iT) == 1 ||  index(iT) == length(tempCsprl) || tmpMin > 0.3 )
				index(iT) = NaN;
			end
		end
		indexFinal = index(~isnan(index));
		time_secs_s = time_secs/24/3600;
		
		if(tempCsprl(1) < tempCsprl(end))
			
			plot(dummyX,time_secs_s,'Color','w');
			yLblFinal = flip(yLbl(~isnan(index)));
			ax.YTick = flip(time_secs_s(indexFinal));
			set(gca,'YDir','reverse');
			ax.YTickLabel = yLblFinal;
		else
			plot(dummyX,time_secs_s,'Color','w');
			yLblFinal = yLbl(~isnan(index));
			ax.YTick = time_secs_s(indexFinal);
			ax.YTickLabel = yLblFinal;
		end
		
		ylabel(sprintf('T (%cC)', char(176)));
		
		
		if ~plotAvg
				title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - 1s Avg']);
		else
			title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - ' num2str(avgTime) 's Avg']);
		end
	

		set(findall(gcf,'-property','FontSize'),'FontSize',26)
		set(tMB,'FontSize',14);
		set(tMT,'FontSize',14);
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			if ~plotAvg
				print([saveDir '/' probe '-ND-Temp_1s/' flight '_' probe '_ND-Temp_1s_S' num2str(ix) outFileAppend],Ftype,Fres)
			else
				print([saveDir '/' probe '-ND-Temp_' num2str(avgTime) 's/' flight '_' probe '_ND-Temp_' num2str(avgTime) 's_S' num2str(ix) outFileAppend],Ftype,Fres)
			end
		end
    end
end

if plotMDtemp
    for ix = loopVctr
		
		if plotAvg
			mass_twc = mass_twc_avg.(sprlNames{ix});
			time_secs = time_secs_avg.(sprlNames{ix});
		else
			mass_twc = mass_twc_orig.(sprlNames{ix});
			time_secs = time_secs_orig.(sprlNames{ix});
		end
		
		tempCsprl = tempC_orig.(sprlNames{ix});
		time_fl = time_secsFL_orig.(sprlNames{ix});
		

        if saveFigs && noDisp
            figure('visible','off','Position', [10,10,1500,1000]);
        else
            figure('Position', [10,10,1500,1000]);
		end
		
		set(gcf,'defaultAxesColorOrder',[0 0 0; 0 0 0]);

        yyaxis left
		ax = gca;
        contourf(bin_mid,time_secs/24/3600,log10(mass_twc),MDLogLim(1):0.1:MDLogLim(2),'LineColor','none');
        xlabel('D [mm]');
        ylabel('Time');
		set(ax,'XMinorTick','on');
		colormap(jetmod);
        c=colorbar;
		set(c,'Location','southoutside');
		ylabel(c,'log_{10}Mass_{TWC}M(D) [g cm^{-4}]');
		set(ax, 'CLim', MDLogLim);
		
		% Plot ML top/bottom locations and annotate with temp
		hold on
		plot([0 1.9], [1 1]*mlTopTime(ix)/24/3600,'k--')
		topStr = sprintf('%.3f %cC',mlTopTemp(ix),char(176));
		tMT = text(0.1,mlTopTime(ix)/24/3600,topStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
		plot([0 1.9], [1 1]*mlBotTime(ix)/24/3600,'k--')
		botStr = sprintf('%.3f %cC',mlBotTemp(ix),char(176));
		tMB = text(0.1,mlBotTime(ix)/24/3600,botStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
		hold off

		% Swap time direction depending on spiral direction and define top/bottom spiral temps
		if tempCsprl(1) < tempCsprl(end)
			set(ax,'YDir','reverse');
			topT = floor(tempCsprl(1));
			botT = ceil(tempCsprl(end));
		else
			topT = floor(tempCsprl(end));
			botT = ceil(tempCsprl(1));
		end
		
		% Set the tick frequency for the time axis and adjust plot accordingly
		delta = 60; % 1 minute
		dtStrt = (floor(time_secs(1)/delta)*delta)/24/3600;
		dtEnd = time_secs(end)/24/3600;
		dtDelta = delta/24/3600; 
		set(ax,'YTick',dtStrt:dtDelta:dtEnd);
		datetick('y','HH:MM','keepticks','keeplimits');
		
		
        yyaxis right
		
		dummyX = zeros(size(time_secs));
		
		if abs(topT-botT) >= 30
			yTemp = botT:-4:topT;
		else
			yTemp = botT:-2:topT;
		end
		
		index = NaN(size(yTemp));
		yLbl = cell(size(yTemp));
		for iT = 1:length(yTemp)
			yLbl{iT} = num2str(yTemp(iT));
			[tmpMin,index(iT)] = min(abs(tempCsprl-yTemp(iT)));
			% fprintf('Deviation from desired %.1f: %.2f\n',yTemp(iT),tmpMin);
			if (index(iT) == 1 ||  index(iT) == length(tempCsprl) || tmpMin > 0.3 )
				index(iT) = NaN;
			end
		end
		indexFinal = index(~isnan(index));
		time_secs_s = time_secs/24/3600;
		
		if(tempCsprl(1) < tempCsprl(end))
			
			plot(dummyX,time_secs_s,'Color','w');
			yLblFinal = flip(yLbl(~isnan(index)));
			ax.YTick = flip(time_secs_s(indexFinal));
			set(gca,'YDir','reverse');
			ax.YTickLabel = yLblFinal;
		else
			plot(dummyX,time_secs_s,'Color','w');
			yLblFinal = yLbl(~isnan(index));
			ax.YTick = time_secs_s(indexFinal);
			ax.YTickLabel = yLblFinal;
		end
		
		ylabel(sprintf('T (%cC)', char(176)));
		
		
		if ~plotAvg
				title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - 1s Avg']);
		else
			title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - ' num2str(avgTime) 's Avg']);
		end
	

		set(findall(gcf,'-property','FontSize'),'FontSize',26)
		set(tMB,'FontSize',14);
		set(tMT,'FontSize',14);
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			if ~plotAvg
				print([saveDir '/' probe '-MD-Temp_1s/' flight '_' probe '_MD-Temp_1s_S' num2str(ix) outFileAppend],Ftype,Fres)
			else
				print([saveDir '/' probe '-MD-Temp_' num2str(avgTime) 's/' flight '_' probe '_MD-Temp_' num2str(avgTime) 's_S' num2str(ix) outFileAppend],Ftype,Fres)
			end
		end
    end
end

if plotNDtempBinned
% 	maxConc = [];
% 	minConc = [];
	for ix = loopVctr
		
		conc_minR_whole = conc_minR_orig.(sprlNames{ix});
		time_secs_whole = time_secs_orig.(sprlNames{ix});
		
		tempCsprl_whole = tempC_orig.(sprlNames{ix});
		
		tempCsprl = ones(length(tempCsprl_whole),1)*NaN;
		for ii=1:length(tempCsprl_whole)
			if tempCsprl_whole(ii) < 0
				tempCsprl(ii) = ceil(tempCsprl_whole(ii));
			elseif tempCsprl_whole(ii) > 0
				tempCsprl(ii) = floor(tempCsprl_whole(ii));
			end
		end
		
		for iii=min(tempCsprl):max(tempCsprl)
			conc_minR = conc_minR_whole(:,tempCsprl == iii);
			time_secs = time_secs_whole(tempCsprl == iii);
			
		
			if saveFigs && noDisp
				figure('visible','off','Position', [10,10,1500,1000]);
			else
				figure('Position', [10,10,1500,1000]);
			end

% 			tempMean = nanmean(conc_minR,2);
% 			maxConc = [maxConc; nanmax(nanmax(tempMean(tempMean>0)))];
% 			minConc = [minConc; nanmin(nanmin(tempMean(tempMean>0)))];
			stairs(bin_min, nanmean(conc_minR,2), 'r', 'LineWidth', 2);

			title(sprintf('%s - Spiral %d - %s - %d%cC avg',flight,ix,probe,iii,char(176)));

			xlabel('D [mm]');
			ylabel('N(D) [cm^{-4}]');
			set(gca,'Yscale','log');
			if ~isempty(diamLim)
				xlim(diamLim);
			end
			if ~isempty(NDtempBinLim)
				ylim(NDtempBinLim);	
			end
			set(gca,'XMinorTick','on','YMinorTick','on');
			set(findall(gcf,'-property','FontSize'),'FontSize',28)
			grid


			if saveFigs
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				print([saveDir '/' probe '-ND-TempBinned/' flight '_' probe '_ND-Temp_S' num2str(ix) '_' num2str(iii) 'degC' outFileAppend],Ftype,Fres)

			end
		end
	end
end

if plotMDtempBinned
	for ix = loopVctr
		mass_twc_whole = mass_twc_orig.(sprlNames{ix});

		time_secs_whole = time_secs_orig.(sprlNames{ix});
		
		tempCsprl_whole = tempC_orig.(sprlNames{ix});
		
		tempCsprl = ones(length(tempCsprl_whole),1)*NaN;
		for ii=1:length(tempCsprl_whole)
			if tempCsprl_whole(ii) < 0
				tempCsprl(ii) = ceil(tempCsprl_whole(ii));
			elseif tempCsprl_whole(ii) > 0
				tempCsprl(ii) = floor(tempCsprl_whole(ii));
			end
		end
		
		for iii=min(tempCsprl):max(tempCsprl)
			mass_twc = mass_twc_whole(tempCsprl == iii,:);
			time_secs = time_secs_whole(tempCsprl == iii);
			
		
			meanMass = nanmean(mass_twc,1);
			
			if saveFigs && noDisp
				figure('visible','off','Position', [10,10,1500,1000]);
			else
				figure('Position', [10,10,1500,1000]);
			end
			
			yyaxis left
			stairs(bin_min, meanMass, 'b', 'LineWidth', 2);

			title(sprintf('%s - Spiral %d - %s - %d%cC avg',flight,ix,probe,iii,char(176)));

			xlabel('D [mm]');
			ylabel('M_{twc}(D) [g cm^{-4}]');
			set(gca,'Yscale','log');
			if ~isempty(diamLim)
				xlim(diamLim);
			end
			if ~isempty(MDtempBinLim)
				ylim(MDtempBinLim);	
			end
			set(gca,'XMinorTick','on','YMinorTick','on');
			set(findall(gcf,'-property','FontSize'),'FontSize',28)
			grid
			
			massCDF = zeros(size(meanMass));
			sumMeanMass = nansum(meanMass);
			for ii=1:length(meanMass)
				massCDF(ii) = nansum(meanMass(1:ii))/sumMeanMass;
			end
		
			yyaxis right
			plot(bin_min,massCDF*100,'LineWidth',2);
			ylabel('Mass_{twc}(D) CDF [%]');


			if saveFigs
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				print([saveDir '/' probe '-MD-TempBinned/' flight '_' probe '_MD-Temp_S' num2str(ix) '_' num2str(iii) 'degC' outFileAppend],Ftype,Fres)

			end
		end
	end
end

if plotNtTemp
    for ix = loopVctr
		
		if plotAvg
			NtSprl = n_avg.(sprlNames{ix})';
			time_secs = time_secs_avg.(sprlNames{ix});
			time_fl = time_secsFL_avg.(sprlNames{ix});
			tempCsprl = tempC_avg.(sprlNames{ix});
		else
			NtSprl = n_orig.(sprlNames{ix})';
			time_secs = time_secs_orig.(sprlNames{ix});
			time_fl = time_secsFL_orig.(sprlNames{ix});
			tempCsprl = tempC_orig.(sprlNames{ix});
		end
		

        if saveFigs && noDisp
            figure('visible','off','Position', [10,10,1500,1000]);
        else
            figure('Position', [10,10,1500,1000]);
		end

				
		if ~plotAvg
			plot(NtSprl,tempCsprl,'b.','MarkerSize',10);
			title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - 1s Avg']);
		else
			plot(NtSprl,tempCsprl,'b-');
			title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - ' num2str(avgTime) 's Avg']);
		end
		
		
		% Plot ML top/bottom locations and annotate with temp
		hold on
		plot(NtRangeAll, [1 1]*mlTopTemp(ix),'k--')
		topStr = sprintf('%.3f %cC',mlTopTemp(ix),char(176));
		tMT = text(50,mlTopTemp(ix),topStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
		plot(NtRangeAll, [1 1]*mlBotTemp(ix),'k--')
		botStr = sprintf('%.3f %cC',mlBotTemp(ix),char(176));
		tMB = text(50,mlBotTemp(ix),botStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
		hold off


		ylabel(sprintf('Temp (%cC)', char(176)));
		xlabel('N_t [cm^{-3}]')
		if ~isempty(tempRangeAll)
			ylim(tempRangeAll);
		end
		if ~isempty(NtRangeAll)
			xlim(NtRangeAll);
		end
		set(gca,'XMinorTick','on','YMinorTick','on','YDir','reverse','Xscale','log');
		set(findall(gcf,'-property','FontSize'),'FontSize',28)
		set(tMB,'FontSize',14);
		set(tMT,'FontSize',14);
		grid;
		
        if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			if ~plotAvg
				print([saveDir '/' probe '-Nt-Temp_1s/' flight '_' probe '_Nt-Temp_1s_S' num2str(ix) outFileAppend],Ftype,Fres)
			else
				print([saveDir '/' probe '-Nt-Temp_' num2str(avgTime) 's/' flight '_' probe '_Nt-Temp_' num2str(avgTime) 's_S' num2str(ix) outFileAppend],Ftype,Fres)
			end
        end
    end
end

if plotNtTempAll && length(sprlNames) > 1
		
	colors = varycolor(length(sprlNames));
	
	if plotAvg
		Nt = n_avg.(sprlNames{1})';
		time_secs = time_secs_avg.(sprlNames{1});
		time_fl = time_secsFL_avg.(sprlNames{1});
		tempCsprl = tempC_avg.(sprlNames{1});
	else
		Nt = n_orig.(sprlNames{1})';
		time_secs = time_secs_orig.(sprlNames{1});
		time_fl = time_secsFL_orig.(sprlNames{1});
		tempCsprl = tempC_orig.(sprlNames{1});
	end
	
	if saveFigs && noDisp
		figure('visible','off','Position', [10,10,1500,1000]);
	else
		figure('Position', [10,10,1500,1000]);
	end

	
	if ~plotAvg
		plot(Nt,tempCsprl,'.','Color',colors(1,:),'MarkerSize',10);
		title([flight ' - All Spirals - ' probe ' - 1s Avg']);
	else
		plot(Nt,tempCsprl,'Color',colors(1,:));
		title([flight ' - All Spirals - ' probe ' - ' num2str(avgTime) 's Avg']);
	end


	ylabel(sprintf('Temp (%cC)', char(176)));
	xlabel('N_t [cm^{-3}]')
	if ~isempty(tempRangeAll)
		ylim(tempRangeAll);
	end
	if ~isempty(NtRangeAll)
		xlim(NtRangeAll);
	end
	set(gca,'XMinorTick','on','YMinorTick','on','YDir','reverse','Xscale','log');
	set(findall(gcf,'-property','FontSize'),'FontSize',28)
	hold on
	
	for ix = 2:length(sprlNames)
		
		if plotAvg
			Nt = n_avg.(sprlNames{ix})';
			time_secs = time_secs_avg.(sprlNames{ix});
			time_fl = time_secsFL_avg.(sprlNames{ix});
			tempCsprl = tempC_avg.(sprlNames{ix});
			plot(Nt,tempCsprl,'Color',colors(ix,:));
		else
			Nt = n_orig.(sprlNames{ix})';
			time_secs = time_secs_orig.(sprlNames{ix});
			time_fl = time_secsFL_orig.(sprlNames{ix});
			tempCsprl = tempC_orig.(sprlNames{ix});
			plot(Nt,tempCsprl,'.','Color',colors(ix,:),'MarkerSize',10);
		end

		
	end
	
	grid;
	
	% Create legends depending on which flight it is
	if strcmp(flight,'20150617')
		legend('S1','S2','S3','S4','S5','S6','S7','Location','eastoutside');
	elseif strcmp(flight,'20150620')
		legend('S1','S2','S3','S4','S5','S6','S7','Location','eastoutside');
	elseif strcmp(flight,'20150702')
		legend('S1','S2','S3','Location','eastoutside');
	elseif strcmp(flight,'20150706')
		legend('S1','S2','S3','S4','S5','S6','S7','S8','Location','eastoutside');
	elseif strcmp(flight,'20150709')
		legend('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10',...
			'S11','S12','S13','S14','S15','S16','Location','eastoutside');
	end
	
	if saveFigs
		set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		if ~plotAvg
			print([saveDir '/' probe '-Nt-Temp_1s/' flight '_' probe '_Nt-Temp_1s_All' outFileAppend],Ftype,Fres)
		else
			print([saveDir '/' probe '-Nt-Temp_' num2str(avgTime) 's/' flight '_' probe '_Nt-Temp_' num2str(avgTime) 's_All' outFileAppend],Ftype,Fres)
		end
	end
end

if plotTWCtemp
    for ix = loopVctr
		
		if plotAvg
			twc = twc_avg.(sprlNames{ix});
			time_secs = time_secs_avg.(sprlNames{ix});
			time_fl = time_secsFL_avg.(sprlNames{ix});
			tempCsprl = tempC_avg.(sprlNames{ix});
		else
			twc = twc_orig.(sprlNames{ix});
			time_secs = time_secs_orig.(sprlNames{ix});
			time_fl = time_secsFL_orig.(sprlNames{ix});
			tempCsprl = tempC_orig.(sprlNames{ix});
		end
		

        if saveFigs && noDisp
            figure('visible','off','Position', [10,10,1500,1000]);
        else
            figure('Position', [10,10,1500,1000]);
		end

		
		
		if ~plotAvg
			plot(twc,tempCsprl,'b.','MarkerSize',10);
			title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - 1s Avg']);
		else
			plot(twc,tempCsprl,'b-');
			title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - ' num2str(avgTime) 's Avg']);
		end
		
		% Plot ML top/bottom locations and annotate with temp
		hold on
		plot(TWCrangeAll, [1 1]*mlTopTemp(ix),'k--')
		topStr = sprintf('%.3f %cC',mlTopTemp(ix),char(176));
		tMT = text(7,mlTopTemp(ix),topStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
		plot(TWCrangeAll, [1 1]*mlBotTemp(ix),'k--')
		botStr = sprintf('%.3f %cC',mlBotTemp(ix),char(176));
		tMB = text(7,mlBotTemp(ix),botStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
		hold off


		ylabel(sprintf('Temp (%cC)', char(176)));
		xlabel('TWC [g m^{-3}]')
		if ~isempty(tempRangeAll)
			ylim(tempRangeAll);
		end
		if ~isempty(TWCrangeAll)
			xlim(TWCrangeAll);
		end
		set(gca,'XMinorTick','on','YMinorTick','on','YDir','reverse','Xscale','log');
		set(findall(gcf,'-property','FontSize'),'FontSize',28)
		set(tMB,'FontSize',14);
		set(tMT,'FontSize',14);
		grid;
		
        if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			if ~plotAvg
				print([saveDir '/' probe '-TWC-Temp_1s/' flight '_' probe '_TWC-Temp_1s_S' num2str(ix) outFileAppend],Ftype,Fres)
			else
				print([saveDir '/' probe '-TWC-Temp_' num2str(avgTime) 's/' flight '_' probe '_TWC-Temp_' num2str(avgTime) 's_S' num2str(ix) outFileAppend],Ftype,Fres)
			end
        end
    end
end

if plotTWCtempAll && length(sprlNames) > 1
	
	colors = varycolor(length(sprlNames));
	
	if plotAvg
		twc = twc_avg.(sprlNames{1});
		time_secs = time_secs_avg.(sprlNames{1});
		time_fl = time_secsFL_avg.(sprlNames{1});
		tempCsprl = tempC_avg.(sprlNames{1});
	else
		twc = twc_orig.(sprlNames{1});
		time_secs = time_secs_orig.(sprlNames{1});
		time_fl = time_secsFL_orig.(sprlNames{1});
		tempCsprl = tempC_orig.(sprlNames{1});
	end
	
	if saveFigs && noDisp
		figure('visible','off','Position', [10,10,1500,1000]);
	else
		figure('Position', [10,10,1500,1000]);
	end

	
	if ~plotAvg
		plot(twc,tempCsprl,'.','Color',colors(1,:),'MarkerSize',10);
		title([flight ' - All Spirals - ' probe ' - 1s Avg']);
	else
		plot(twc,tempCsprl,'Color',colors(1,:));
		title([flight ' - All Spirals - ' probe ' - ' num2str(avgTime) 's Avg']);
	end


	ylabel(sprintf('Temp (%cC)', char(176)));
	xlabel('TWC [g m^{-3}]')
	if ~isempty(tempRangeAll)
		ylim(tempRangeAll);
	end
	if ~isempty(TWCrangeAll)
		xlim(TWCrangeAll);
	end
	set(gca,'XMinorTick','on','YMinorTick','on','YDir','reverse','Xscale','log');
	set(findall(gcf,'-property','FontSize'),'FontSize',28)
	hold on
	
	for ix = 2:length(sprlNames)
		
		if plotAvg
			twc = twc_avg.(sprlNames{ix});
			time_secs = time_secs_avg.(sprlNames{ix});
			time_fl = time_secsFL_avg.(sprlNames{ix});
			tempCsprl = tempC_avg.(sprlNames{ix});
			plot(twc,tempCsprl,'Color',colors(ix,:));
		else
			twc = twc_orig.(sprlNames{ix});
			time_secs = time_secs_orig.(sprlNames{ix});
			time_fl = time_secsFL_orig.(sprlNames{ix});
			tempCsprl = tempC_orig.(sprlNames{ix});
			plot(twc,tempCsprl,'.','Color',colors(ix,:),'MarkerSize',10);
		end

		
	end
	
	grid;
	
	% Create legends depending on which flight it is
	if strcmp(flight,'20150617')
		legend('S1','S2','S3','S4','S5','S6','S7','Location','eastoutside');
	elseif strcmp(flight,'20150620')
		legend('S1','S2','S3','S4','S5','S6','S7','Location','eastoutside');
	elseif strcmp(flight,'20150702')
		legend('S1','S2','S3','Location','eastoutside');
	elseif strcmp(flight,'20150706')
		legend('S1','S2','S3','S4','S5','S6','S7','S8','Location','eastoutside');
	elseif strcmp(flight,'20150709')
		legend('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10',...
			'S11','S12','S13','S14','S15','S16','Location','eastoutside');
	end
	
	if saveFigs
		set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		if ~plotAvg
			print([saveDir '/' probe '-TWC-Temp_1s/' flight '_' probe '_TWC-Temp_1s_All' outFileAppend],Ftype,Fres)
		else
			print([saveDir '/' probe '-TWC-Temp_' num2str(avgTime) 's/' flight '_' probe '_TWC-Temp_' num2str(avgTime) 's_All' outFileAppend],Ftype,Fres)
		end
	end
end

if plotDmmTemp
    for ix = loopVctr
		
		if plotAvg
			Dmm = Dmm_twc_avg.(sprlNames{ix});
			time_secs = time_secs_avg.(sprlNames{ix});
			time_fl = time_secsFL_avg.(sprlNames{ix});
			tempCsprl = tempC_avg.(sprlNames{ix});
		else
			Dmm = Dmm_twc_orig.(sprlNames{ix});
			time_secs = time_secs_orig.(sprlNames{ix});
			time_fl = time_secsFL_orig.(sprlNames{ix});
			tempCsprl = tempC_orig.(sprlNames{ix});
		end
		
		
        if saveFigs && noDisp
            figure('visible','off','Position', [10,10,1500,1000]);
        else
            figure('Position', [10,10,1500,1000]);
		end

		plot(Dmm,tempCsprl,'b-');
		
		% Plot ML top/bottom locations and annotate with temp
		hold on
		plot(DmmRangeAll, [1 1]*mlTopTemp(ix),'k--')
		topStr = sprintf('%.3f %cC',mlTopTemp(ix),char(176));
		tMT = text(0.1,mlTopTemp(ix),topStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
		plot(DmmRangeAll, [1 1]*mlBotTemp(ix),'k--')
		botStr = sprintf('%.3f %cC',mlBotTemp(ix),char(176));
		tMB = text(0.1,mlBotTemp(ix),botStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
		hold off
		
		if ~plotAvg
			title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - 1s Avg']);
		else
			title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - ' num2str(avgTime) 's Avg']);
		end


		ylabel(sprintf('Temp (%cC)', char(176)));
		xlabel('D_{mm} [mm]')
		if ~isempty(tempRangeAll)
			ylim(tempRangeAll);
		end
		if ~isempty(DmmRangeAll)
			xlim(DmmRangeAll);
		end
		set(gca,'XMinorTick','on','YMinorTick','on','YDir','reverse');
		set(findall(gcf,'-property','FontSize'),'FontSize',28)
		set(tMB,'FontSize',14);
		set(tMT,'FontSize',14);
		grid;
		
        if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			if ~plotAvg
				print([saveDir '/' probe '-Dmm-Temp_1s/' flight '_' probe '_Dmm-Temp_1s_S' num2str(ix) outFileAppend],Ftype,Fres)
			else
				print([saveDir '/' probe '-Dmm-Temp_' num2str(avgTime) 's/' flight '_' probe '_Dmm-Temp_' num2str(avgTime) 's_S' num2str(ix) outFileAppend],Ftype,Fres)
			end
        end
    end
end

if plotDmmTempAll && length(sprlNames) > 1
	
	colors = varycolor(length(sprlNames));
	
	if plotAvg
		Dmm = Dmm_twc_avg.(sprlNames{1});
		time_secs = time_secs_avg.(sprlNames{1});
		time_fl = time_secsFL_avg.(sprlNames{1});
		tempCsprl = tempC_avg.(sprlNames{1});
	else
		Dmm = Dmm_twc_orig.(sprlNames{1});
		time_secs = time_secs_orig.(sprlNames{1});
		time_fl = time_secsFL_orig.(sprlNames{1});
		tempCsprl = tempC_orig.(sprlNames{1});
	end
	
	if saveFigs && noDisp
		figure('visible','off','Position', [10,10,1500,1000]);
	else
		figure('Position', [10,10,1500,1000]);
	end

	plot(Dmm,tempCsprl,'Color',colors(1,:));
	if ~plotAvg
		title([flight ' - All Spirals - ' probe ' - 1s Avg']);
	else
		title([flight ' - All Spirals - ' probe ' - ' num2str(avgTime) 's Avg']);
	end


	ylabel(sprintf('Temp (%cC)', char(176)));
	xlabel('D_{mm} [mm]')
	if ~isempty(tempRangeAll)
		ylim(tempRangeAll);
	end
	if ~isempty(DmmRangeAll)
		xlim(DmmRangeAll);
	end
	set(gca,'XMinorTick','on','YMinorTick','on','YDir','reverse');
	set(findall(gcf,'-property','FontSize'),'FontSize',28)
	hold on
	
	for ix = 2:length(sprlNames)
		
		if plotAvg
			Dmm = Dmm_twc_avg.(sprlNames{ix});
			time_secs = time_secs_avg.(sprlNames{ix});
			time_fl = time_secsFL_avg.(sprlNames{ix});
			tempCsprl = tempC_avg.(sprlNames{ix});
		else
			Dmm = Dmm_twc_orig.(sprlNames{ix});
			time_secs = time_secs_orig.(sprlNames{ix});
			time_fl = time_secsFL_orig.(sprlNames{ix});
			tempCsprl = tempC_orig.(sprlNames{ix});
		end

		plot(Dmm,tempCsprl,'Color',colors(ix,:));
	end
	
	% Create legends depending on which flight it is
	if strcmp(flight,'20150617')
		legend('S1','S2','S3','S4','S5','S6','S7','Location','eastoutside');
	elseif strcmp(flight,'20150620')
		legend('S1','S2','S3','S4','S5','S6','S7','Location','eastoutside');
	elseif strcmp(flight,'20150702')
		legend('S1','S2','S3','Location','eastoutside');
	elseif strcmp(flight,'20150706')
		legend('S1','S2','S3','S4','S5','S6','S7','S8','Location','eastoutside');
	elseif strcmp(flight,'20150709')
		legend('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10',...
			'S11','S12','S13','S14','S15','S16','Location','eastoutside');
	end
	
	grid;
	
	if saveFigs
		set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		if ~plotAvg
			print([saveDir '/' probe '-Dmm-Temp_1s/' flight '_' probe '_Dmm-Temp_1s_All' outFileAppend],Ftype,Fres)
		else
			print([saveDir '/' probe '-Dmm-Temp_' num2str(avgTime) 's/' flight '_' probe '_Dmm-Temp_' num2str(avgTime) 's_All' outFileAppend],Ftype,Fres)
		end
	end
end

if plotNtTempSprd
	
	if plotAvg
		NtSprl = n_avg.(sprlNames{1})';
		tempCsprl = tempC_avg.(sprlNames{1});
	else
		NtSprl = n_orig.(sprlNames{1})';
		tempCsprl = tempC_orig.(sprlNames{1});
	end
	
	allSprlNt = NtSprl;
	allSprlTemp = tempCsprl;

	% Determine the subsequent spiral profiles
	for ix = 2:length(startT)
		if plotAvg
			NtSprl = n_avg.(sprlNames{ix})';
			tempCsprl = tempC_avg.(sprlNames{ix});
		else
			NtSprl = n_orig.(sprlNames{ix})';
			tempCsprl = tempC_orig.(sprlNames{ix});
		end
		
		allSprlNt = vertcat(allSprlNt, NtSprl);
		allSprlTemp = vertcat(allSprlTemp, tempCsprl);
	end

% 	edgesTemp = (-15.125:.25:20.125);
% 	edgesTemp = (-19.125:.25:20.125);
	edgesTemp = (-19.25:.5:20.25);
% 	edgesTemp = (-19.5:1.0:20.5);
	numBins = length(edgesTemp)-1;

	bin_mid = (edgesTemp(1:end-1) + edgesTemp(2:end))/2;
	
	whichBinTemp = discretize(allSprlTemp,edgesTemp);

	for ix=1:numBins
		flagBinMemTemp = (whichBinTemp == ix);
		binNt = allSprlNt(flagBinMemTemp);
		if (isempty(binNt))
			binMean(ix) = NaN;
			binMedian(ix) = NaN;
			binMax(ix) = NaN;
			binMin(ix) = NaN;
			bin25pct(ix) = NaN;
			bin75pct(ix) = NaN;
		else
			binMean(ix) = nanmean(binNt);
			binMedian(ix) = nanmedian(binNt);
			binMax(ix) = max(binNt); % Max/min omit NaN by default
			binMin(ix) = min(binNt);
			bin25pct(ix) = prctile(binNt,25);
			bin75pct(ix) = prctile(binNt,75);
		end
	end
	
	if saveFigs && noDisp
		figure('visible','off','Position', [10,10,1500,1500]);
	else
		figure('Position', [10,10,1500,1500]);
	end
	
	if strcmp(cntrLine,'Mean')
		plot(bin_mid,binMean,'k--');
	elseif strcmp(cntrLine,'Median')
		plot(bin_mid,binMedian,'k--');
	end
	hold on
	plot(bin_mid,binMin,'k');
	plot(bin_mid,binMax,'k');
	plot(bin_mid,bin25pct,'Color',[0.82 0.82 0.82]);
	plot(bin_mid,bin75pct,'Color',[0.82 0.82 0.82]);
	set(gca,'XMinorTick','on','YMinorTick','on','XDir','reverse','Yscale','log');
	if ~isempty(tempRangeAll)
		xlim(tempRangeAll);
	end
	if ~isempty(NtRangeAll)
		ylim(NtRangeAll)
	end
	view([90 -90])
	grid
	if ~plotAvg
		title([flight ' - Nt vs. Temp - ' cntrLine ' and Spread all Spirals - ' probe ' - 1s Avg']);
	else
		title([flight ' - Nt vs. Temp - ' cntrLine ' and Spread all Spirals - ' probe ' - ' num2str(avgTime) 's Avg']);
	end
	ylabel('N_t [cm^{-3}]')
	xlabel(sprintf('Temp (%cC)', char(176)));
	set(findall(gcf,'-property','FontSize'),'FontSize',28);

	if saveFigs
		set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		if ~plotAvg
			print([saveDir '/' probe '-Nt-Temp_1s/' flight '_' probe '_Nt-Temp_' cntrLine '-Spread_1s' outFileAppend],Ftype,Fres)
		else
			print([saveDir '/' probe '-Nt-Temp_' num2str(avgTime) 's/' flight '_' probe '_Nt-Temp_' cntrLine '-Spread_' num2str(avgTime) 's' outFileAppend],Ftype,Fres)
		end
	end
end

if plotTWCtempSprd
	
	if plotAvg
		TWCsprl = twc_avg.(sprlNames{1});
		tempCsprl = tempC_avg.(sprlNames{1});
	else
		TWCsprl = twc_orig.(sprlNames{1});
		tempCsprl = tempC_orig.(sprlNames{1});
	end
	
	allSprlTWC = TWCsprl;
	allSprlTemp = tempCsprl;

	% Determine the subsequent spiral profiles
	for ix = 2:length(startT)
		if plotAvg
			TWCsprl = twc_avg.(sprlNames{ix});
			tempCsprl = tempC_avg.(sprlNames{ix});
		else
			TWCsprl = twc_orig.(sprlNames{ix});
			tempCsprl = tempC_orig.(sprlNames{ix});
		end
		
		allSprlTWC = vertcat(allSprlTWC, TWCsprl);
		allSprlTemp = vertcat(allSprlTemp, tempCsprl);
	end

	edgesTemp = (-15.125:.25:20.125);
	numBins = length(edgesTemp)-1;

	bin_mid = (edgesTemp(1:end-1) + edgesTemp(2:end))/2;
	
	whichBinTemp = discretize(allSprlTemp,edgesTemp);

	for ix=1:numBins
		flagBinMemTemp = (whichBinTemp == ix);
		binTWC = allSprlTWC(flagBinMemTemp);
		if (isempty(binTWC))
			binMean(ix) = NaN;
			binMedian(ix) = NaN;
			binMax(ix) = NaN;
			binMin(ix) = NaN;
			bin25pct(ix) = NaN;
			bin75pct(ix) = NaN;
		else
			binMean(ix) = nanmean(binTWC);
			binMedian(ix) = nanmedian(binTWC);
			binMax(ix) = max(binTWC);
			binMin(ix) = min(binTWC);
			bin25pct(ix) = prctile(binTWC,25);
			bin75pct(ix) = prctile(binTWC,75);
		end
	end
	
	if saveFigs && noDisp
		figure('visible','off','Position', [10,10,1500,1500]);
	else
		figure('Position', [10,10,1500,1500]);
	end
	
	if strcmp(cntrLine,'Mean')
		plot(bin_mid,binMean,'k--');
	elseif strcmp(cntrLine,'Median')
		plot(bin_mid,binMedian,'k--');
	end
	hold on
	plot(bin_mid,binMin,'k');
	plot(bin_mid,binMax,'k');
	plot(bin_mid,bin25pct,'Color',[0.82 0.82 0.82]);
	plot(bin_mid,bin75pct,'Color',[0.82 0.82 0.82])
	set(gca,'XMinorTick','on','YMinorTick','on','XDir','reverse','Yscale','log');
	if ~isempty(tempRangeAll)
		xlim(tempRangeAll);
	end
	if ~isempty(TWCrangeAll)
%		ylim(TWCrangeAll)
	end
	view([90 -90])
	grid
	if ~plotAvg
		title([flight ' - TWC vs. Temp - ' cntrLine ' and Spread all Spirals - ' probe ' - 1s Avg']);
	else
		title([flight ' - TWC vs. Temp - ' cntrLine ' and Spread all Spirals - ' probe ' - ' num2str(avgTime) 's Avg']);
	end
	ylabel('TWC [g m^{-3}]');
	xlabel(sprintf('Temp (%cC)', char(176)));
	set(findall(gcf,'-property','FontSize'),'FontSize',28);

	if saveFigs
		set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		if ~plotAvg
			print([saveDir '/' probe '-TWC-Temp_1s/' flight '_' probe '_TWC-Temp_' cntrLine '-Spread_1s' outFileAppend],Ftype,Fres)
		else
			print([saveDir '/' probe '-TWC-Temp_' num2str(avgTime) 's/' flight '_' probe '_TWC-Temp_' cntrLine '-Spread_' num2str(avgTime) 's' outFileAppend],Ftype,Fres)
		end
	end
end

if plotDmmTempSprd
	
	if plotAvg
		DmmSprl = Dmm_twc_avg.(sprlNames{1});
		tempCsprl = tempC_avg.(sprlNames{1});
	else
		DmmSprl = Dmm_twc_orig.(sprlNames{1});
		tempCsprl = tempC_orig.(sprlNames{1});
	end
	
	allSprlDmm = DmmSprl;
	allSprlTemp = tempCsprl;

	% Determine the subsequent spiral profiles
	for ix = 2:length(startT)
		if plotAvg
			DmmSprl = Dmm_twc_avg.(sprlNames{ix});
			tempCsprl = tempC_avg.(sprlNames{ix});
		else
			DmmSprl = Dmm_twc_orig.(sprlNames{ix});
			tempCsprl = tempC_orig.(sprlNames{ix});
		end
		
		allSprlDmm = vertcat(allSprlDmm, DmmSprl);
		allSprlTemp = vertcat(allSprlTemp, tempCsprl);
	end

	edgesTemp = (-15.125:.25:20.125);
	numBins = length(edgesTemp)-1;

% 	bin_min = edgesTemp(1:end-1); 
% 	bin_max = edgesTemp(2:end); 
	bin_mid = (edgesTemp(1:end-1) + edgesTemp(2:end))/2;
	
	whichBinTemp = discretize(allSprlTemp,edgesTemp);

	for ix=1:numBins
		flagBinMemTemp = (whichBinTemp == ix);
		binDmm = allSprlDmm(flagBinMemTemp);
		if (isempty(binDmm))
			binMean(ix) = NaN;
			binMedian(ix) = NaN;
			binMax(ix) = NaN;
			binMin(ix) = NaN;
			bin25pct(ix) = NaN;
			bin75pct(ix) = NaN;
		else
			binMean(ix) = nanmean(binDmm);
			binMedian(ix) = nanmedian(binDmm);
			binMax(ix) = max(binDmm);
			binMin(ix) = min(binDmm);
			bin25pct(ix) = prctile(binDmm,25);
			bin75pct(ix) = prctile(binDmm,75);
		end
	end
	
	if saveFigs && noDisp
		figure('visible','off','Position', [10,10,1500,1500]);
	else
		figure('Position', [10,10,1500,1500]);
	end
	
	if strcmp(cntrLine,'Mean')
		plot(bin_mid,binMean,'k--');
	elseif strcmp(cntrLine,'Median')
		plot(bin_mid,binMedian,'k--');
	end
	hold on
	plot(bin_mid,binMin,'k');
	plot(bin_mid,binMax,'k');
	plot(bin_mid,bin25pct,'Color',[0.82 0.82 0.82]);
	plot(bin_mid,bin75pct,'Color',[0.82 0.82 0.82])
	set(gca,'XMinorTick','on','YMinorTick','on','XDir','reverse');
	if ~isempty(tempRangeAll)
		xlim(tempRangeAll);
	end
	if ~isempty(DmmRangeAll)
		ylim(DmmRangeAll)
	end
	view([90 -90])
	grid
	if ~plotAvg
		title([flight ' - Dmm vs. Temp - ' cntrLine ' and Spread all Spirals - ' probe ' - 1s Avg']);
	else
		title([flight ' - Dmm vs. Temp - ' cntrLine ' and Spread all Spirals - ' probe ' - ' num2str(avgTime) 's Avg']);
	end
	ylabel('D_{mm} [mm]');
	xlabel(sprintf('Temp (%cC)', char(176)));
	set(findall(gcf,'-property','FontSize'),'FontSize',28);

	if saveFigs
		set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		if ~plotAvg
			print([saveDir '/' probe '-Dmm-Temp_1s/' flight '_' probe '_Dmm-Temp_' cntrLine '-Spread_1s' outFileAppend],Ftype,Fres)
		else
			print([saveDir '/' probe '-Dmm-Temp_' num2str(avgTime) 's/' flight '_' probe '_Dmm-Temp_' cntrLine '-Spread_' num2str(avgTime) 's' outFileAppend],Ftype,Fres)
		end
	end
end

if plotNtTempSprdF
	
	if ~isempty(tzSprls)
	
		if plotAvg
			NtSprl = n_avg.(sprlNames{tzSprls(1)})';
			tempCsprl = tempC_avg.(sprlNames{tzSprls(1)});
		else
			NtSprl = n_orig.(sprlNames{tzSprls(1)})';
			tempCsprl = tempC_orig.(sprlNames{tzSprls(1)});
		end

		tzNt = NtSprl;
		tzTemp = tempCsprl;

		% Determine the subsequent spiral profiles
		if length(tzSprls) > 1
			for ix = 2:length(tzSprls)
				if plotAvg
					NtSprl = n_avg.(sprlNames{tzSprls(ix)})';
					tempCsprl = tempC_avg.(sprlNames{tzSprls(ix)});
				else
					NtSprl = n_orig.(sprlNames{tzSprls(ix)})';
					tempCsprl = tempC_orig.(sprlNames{tzSprls(ix)});
				end

				tzNt = vertcat(tzNt, NtSprl);
				tzTemp = vertcat(tzTemp, tempCsprl);
			end
		end
		
		tzNt_orig = tzNt;
		tzNt(tzNt == 0) = 1e-10; % Set zeros to a very small value to allow plotting/fill to work properly
		tzNt(isnan(tzNt)) = 1e-10;
	end
	
	
	if ~isempty(esrSprls)
	
		if plotAvg
			NtSprl = n_avg.(sprlNames{esrSprls(1)})';
			tempCsprl = tempC_avg.(sprlNames{esrSprls(1)});
		else
			NtSprl = n_orig.(sprlNames{esrSprls(1)})';
			tempCsprl = tempC_orig.(sprlNames{esrSprls(1)});
		end

		esrNt = NtSprl;
		esrTemp = tempCsprl;

		% Determine the subsequent spiral profiles
		if length(esrSprls) > 1
			for ix = 2:length(esrSprls)
				if plotAvg
					NtSprl = n_avg.(sprlNames{esrSprls(ix)})';
					tempCsprl = tempC_avg.(sprlNames{esrSprls(ix)});
				else
					NtSprl = n_orig.(sprlNames{esrSprls(ix)})';
					tempCsprl = tempC_orig.(sprlNames{esrSprls(ix)});
				end

				esrNt = vertcat(esrNt, NtSprl);
				esrTemp = vertcat(esrTemp, tempCsprl);
			end
		end
		
		esrNt_orig = esrNt;
		esrNt(esrNt == 0) = 1e-10; % Set zeros to a very small value to allow plotting/fill to work properly
		esrNt(isnan(esrNt)) = 1e-10;
	end
	
	
	if ~isempty(raSprls)
	
		if plotAvg
			NtSprl = n_avg.(sprlNames{raSprls(1)})';
			tempCsprl = tempC_avg.(sprlNames{raSprls(1)});
		else
			NtSprl = n_orig.(sprlNames{raSprls(1)})';
			tempCsprl = tempC_orig.(sprlNames{raSprls(1)});
		end

		raNt = NtSprl;
		raTemp = tempCsprl;

		% Determine the subsequent spiral profiles
		if length(raSprls) > 1
			for ix = 2:length(raSprls)
				if plotAvg
					NtSprl = n_avg.(sprlNames{raSprls(ix)})';
					tempCsprl = tempC_avg.(sprlNames{raSprls(ix)});
				else
					NtSprl = n_orig.(sprlNames{raSprls(ix)})';
					tempCsprl = tempC_orig.(sprlNames{raSprls(ix)});
				end

				raNt = vertcat(raNt, NtSprl);
				raTemp = vertcat(raTemp, tempCsprl);
			end
		end
	end

	edgesTemp = (-15.125:.25:20.125);
	numBins = length(edgesTemp)-1;

	
	if saveFigs && noDisp
		figure('visible','off','Position', [10,10,1500,1500]);
	else
		figure('Position', [10,10,1500,1500]);
	end
	hold on
	
	
	if ~isempty(tzSprls)
		tzNtBin_mid = (edgesTemp(1:end-1) + edgesTemp(2:end))/2;

		tzNtTempBinIx = discretize(tzTemp,edgesTemp);

		tzNtFlgRmvBin = [];
		tzNtBinMean = [];
		tzNtBinMedian = [];
		tzNtBinMax = [];
		tzNtBinMin = [];
		tzNtBin25p = [];
		tzNtBin75p = [];

		for ix=1:numBins
			tzNtBinMem = (tzNtTempBinIx == ix);
			tzBinNt = tzNt(tzNtBinMem);
			tzBinNt_orig = tzNt_orig(tzNtBinMem);
			
			if ~isempty(tzBinNt)
				tzNtBinMean = horzcat(tzNtBinMean,nanmean(tzBinNt_orig));
				tzNtBinMedian = horzcat(tzNtBinMedian,nanmedian(tzBinNt));
				tzNtBinMax = horzcat(tzNtBinMax,max(tzBinNt));
				tzNtBinMin = horzcat(tzNtBinMin,min(tzBinNt));
				tzNtBin25p = horzcat(tzNtBin25p,prctile(tzBinNt,25));
				tzNtBin75p = horzcat(tzNtBin75p,prctile(tzBinNt,75));
				
			else
				tzNtFlgRmvBin = horzcat(tzNtFlgRmvBin,ix);
			end
		end

		tzNtBin_mid(tzNtFlgRmvBin) = [];

		tzNtBin_midFlip = [tzNtBin_mid,fliplr(tzNtBin_mid)];
% 		tzNtSpread = [log10(tzNtBinMin), fliplr(log10(tzNtBinMax))];
		tzNtSpread = [log10(tzNtBin25p), fliplr(log10(tzNtBin75p))];

		f1 = fill(tzNtBin_midFlip,tzNtSpread,'b','FaceAlpha', 0.5);

		if strcmp(cntrLine,'Mean')
			m1 = plot(tzNtBin_mid,log10(tzNtBinMean),'Color',[27/255 27/255 100/255],'LineWidth', 3);
		elseif strcmp(cntrLine,'Median')
			m1 = plot(tzNtBin_mid,log10(tzNtBinMedian),'Color',[27/255 27/255 100/255],'LineWidth', 3);
		end
		
	end
	
	if ~isempty(esrSprls)
		esrNtBin_mid = (edgesTemp(1:end-1) + edgesTemp(2:end))/2;

		esrNtTempBinIx = discretize(esrTemp,edgesTemp);

		esrNtFlgRmvBin = [];
		esrNtBinMean = [];
		esrNtBinMedian = [];
		esrNtBinMax = [];
		esrNtBinMin = [];
		esrNtBin25p = [];
		esrNtBin75p = [];
		
		for ix=1:numBins
			esrNtBinMem = (esrNtTempBinIx == ix);
			esrBinNt = esrNt(esrNtBinMem);

			if ~isempty(esrBinNt)
				esrNtBinMean = horzcat(esrNtBinMean,nanmean(esrBinNt));
				esrNtBinMedian = horzcat(esrNtBinMedian,nanmedian(esrBinNt));
				esrNtBinMax = horzcat(esrNtBinMax,max(esrBinNt));
				esrNtBinMin = horzcat(esrNtBinMin,min(esrBinNt));
				esrNtBin25p = horzcat(esrNtBin25p,prctile(esrBinNt,25));
				esrNtBin75p = horzcat(esrNtBin75p,prctile(esrBinNt,75));
				
			else
				esrNtFlgRmvBin = horzcat(esrNtFlgRmvBin,ix);
			end
		end

		esrNtBin_mid(esrNtFlgRmvBin) = [];

		esrNtBin_midFlip = [esrNtBin_mid,fliplr(esrNtBin_mid)];
% 		esrNtSpread = [log10(esrNtBinMin), fliplr(log10(esrNtBinMax))];
		esrNtSpread = [log10(esrNtBin25p), fliplr(log10(esrNtBin75p))];

		f2 = fill(esrNtBin_midFlip,esrNtSpread,'r','FaceAlpha', 0.5);

		if strcmp(cntrLine,'Mean')
			m2 = plot(esrNtBin_mid,log10(esrNtBinMean),'Color',[100/255 27/255 27/255],'LineWidth', 3);
		elseif strcmp(cntrLine,'Median')
			m2 = plot(esrNtBin_mid,log10(esrNtBinMedian),'Color','y','LineWidth', 3);
		end
	end
	
	
	if ~plotAvg
		title([flight ' - Nt vs. Temp - ' cntrLine ' and Spread - ' probe ' - 1s Avg']);
	else
		title([flight ' - Nt vs. Temp - ' cntrLine ' and Spread - ' probe ' - ' num2str(avgTime) 's Avg']);
	end
	set(gca,'XMinorTick','on','YMinorTick','on','XDir','reverse');
	if ~isempty(tempRangeAll)
		xlim(tempRangeAll);
	end
	if ~isempty(NtRangeFill)
		ylim(NtRangeFill)
	end
	ylabel('log_{10}(N_t) [cm^{-3}]')
	xlabel(sprintf('Temp (%cC)', char(176)));
	set(findall(gcf,'-property','FontSize'),'FontSize',28);
	view([90 -90])
	
	legend([f1 f2],{'Transition Zone','Enhanced Stratiform'},'Location','eastoutside');
	

	if saveFigs
		set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		if ~plotAvg
			print([saveDir '/' probe '-Nt-Temp_1s/' flight '_' probe '_Nt-Temp_' cntrLine '-Spread-Fill_1s' outFileAppend],Ftype,Fres)
		else
			print([saveDir '/' probe '-Nt-Temp_' num2str(avgTime) 's/' flight '_' probe '_Nt-Temp_' cntrLine '-Spread-Fill_' num2str(avgTime) 's' outFileAppend],Ftype,Fres)
		end
	end
end

if plotTWCTempSprdF
	
	if ~isempty(tzSprls)
	
		if plotAvg
			TWCSprl = twc_avg.(sprlNames{tzSprls(1)});
			tempCsprl = tempC_avg.(sprlNames{tzSprls(1)});
		else
			TWCSprl = twc_orig.(sprlNames{tzSprls(1)});
			tempCsprl = tempC_orig.(sprlNames{tzSprls(1)});
		end

		tzTWC = TWCSprl;
		tzTemp = tempCsprl;

		% Determine the subsequent spiral profiles
		if length(tzSprls) > 1
			for ix = 2:length(tzSprls)
				if plotAvg
					TWCSprl = twc_avg.(sprlNames{tzSprls(ix)});
					tempCsprl = tempC_avg.(sprlNames{tzSprls(ix)});
				else
					TWCSprl = twc_orig.(sprlNames{tzSprls(ix)});
					tempCsprl = tempC_orig.(sprlNames{tzSprls(ix)});
				end

				tzTWC = vertcat(tzTWC, TWCSprl);
				tzTemp = vertcat(tzTemp, tempCsprl);
			end
		end
		
		tzTWC_orig = tzTWC;
		tzTWC(tzTWC == 0) = 1e-8; % Set zeros to a very small value to allow plotting/fill to work properly
		tzTWC(isnan(tzTWC)) = 1e-8;
	end
	
	
	if ~isempty(esrSprls)
	
		if plotAvg
			TWCSprl = twc_avg.(sprlNames{esrSprls(1)});
			tempCsprl = tempC_avg.(sprlNames{esrSprls(1)});
		else
			TWCSprl = twc_orig.(sprlNames{esrSprls(1)});
			tempCsprl = tempC_orig.(sprlNames{esrSprls(1)});
		end

		esrTWC = TWCSprl;
		esrTemp = tempCsprl;

		% Determine the subsequent spiral profiles
		if length(esrSprls) > 1
			for ix = 2:length(esrSprls)
				if plotAvg
					TWCSprl = twc_avg.(sprlNames{esrSprls(ix)});
					tempCsprl = tempC_avg.(sprlNames{esrSprls(ix)});
				else
					TWCSprl = twc_orig.(sprlNames{esrSprls(ix)});
					tempCsprl = tempC_orig.(sprlNames{esrSprls(ix)});
				end

				esrTWC = vertcat(esrTWC, TWCSprl);
				esrTemp = vertcat(esrTemp, tempCsprl);
			end
		end
		
		esrTWC_orig = esrTWC;
		esrTWC(esrTWC == 0) = 1e-8; % Set zeros to a very small value to allow plotting/fill to work properly
		esrTWC(isnan(esrTWC)) = 1e-8;
	end
	
	
	if ~isempty(raSprls)
	
		if plotAvg
			TWCSprl = twc_avg.(sprlNames{raSprls(1)})';
			tempCsprl = tempC_avg.(sprlNames{raSprls(1)});
		else
			TWCSprl = twc_orig.(sprlNames{raSprls(1)})';
			tempCsprl = tempC_orig.(sprlNames{raSprls(1)});
		end

		raTWC = TWCSprl;
		raTemp = tempCsprl;

		% Determine the subsequent spiral profiles
		if length(raSprls) > 1
			for ix = 2:length(raSprls)
				if plotAvg
					TWCSprl = twc_avg.(sprlNames{raSprls(ix)})';
					tempCsprl = tempC_avg.(sprlNames{raSprls(ix)});
				else
					TWCSprl = twc_orig.(sprlNames{raSprls(ix)})';
					tempCsprl = tempC_orig.(sprlNames{raSprls(ix)});
				end

				raTWC = vertcat(raTWC, TWCSprl);
				raTemp = vertcat(raTemp, tempCsprl);
			end
		end
	end

	edgesTemp = (-15.125:.25:20.125);
	numBins = length(edgesTemp)-1;

	
	if saveFigs && noDisp
		figure('visible','off','Position', [10,10,1500,1500]);
	else
		figure('Position', [10,10,1500,1500]);
	end
	hold on
	
	
	if ~isempty(tzSprls)
		tzTWCBin_mid = (edgesTemp(1:end-1) + edgesTemp(2:end))/2;

		tzTWCTempBinIx = discretize(tzTemp,edgesTemp);

		tzTWCFlgRmvBin = [];
		tzTWCBinMean = [];
		tzTWCBinMedian = [];
		tzTWCBinMax = [];
		tzTWCBinMin = [];
		tzTWCBin25p = [];
		tzTWCBin75p = [];

		for ix=1:numBins
			tzTWCBinMem = (tzTWCTempBinIx == ix);
			tzBinTWC = tzTWC(tzTWCBinMem);
			tzBinTWC_orig = tzTWC_orig(tzTWCBinMem);
			
			if ~isempty(tzBinTWC)
				tzTWCBinMean = horzcat(tzTWCBinMean,nanmean(tzBinTWC_orig));
				tzTWCBinMedian = horzcat(tzTWCBinMedian,nanmedian(tzBinTWC));
				tzTWCBinMax = horzcat(tzTWCBinMax,max(tzBinTWC));
				tzTWCBinMin = horzcat(tzTWCBinMin,min(tzBinTWC));
				tzTWCBin25p = horzcat(tzTWCBin25p,prctile(tzBinTWC,25));
				tzTWCBin75p = horzcat(tzTWCBin75p,prctile(tzBinTWC,75));
				
			else
				tzTWCFlgRmvBin = horzcat(tzTWCFlgRmvBin,ix);
			end
		end

		tzTWCBin_mid(tzTWCFlgRmvBin) = [];

		tzTWCBin_midFlip = [tzTWCBin_mid,fliplr(tzTWCBin_mid)];
% 		tzTWCSpread = [log10(tzTWCBinMin), fliplr(log10(tzTWCBinMax))];
		tzTWCSpread = [log10(tzTWCBin25p), fliplr(log10(tzTWCBin75p))];

		f1 = fill(tzTWCBin_midFlip,tzTWCSpread,'b','FaceAlpha', 0.5);

		if strcmp(cntrLine,'Mean')
			m1 = plot(tzTWCBin_mid,log10(tzTWCBinMean),'Color',[27/255 27/255 100/255],'LineWidth', 3);
		elseif strcmp(cntrLine,'Median')
			m1 = plot(tzTWCBin_mid,log10(tzTWCBinMedian),'Color',[27/255 27/255 100/255],'LineWidth', 3);
		end	
	end
	
	if ~isempty(esrSprls)
		esrTWCBin_mid = (edgesTemp(1:end-1) + edgesTemp(2:end))/2;

		esrTWCTempBinIx = discretize(esrTemp,edgesTemp);

		esrTWCFlgRmvBin = [];
		esrTWCBinMean = [];
		esrTWCBinMedian = [];
		esrTWCBinMax = [];
		esrTWCBinMin = [];
		esrTWCBin25p = [];
		esrTWCBin75p = [];
		
		for ix=1:numBins
			esrTWCBinMem = (esrTWCTempBinIx == ix);
			esrBinTWC = esrTWC(esrTWCBinMem);

			if ~isempty(esrBinTWC)
				esrTWCBinMean = horzcat(esrTWCBinMean,nanmean(esrBinTWC));
				esrTWCBinMedian = horzcat(esrTWCBinMedian,nanmedian(esrBinTWC));
				esrTWCBinMax = horzcat(esrTWCBinMax,max(esrBinTWC));
				esrTWCBinMin = horzcat(esrTWCBinMin,min(esrBinTWC));
				esrTWCBin25p = horzcat(esrTWCBin25p,prctile(esrBinTWC,25));
				esrTWCBin75p = horzcat(esrTWCBin75p,prctile(esrBinTWC,75));
				
			else
				esrTWCFlgRmvBin = horzcat(esrTWCFlgRmvBin,ix);
			end
		end

		esrTWCBin_mid(esrTWCFlgRmvBin) = [];

		esrTWCBin_midFlip = [esrTWCBin_mid,fliplr(esrTWCBin_mid)];
% 		esrTWCSpread = [log10(esrTWCBinMin), fliplr(log10(esrTWCBinMax))];
		esrTWCSpread = [log10(esrTWCBin25p), fliplr(log10(esrTWCBin75p))];

		f2 = fill(esrTWCBin_midFlip,esrTWCSpread,'r','FaceAlpha', 0.5);

		if strcmp(cntrLine,'Mean')
			m2 = plot(esrTWCBin_mid,log10(esrTWCBinMean),'Color',[100/255 27/255 27/255],'LineWidth', 3);
		elseif strcmp(cntrLine,'Median')
			m2 = plot(esrTWCBin_mid,log10(esrTWCBinMedian),'Color','y','LineWidth', 3);
		end	
	end
	
	
	if ~plotAvg
		title([flight ' - TWC vs. Temp - ' cntrLine ' and Spread - ' probe ' - 1s Avg']);
	else
		title([flight ' - TWC vs. Temp - ' cntrLine ' and Spread - ' probe ' - ' num2str(avgTime) 's Avg']);
	end
	set(gca,'XMinorTick','on','YMinorTick','on','XDir','reverse');
	if ~isempty(tempRangeAll)
		xlim(tempRangeAll);
	end
	if ~isempty(TWCRangeFill)
		ylim(TWCRangeFill)
	end
	ylabel('log_{10}(TWC) [g m^{-3}]')
	xlabel(sprintf('Temp (%cC)', char(176)));
	set(findall(gcf,'-property','FontSize'),'FontSize',28);
	view([90 -90])
	
	legend([f1 f2],{'Transition Zone','Enhanced Stratiform'},'Location','eastoutside');
	

	if saveFigs
		set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		if ~plotAvg
			print([saveDir '/' probe '-TWC-Temp_1s/' flight '_' probe '_TWC-Temp_' cntrLine '-Spread-Fill_1s' outFileAppend],Ftype,Fres)
		else
			print([saveDir '/' probe '-TWC-Temp_' num2str(avgTime) 's/' flight '_' probe '_TWC-Temp_' cntrLine '-Spread-Fill_' num2str(avgTime) 's' outFileAppend],Ftype,Fres)
		end
	end
end

if plotDmmTempSprdF
	
	if ~isempty(tzSprls)
	
		if plotAvg
			DmmSprl = Dmm_twc_avg.(sprlNames{tzSprls(1)});
			tempCsprl = tempC_avg.(sprlNames{tzSprls(1)});
		else
			DmmSprl = Dmm_twc_orig.(sprlNames{tzSprls(1)});
			tempCsprl = tempC_orig.(sprlNames{tzSprls(1)});
		end

		tzDmm = DmmSprl;
		tzTemp = tempCsprl;

		% Determine the subsequent spiral profiles
		if length(tzSprls) > 1
			for ix = 2:length(tzSprls)
				if plotAvg
					DmmSprl = Dmm_twc_avg.(sprlNames{tzSprls(ix)});
					tempCsprl = tempC_avg.(sprlNames{tzSprls(ix)});
				else
					DmmSprl = Dmm_twc_orig.(sprlNames{tzSprls(ix)});
					tempCsprl = tempC_orig.(sprlNames{tzSprls(ix)});
				end

				tzDmm = vertcat(tzDmm, DmmSprl);
				tzTemp = vertcat(tzTemp, tempCsprl);
			end
		end
		
		tzDmm_orig = tzDmm;
		tzDmm(tzDmm == 0) = 1e-8; % Set zeros to a very small value to allow plotting/fill to work properly
		tzDmm(isnan(tzDmm)) = 1e-8;
	end
	
	
	if ~isempty(esrSprls)
	
		if plotAvg
			DmmSprl = Dmm_twc_avg.(sprlNames{esrSprls(1)});
			tempCsprl = tempC_avg.(sprlNames{esrSprls(1)});
		else
			DmmSprl = Dmm_twc_orig.(sprlNames{esrSprls(1)});
			tempCsprl = tempC_orig.(sprlNames{esrSprls(1)});
		end

		esrDmm = DmmSprl;
		esrTemp = tempCsprl;

		% Determine the subsequent spiral profiles
		if length(esrSprls) > 1
			for ix = 2:length(esrSprls)
				if plotAvg
					DmmSprl = Dmm_twc_avg.(sprlNames{esrSprls(ix)});
					tempCsprl = tempC_avg.(sprlNames{esrSprls(ix)});
				else
					DmmSprl = Dmm_twc_orig.(sprlNames{esrSprls(ix)});
					tempCsprl = tempC_orig.(sprlNames{esrSprls(ix)});
				end

				esrDmm = vertcat(esrDmm, DmmSprl);
				esrTemp = vertcat(esrTemp, tempCsprl);
			end
		end
		
		esrDmm_orig = esrDmm;
		esrDmm(esrDmm == 0) = 1e-8; % Set zeros to a very small value to allow plotting/fill to work properly
		esrDmm(isnan(esrDmm)) = 1e-8;
	end
	
	
	if ~isempty(raSprls)
	
		if plotAvg
			DmmSprl = Dmm_twc_avg.(sprlNames{raSprls(1)})';
			tempCsprl = tempC_avg.(sprlNames{raSprls(1)});
		else
			DmmSprl = Dmm_twc_orig.(sprlNames{raSprls(1)})';
			tempCsprl = tempC_orig.(sprlNames{raSprls(1)});
		end

		raDmm = DmmSprl;
		raTemp = tempCsprl;

		% Determine the subsequent spiral profiles
		if length(raSprls) > 1
			for ix = 2:length(raSprls)
				if plotAvg
					DmmSprl = Dmm_twc_avg.(sprlNames{raSprls(ix)})';
					tempCsprl = tempC_avg.(sprlNames{raSprls(ix)});
				else
					DmmSprl = Dmm_twc_orig.(sprlNames{raSprls(ix)})';
					tempCsprl = tempC_orig.(sprlNames{raSprls(ix)});
				end

				raDmm = vertcat(raDmm, DmmSprl);
				raTemp = vertcat(raTemp, tempCsprl);
			end
		end
	end

	edgesTemp = (-15.125:.25:20.125);
	numBins = length(edgesTemp)-1;

	
	if saveFigs && noDisp
		figure('visible','off','Position', [10,10,1500,1500]);
	else
		figure('Position', [10,10,1500,1500]);
	end
	hold on
	
	
	if ~isempty(tzSprls)
		tzDmmBin_mid = (edgesTemp(1:end-1) + edgesTemp(2:end))/2;

		tzDmmTempBinIx = discretize(tzTemp,edgesTemp);

		tzDmmFlgRmvBin = [];
		tzDmmBinMean = [];
		tzDmmBinMedian = [];
		tzDmmBinMax = [];
		tzDmmBinMin = [];
		tzDmmBin25p = [];
		tzDmmBin75p = [];

		for ix=1:numBins
			tzDmmBinMem = (tzDmmTempBinIx == ix);
			tzBinDmm = tzDmm(tzDmmBinMem);
			tzBinDmm_orig = tzDmm_orig(tzDmmBinMem);
			
			if ~isempty(tzBinDmm)
				tzDmmBinMean = horzcat(tzDmmBinMean,nanmean(tzBinDmm_orig));
				tzDmmBinMedian = horzcat(tzDmmBinMedian,nanmedian(tzBinDmm));
				tzDmmBinMax = horzcat(tzDmmBinMax,max(tzBinDmm));
				tzDmmBinMin = horzcat(tzDmmBinMin,min(tzBinDmm));
				tzDmmBin25p = horzcat(tzDmmBin25p,prctile(tzBinDmm,25));
				tzDmmBin75p = horzcat(tzDmmBin75p,prctile(tzBinDmm,75));
				
			else
				tzDmmFlgRmvBin = horzcat(tzDmmFlgRmvBin,ix);
			end
		end

		tzDmmBin_mid(tzDmmFlgRmvBin) = [];

		tzDmmBin_midFlip = [tzDmmBin_mid,fliplr(tzDmmBin_mid)];
		
% 		tzDmmSpread = [tzDmmBinMin, fliplr(tzDmmBinMax)];
		tzDmmSpread = [tzDmmBin25p, fliplr(tzDmmBin75p)];

		f1 = fill(tzDmmBin_midFlip,tzDmmSpread,'b','FaceAlpha', 0.5);

		if strcmp(cntrLine,'Mean')
			m1 = plot(tzDmmBin_mid,tzDmmBinMean,'Color',[27/255 27/255 100/255],'LineWidth', 3);
		elseif strcmp(cntrLine,'Median')
			m1 = plot(tzDmmBin_mid,tzDmmBinMedian,'Color',[27/255 27/255 100/255],'LineWidth', 3);
		end	
	end
	
	if ~isempty(esrSprls)
		esrDmmBin_mid = (edgesTemp(1:end-1) + edgesTemp(2:end))/2;

		esrDmmTempBinIx = discretize(esrTemp,edgesTemp);

		esrDmmFlgRmvBin = [];
		esrDmmBinMean = [];
		esrDmmBinMedian = [];
		esrDmmBinMax = [];
		esrDmmBinMin = [];
		esrDmmBin25p = [];
		esrDmmBin75p = [];
		
		for ix=1:numBins
			esrDmmBinMem = (esrDmmTempBinIx == ix);
			esrBinDmm = esrDmm(esrDmmBinMem);

			if ~isempty(esrBinDmm)
				esrDmmBinMean = horzcat(esrDmmBinMean,nanmean(esrBinDmm));
				esrDmmBinMedian = horzcat(esrDmmBinMedian,nanmedian(esrBinDmm));
				esrDmmBinMax = horzcat(esrDmmBinMax,max(esrBinDmm));
				esrDmmBinMin = horzcat(esrDmmBinMin,min(esrBinDmm));
				esrDmmBin25p = horzcat(esrDmmBin25p,prctile(esrBinDmm,25));
				esrDmmBin75p = horzcat(esrDmmBin75p,prctile(esrBinDmm,75));
				
			else
				esrDmmFlgRmvBin = horzcat(esrDmmFlgRmvBin,ix);
			end
		end

		esrDmmBin_mid(esrDmmFlgRmvBin) = [];

		esrDmmBin_midFlip = [esrDmmBin_mid,fliplr(esrDmmBin_mid)];
		
% 		esrDmmSpread = [esrDmmBinMin, fliplr(esrDmmBinMax)];
		esrDmmSpread = [esrDmmBin25p, fliplr(esrDmmBin75p)];

		f2 = fill(esrDmmBin_midFlip,esrDmmSpread,'r','FaceAlpha', 0.5);

		if strcmp(cntrLine,'Mean')
			m2 = plot(esrDmmBin_mid,esrDmmBinMean,'Color',[100/255 27/255 27/255],'LineWidth', 3);
		elseif strcmp(cntrLine,'Median')
			m2 = plot(esrDmmBin_mid,esrDmmBinMedian,'Color','y','LineWidth', 3);
		end
	end
	
	
	if ~plotAvg
		title([flight ' - D_{mm} vs. Temp - ' cntrLine ' and Spread - ' probe ' - 1s Avg']);
	else
		title([flight ' - D_{mm} vs. Temp - ' cntrLine ' and Spread - ' probe ' - ' num2str(avgTime) 's Avg']);
	end
	set(gca,'XMinorTick','on','YMinorTick','on','XDir','reverse');
	if ~isempty(tempRangeAll)
		xlim(tempRangeAll);
	end
	if ~isempty(DmmRangeAll)
		ylim(DmmRangeAll)
	end
	ylabel('D_{mm} [mm]')
	xlabel(sprintf('Temp (%cC)', char(176)));
	set(findall(gcf,'-property','FontSize'),'FontSize',28);
	view([90 -90])
	
	legend([f1 f2],{'Transition Zone','Enhanced Stratiform'},'Location','eastoutside');
	

	if saveFigs
		set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		if ~plotAvg
			print([saveDir '/' probe '-Dmm-Temp_1s/' flight '_' probe '_Dmm-Temp_' cntrLine '-Spread-Fill_1s' outFileAppend],Ftype,Fres)
		else
			print([saveDir '/' probe '-Dmm-Temp_' num2str(avgTime) 's/' flight '_' probe '_Dmm-Temp_' cntrLine '-Spread-Fill_' num2str(avgTime) 's' outFileAppend],Ftype,Fres)
		end
	end
end