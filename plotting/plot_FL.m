% Make various plots using flight-level data from PECAN

clearvars; close all;


savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';
dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

flight = '20150709';

saveFigs	= 1;
noDisp		= 1;
Ftype		= '-dpdf';
% Ftype		= '-dpng';

if strcmp(Ftype,'-dpdf')
	Fres = '-painters'; % Use this for fully vectorized files
else
	Fres = Fres; % Use this for smaller files - saves figure with same resolution/size as displayed on screen
end

allSpirals = 1; % Plot every spiral for given flight - otherwise specified spirals given below are plotted
showMarkers = 0; % Plot markers at given (below) locations on line plots

% One plot for each spiral
plotTempRHAlt		= 1;
plotRHTemp			= 0;
plotTempAlt			= 0;
plotTempTdAlt		= 0;
plotTimeTemp		= 0;
plotTimeAlt			= 0;
plotTimeRH			= 0;

% Single plots with data from all spirals
plotAllRHTemp		= 0;
plotAllTempAlt		= 0;

% Single plots with mean/median, max/min spread and percentiles of all spirals
cntrLine = 'Median'; % Do we plot median or mean on spread plots?
% cntrLine = 'Mean';
plotRHTempSprd		= 0;
plotTempAltSprd		= 0;

% Single plots similar to spread plots above, but with color-filled
% spreads separated by MCS region
plotRHTempSprdF		= 0;

% Timeseries of temperature over whole flight
plotTempTS			= 0;

% Range of plotted variables over subset/all flights/spirals
% Useful if it is desired to directly compared plots
tempRangeAll = [-18.5 22];
RHrangeAll = [0 120];
altRangeAll = [1200 7500];

%% Get flight-specific parameters
startT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'startT');
endT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'endT');
mlBotTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlBotTime');
mlTopTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlTopTime');
mlBotTemp = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlBotTemp');
mlTopTemp = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlTopTemp');
% mcsStg = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mcsStg');
% sprlZone = nc_varget([dataPath '/' flight '_PECANparams.nc'],'sprlZone');
% 
% zoneUndef = (numel(unique(sprlZone)) == 1 && unique(sprlZone) == 'U');

fltLvlFile = ['/Users/danstechman/GoogleDrive/PECAN-Data/FlightLevelData/Processed/' flight '_FltLvl_Processed.mat'];

%% Load in processed flight-level data

load(fltLvlFile, '-regexp', '^(?!flight$|dataPath$|lat$|lon$)\w');

%% Determine which spirals to loop over
if ~allSpirals
	switch flight
		case '20150706'
			loopVctr = [1, 6]; % Spiral numbers to plot
			if showMarkers
				loopFLix = [278, 28];; %Index of 10 sec avg'd SD data where marker is desired
			end
			
		case '20150709'
			loopVctr = [2, 5, 9];
			if showMarkers
				loopFLix = [229, 19, 719]; %10 sec avg
			end
	end
else
	loopVctr = 1:length(startT);
end

%% Sort spirals based on spiral zone and MCS evolution
% tzSprls = find(sprlZone == 'T');
% esrSprls = find(sprlZone == 'S');
% raSprls = find(sprlZone == 'A');
% 
% formSprls = find(mcsStg == 'F');
% matureSprls = find(mcsStg == 'M');
% weakSprls = find(mcsStg == 'W');

%% Determine data indices for spiral start/end
startFLix = zeros(length(startT),1);
endFLix = zeros(length(startT),1);
mlTopIx = NaN(length(startT),1);
mlBotIx = NaN(length(startT),1);
for sprl=1:length(startT)
	startFLix(sprl) = find(time_secs_FL == startT(sprl));
	endFLix(sprl) = find(time_secs_FL == endT(sprl));
	if ~isnan(mlTopTime(sprl))
		mlTopIx(sprl) = find(time_secs_FL == mlTopTime(sprl));
	end
	if ~isnan(mlBotTime(sprl))
		mlBotIx(sprl) = find(time_secs_FL == mlBotTime(sprl));
	end
end

%% Create output directories
if saveFigs
    saveDir = [savePath flight];
	if (exist(saveDir, 'dir') ~= 7)
        mkdir(saveDir)
	end
	if ( (plotRHTemp || plotAllRHTemp || plotRHTempSprd || plotRHTempSprdF) && exist([saveDir '/RH-Temp'], 'dir') ~= 7 )
        mkdir([saveDir '/RH-Temp'])
	end
	if ( (plotTempAlt || plotAllTempAlt || plotTempAltSprd) && exist([saveDir '/Temp-Alt'], 'dir') ~= 7)
        mkdir([saveDir '/Temp-Alt'])
	end
	if ( plotTempTdAlt && exist([saveDir '/TempTd-Alt'], 'dir') ~= 7)
        mkdir([saveDir '/TempTd-Alt'])
	end
	if (plotTempRHAlt && exist([saveDir '/Temp-RH-Alt'], 'dir') ~= 7)
        mkdir([saveDir '/Temp-RH-Alt'])
	end
end


%% Plot creation

%%% RH & Temp vs. Alt - each spiral
% Each on their own x-axis (using plotxx function from the MATLAB file exchange). 
if plotTempRHAlt
	j = 1; % Only used if showMarkers is True
	for ix = loopVctr
		if saveFigs && noDisp
            figure('visible','off','Position', [10,10,1200,1200]);
		else
            figure('Position', [10,10,1200,1200]);
		end
		
		TAtmp = TA(startFLix(ix):endFLix(ix));
		AltTmp = Alt(startFLix(ix):endFLix(ix));
		RHtmp = RH_hybrid(startFLix(ix):endFLix(ix));
		
        [hAx,hTemp,hRH] = plotxx(TAtmp,AltTmp,RHtmp,AltTmp,...
            {sprintf('Temperature (%cC)', char(176)),'Relative Humidity (%)'},{'Altitude (m MSL)'});
        hTemp.Color = 'r';
        hRH.Color = 'b';
        set(hAx(1),'xlim',tempRangeAll);
        set(hAx(2),'xlim',RHrangeAll);
        set(hAx(1),'ylim',altRangeAll);
        set(hAx(2),'ylim',altRangeAll);
        hAx(1).XColor = 'r';
        hAx(2).XColor = 'b'; 
		hTemp.LineWidth = 3;
		hRH.LineWidth = 3;

        title([flight ' - Spiral ' num2str(ix) ' - RH & Temp vs. Alt']);
		hold(hAx(1),'on')
		hold(hAx(2),'on')

        plot(hAx(1),[0 0],ylim,'r--','LineWidth',1); % Vertical line at 0 deg C
		plot(hAx(2),[100 100],ylim,'b--','LineWidth',1); % Vertical line at 100% RH
		if (showMarkers && ~allSpirals)
			plot(hAx(1),TAtmp(loopFLix(j)),AltTmp(loopFLix(j)),'Marker','o','MarkerFaceColor','black',...
				'MarkerEdgeColor','black','Markersize',10);
			plot(hAx(2),RHtmp(loopFLix(j)),AltTmp(loopFLix(j)),'Marker','o','MarkerFaceColor','black',...
				'MarkerEdgeColor','black','Markersize',10);
			j = j+1;
		end
		
		txtLoc = (tempRangeAll(2)-tempRangeAll(1))*0.5;
		if ~isnan(mlTopIx(ix))
			plot(hAx(1),tempRangeAll,[1 1]*Alt(mlTopIx(ix)),'k--')
			topStr = sprintf('%.3f %cC',mlTopTemp(ix),char(176));
			tMT = text(hAx(1),txtLoc,Alt(mlTopIx(ix)),topStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
		end
		if ~isnan(mlBotIx(ix))
			plot(hAx(1),tempRangeAll,[1 1]*Alt(mlBotIx(ix)),'k--')
			botStr = sprintf('%.3f %cC',mlBotTemp(ix),char(176));
			tMB = text(hAx(1),txtLoc,Alt(mlBotIx(ix)),botStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
		end
		
        set(gca,'ygrid','on')
        set(findall(gcf,'-property','FontSize'),'FontSize',28)
		
		if ~isnan(mlBotIx(ix))
			set(tMB,'FontSize',14);
		end
		if ~isnan(mlTopIx(ix))
			set(tMT,'FontSize',14);
		end
		
		tightfig(gcf);
		
		if saveFigs
            set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/Temp-RH-Alt/' flight '_TempRHAlt_S' num2str(ix)],Ftype,Fres)
		end
	end
end

if plotRHTemp
	j = 1;
	for ix = loopVctr
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1500,1500]);
		else
			figure('Position', [10,10,1500,1500]);
		end
		
		TAtmp = TA(startFLix(ix):endFLix(ix));
		RHtmp = RH_hybrid(startFLix(ix):endFLix(ix));
		
        plot(RHtmp,TAtmp,'b-');
		
		if (showMarkers && ~allSpirals)
			hold on;
			plot(RHtmp(loopFLix(j)),TAtmp(loopFLix(j)),'Marker','o','MarkerFaceColor','black',...
				'MarkerEdgeColor','black','Markersize',10);
			j = j+1;
		end
		
        set(gca,'XMinorTick','on','YMinorTick','on','YDir','reverse');
        grid
        title([flight ' - ' 'Spiral ' num2str(ix) ' - RH vs. Temp']);
        xlabel('RH (%)');
        ylabel('Temperature (deg C)');
        xlim(RHrangeAll);
        ylim(tempRangeAll);
        set(findall(gcf,'-property','FontSize'),'FontSize',28);
        
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/RH-Temp/' flight '_RH-Temp_S' num2str(ix)],Ftype,Fres)
		end

	end
end

if plotTempAlt
	j = 1;
	for ix = loopVctr
        if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1500,1500]);
		else
			figure('Position', [10,10,1500,1500]);
		end
		
		TAtmp = TA(startFLix(ix):endFLix(ix));
		AltTmp = Alt(startFLix(ix):endFLix(ix));
		
        plot(TAtmp,AltTmp,'b-');
		
		if (showMarkers && ~allSpirals)
			hold on;
			plot(TAtmp(loopFLix(j)),AltTmp(loopFLix(j)),'Marker','o','MarkerFaceColor','black',...
				'MarkerEdgeColor','black','Markersize',10);
			j = j+1;
		end
		
        set(gca,'XMinorTick','on','YMinorTick','on');
        grid
        title([flight ' - ' 'Spiral ' num2str(ix) ' - Temp vs. Alt']);
        xlabel('Temperature (deg C)');
        ylabel('Altitude (m MSL)');
        xlim(tempRangeAll);
		ylim(altRangeAll);
        set(findall(gcf,'-property','FontSize'),'FontSize',28);
        
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/Temp-Alt/' flight '_Temp-Alt_S' num2str(ix)],Ftype,Fres)
		end

	end
end

if plotTempTdAlt
	j = 1;
	for ix = loopVctr
        if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1500,1500]);
		else
			figure('Position', [10,10,1500,1500]);
		end
		
		TAtmp = TA(startFLix(ix):endFLix(ix));
		TdTmp = TD(startFLix(ix):endFLix(ix));
		AltTmp = Alt(startFLix(ix):endFLix(ix));
		
        plot(TAtmp,AltTmp,'r-');
		hold on
		plot(TdTmp,AltTmp,'b-x');
		
		if (showMarkers && ~allSpirals)
			hold on;
			plot(TAtmp(loopFLix(j)),AltTmp(loopFLix(j)),'Marker','o','MarkerFaceColor','black',...
				'MarkerEdgeColor','black','Markersize',10);
			j = j+1;
		end
		
        set(gca,'XMinorTick','on','YMinorTick','on');
        grid
        title([flight ' - ' 'Spiral ' num2str(ix) ' - Temp & Td vs. Alt']);
        xlabel('Temperature (deg C)');
        ylabel('Altitude (m MSL)');
%         xlim(tempRangeAll);
		ylim(altRangeAll);
        set(findall(gcf,'-property','FontSize'),'FontSize',28);
        
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/TempTd-Alt/' flight '_TempTd-Alt_S' num2str(ix)],Ftype,Fres)
		end

	end
end

if plotTimeTemp
	j = 1;
	for ix = loopVctr
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1500,1500]);
		else
			figure('Position', [10,10,1500,1500]);
		end
		
		TAtmp = TA(startFLix(ix):endFLix(ix));
		RHtmp = RH_hybrid(startFLix(ix):endFLix(ix));
		
        plot(RHtmp,TAtmp,'b-');
		
		if (showMarkers && ~allSpirals)
			hold on;
			plot(RHtmp(loopFLix(j)),TAtmp(loopFLix(j)),'Marker','o','MarkerFaceColor','black',...
				'MarkerEdgeColor','black','Markersize',10);
			j = j+1;
		end
		
        set(gca,'XMinorTick','on','YMinorTick','on','YDir','reverse');
        grid
        title([flight ' - ' 'Spiral ' num2str(ix) ' - RH vs. Temp']);
        xlabel('RH (%)');
        ylabel('Temperature (deg C)');
        xlim(RHrangeAll);
        ylim(tempRangeAll);
        set(findall(gcf,'-property','FontSize'),'FontSize',28);
        
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/RH-Temp/' flight '_RH-Temp_S' num2str(ix)],Ftype,Fres)
		end

	end
end

if plotAllRHTemp
    
	colors = varycolor(length(startT));

    RHsprl = RH_hybrid(startFLix(1):endFLix(1));
    tempSprl = TA(startFLix(1):endFLix(1));
    
	if saveFigs && noDisp
		figure('visible','off','Position', [10,10,1500,1500]);
	else
		figure('Position', [10,10,1500,1500]);
	end

	% Plot the first spiral profile
	plot(RHsprl,tempSprl,'Color',colors(1,:));
	set(gca,'XMinorTick','on','YMinorTick','on','YDir','reverse');
	grid
	title([flight ' - RH vs. Temp - All Spirals']);
	xlabel('RH (%)');
	ylabel('Temperature (deg C)');
	set(findall(gcf,'-property','FontSize'),'FontSize',28);
	hold on

	% Determine and plot the subsequent spiral profiles
	for ix = 2:length(startT)
        RHsprl = RH_hybrid(startFLix(ix):endFLix(ix));
        tempSprl = TA(startFLix(ix):endFLix(ix));
		
		plot(RHsprl,tempSprl,'Color',colors(ix,:));
	end

	% Create legends depending on which flight it is
	if strcmp(flight,'20150706')
		legend('S1','S2','S3','S4','S5','S6','S7','Location','eastoutside');
	elseif strcmp(flight,'20150709')
		legend('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10',...
			'S11','S12','S13','S14','S15','S16','Location','eastoutside');
	end

    
	if saveFigs
		set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		print([saveDir '/RH-Temp/' flight '_RH-Temp_All'],Ftype,Fres)
	end
   
end

if plotAllTempAlt
    
	colors = varycolor(length(startT));
	
    % Create a matrix with columns representing (spiral#, Temp, Altitude)
	% and fill the values for the first spiral
    altSprl = Alt(startFLix(1):endFLix(1));
    tempSprl = TA(startFLix(1):endFLix(1));
    
	if saveFigs && noDisp
		figure('visible','off','Position', [10,10,1500,1500]);
	else
		figure('Position', [10,10,1500,1500]);
	end

	% Plot the first spiral profile
	plot(tempSprl,altSprl,'Color',colors(1,:));
	set(gca,'XMinorTick','on','YMinorTick','on');
	grid
	title([flight ' - Temp vs. Altitude - All Spirals']);
	xlabel('Temperature (deg C)');
	ylabel('Altitude (m)');
	set(findall(gcf,'-property','FontSize'),'FontSize',28);
	xlim(tempRangeAll);
	ylim(altRangeAll);
	hold on

	% Determine and plot the subsequent spiral profiles
    for ix = 2:length(startT)
        altSprl = Alt(startFLix(ix):endFLix(ix));
        tempSprl = TA(startFLix(ix):endFLix(ix));
		
		plot(tempSprl,altSprl,'Color',colors(ix,:));
	end

	% Create legends depending on which flight it is
	if strcmp(flight,'20150706')
		legend('S1','S2','S3','S4','S5','S6','S7','Location','eastoutside');
	elseif strcmp(flight,'20150709')
		legend('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10',...
			'S11','S12','S13','S14','S15','S16','Location','eastoutside');
	end

    
	if saveFigs
		set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		print([saveDir '/Temp-Alt/' flight '_Temp-Alt_All'],Ftype,Fres)
	end
end

if plotTempTS
	if saveFigs && noDisp
		figure('visible','off','Position', [10,10,1500,1500]);
	else
		figure('Position', [10,10,1500,1500]);
	end
	
    serialTime = time_secs_FL/3600/24;
    plot(serialTime,TA,'r');
    hold on
    
	for ix = 1:length(startT)
        plot(serialTime(startFLix(ix):endFLix(ix)),TA(startFLix(ix):endFLix(ix)),'b');
    end
    title([flight ' - Temperature']);
    set(gca,'ygrid','on','YDir','reverse')
    set(findall(gcf,'-property','FontSize'),'FontSize',28)
	
	% Plot 20 min before first spiral through 20 min past last spiral
	pltStartT = (startT(1)/3600/24)-(1200/3600/24);
	pltEndT = (endT(end)/3600/24)+(1200/3600/24);
    xlim([pltStartT pltEndT]) 
    set(gca,'XTick',(pltStartT:0.01:pltEndT))
    datetick('x','HH:MM:SS','keeplimits','keepticks');
    set(gca,'XTickLabelRotation',45);
    set(gca,'XMinorTick','on','YMinorTick','on')
    
    if saveFigs
        set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		print([saveDir '/' flight '_TempTimeSeries'],Ftype,Fres)
    end
end

if plotRHTempSprd
	
    RHsprl = RH_hybrid(startFLix(1):endFLix(1));
    tempSprl = TA(startFLix(1):endFLix(1));
	
	allSprlRH = RHsprl;
	allSprlTemp = tempSprl;

	% Determine the subsequent spiral profiles
	for ix = 2:length(startT)
        RHsprl = RH_hybrid(startFLix(ix):endFLix(ix));
        tempSprl = TA(startFLix(ix):endFLix(ix));
		
		allSprlRH = vertcat(allSprlRH, RHsprl);
		allSprlTemp = vertcat(allSprlTemp, tempSprl);
	end

	edgesTemp = (-15.125:.25:20.125);
	numBins = length(edgesTemp)-1;

	bin_min = edgesTemp(1:end-1); 
	bin_max = edgesTemp(2:end); 
	bin_mid = (bin_min+bin_max)/2;
	
	whichBinTemp = discretize(allSprlTemp,edgesTemp);

	for ix=1:numBins
		flagBinMemTemp = (whichBinTemp == ix);
		binRH = allSprlRH(flagBinMemTemp);
		if (isempty(binRH))
			binMean(ix) = NaN;
			binMedian(ix) = NaN;
			binMax(ix) = NaN;
			binMin(ix) = NaN;
			bin25pct(ix) = NaN;
			bin75pct(ix) = NaN;
		else
			binMean(ix) = nanmean(binRH);
			binMedian(ix) = nanmedian(binRH);
			binMax(ix) = max(binRH);
			binMin(ix) = min(binRH);
			bin25pct(ix) = prctile(binRH,25);
			bin75pct(ix) = prctile(binRH,75);
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
	xlim(tempRangeAll)
	ylim(RHrangeAll)
	view([90 -90])
	grid
	title([flight ' - RH vs. Temp - ' cntrLine ' and Spread all Spirals']);
	ylabel('RH (%)');
	xlabel('Temperature (deg C)');
	set(findall(gcf,'-property','FontSize'),'FontSize',28);

	if saveFigs
		set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		print([saveDir '/RH-Temp/' flight '_RH-Temp_' cntrLine '-Spread'],Ftype,Fres)
	end
end

if plotTempAltSprd
	
    altSprl = Alt(startFLix(1):endFLix(1));
    tempSprl = TA(startFLix(1):endFLix(1));
	
	allSprlAlt = altSprl;
	allSprlTemp = tempSprl;

	% Determine the subsequent spiral profiles
	for ix = 2:length(startT)
        altSprl = Alt(startFLix(ix):endFLix(ix));
        tempSprl = TA(startFLix(ix):endFLix(ix));
		
		allSprlAlt = vertcat(allSprlAlt, altSprl);
		allSprlTemp = vertcat(allSprlTemp, tempSprl);
	end
	
	edgesAlt = (1887.5:25:6812.5); % Bin edges for altitude/altitude proxy (temp)
	numBins = length(edgesAlt)-1;
	bin_min = edgesAlt(1:end-1); 
	bin_max = edgesAlt(2:end); 
	bin_mid = (bin_min+bin_max)/2;
	
	whichBinAlt = discretize(allSprlAlt,edgesAlt);

	for ix=1:numBins
		flagBinMemAlt = (whichBinAlt == ix);
		binTemps = allSprlTemp(flagBinMemAlt);
		if (isempty(binTemps))
			binMean(ix) = NaN;
			binMedian(ix) = NaN;
			binMax(ix) = NaN;
			binMin(ix) = NaN;
			bin25pct(ix) = NaN;
			bin75pct(ix) = NaN;
		else
			binMean(ix) = nanmean(binTemps);
			binMedian(ix) = nanmedian(binTemps);
			binMax(ix) = max(binTemps);
			binMin(ix) = min(binTemps);
			bin25pct(ix) = prctile(binTemps,25);
			bin75pct(ix) = prctile(binTemps,75);
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
% 	set(gca,'xtick', 0:6:length(edgesAlt), 'xticklabel', 1900:140:6800, 'XMinorTick','on',...
% 		'YMinorTick','on');
	set(gca,'XMinorTick','on','YMinorTick','on')
	xlim(altRangeAll)
	ylim(tempRangeAll)
	view([90 -90])
	grid
	title([flight ' - Temp vs. Altitude - ' cntrLine ' and Spread all Spirals']);
	xlabel('Altitude (m MSL)');
	ylabel('Temperature (deg C)');
	set(findall(gcf,'-property','FontSize'),'FontSize',28);

	if saveFigs
		set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		print([saveDir '/Temp-Alt/' flight '_Temp-Alt_' cntrLine '-Spread'],Ftype,Fres)
	end
end

if plotRHTempSprdF
	if zoneUndef
		disp('MCS zones have not been defined in the PECAN parameters files - RH-Temp Spread not plotted')
	else
		
		if ~isempty(tzSprls)
			
			
			RHSprl = RH_hybrid(startFLix(tzSprls(1)):endFLix(tzSprls(1)));
			tempCsprl = TA(startFLix(tzSprls(1)):endFLix(tzSprls(1)));
			
			tzRH = RHSprl;
			tzTemp = tempCsprl;
			
			% Determine the subsequent spiral profiles
			if length(tzSprls) > 1
				for ix = 2:length(tzSprls)
					RHSprl = RH_hybrid(startFLix(tzSprls(ix)):endFLix(tzSprls(ix)));
					tempCsprl = TA(startFLix(tzSprls(ix)):endFLix(tzSprls(ix)));
					
					tzRH = vertcat(tzRH, RHSprl);
					tzTemp = vertcat(tzTemp, tempCsprl);
				end
			end
			
			tzRH_orig = tzRH;
			tzRH(isnan(tzRH)) = 0;
		end
		
		
		if ~isempty(esrSprls)
			
			RHSprl = RH_hybrid(startFLix(esrSprls(1)):endFLix(esrSprls(1)));
			tempCsprl = TA(startFLix(esrSprls(1)):endFLix(esrSprls(1)));
			
			esrRH = RHSprl;
			esrTemp = tempCsprl;
			
			% Determine the subsequent spiral profiles
			if length(esrSprls) > 1
				for ix = 2:length(esrSprls)
					RHSprl = RH_hybrid(startFLix(esrSprls(ix)):endFLix(esrSprls(ix)));
					tempCsprl = TA(startFLix(esrSprls(ix)):endFLix(esrSprls(ix)));
					
					esrRH = vertcat(esrRH, RHSprl);
					esrTemp = vertcat(esrTemp, tempCsprl);
				end
			end
			
			esrRH_orig = esrRH;
			% 		esrRH(isnan(esrRH)) = 0;
		end
		
		
		if ~isempty(raSprls)
			
			RHSprl = RH_hybrid(startFLix(raSprls(1)):endFLix(raSprls(1)));
			tempCsprl = TA(startFLix(raSprls(1)):endFLix(raSprls(1)));
			
			raRH = RHSprl;
			raTemp = tempCsprl;
			
			% Determine the subsequent spiral profiles
			if length(raSprls) > 1
				for ix = 2:length(raSprls)
					RHSprl = RH_hybrid(startFLix(raSprls(ix)):endFLix(raSprls(ix)));
					tempCsprl = TA(startFLix(raSprls(ix)):endFLix(raSprls(ix)));
					
					raRH = vertcat(raRH, RHSprl);
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
			tzRHBin_mid = (edgesTemp(1:end-1) + edgesTemp(2:end))/2;
			
			tzRHTempBinIx = discretize(tzTemp,edgesTemp);
			
			tzRHFlgRmvBin = [];
			tzRHBinMean = [];
			tzRHBinMedian = [];
			tzRHBinMax = [];
			tzRHBinMin = [];
			tzRHBin25p = [];
			tzRHBin75p = [];
			
			for ix=1:numBins
				tzRHBinMem = (tzRHTempBinIx == ix);
				tzBinRH = tzRH(tzRHBinMem);
				tzBinRH_orig = tzRH_orig(tzRHBinMem);
				
				if ~isempty(tzBinRH)
					tzRHBinMean = horzcat(tzRHBinMean,nanmean(tzBinRH_orig));
					tzRHBinMedian = horzcat(tzRHBinMedian,nanmedian(tzBinRH_orig));
					tzRHBinMax = horzcat(tzRHBinMax,max(tzBinRH));
					tzRHBinMin = horzcat(tzRHBinMin,min(tzBinRH));
					tzRHBin25p = horzcat(tzRHBin25p,prctile(tzBinRH,25));
					tzRHBin75p = horzcat(tzRHBin75p,prctile(tzBinRH,75));
					
				else
					tzRHFlgRmvBin = horzcat(tzRHFlgRmvBin,ix);
				end
			end
			
			tzRHBin_mid(tzRHFlgRmvBin) = [];
			
			tzRHBin_midFlip = [tzRHBin_mid,fliplr(tzRHBin_mid)];
			% 		tzRHSpread = [tzRHBinMin, fliplr(tzRHBinMax)];
			tzRHSpread = [tzRHBin25p, fliplr(tzRHBin75p)];
			
			f1 = fill(tzRHBin_midFlip,tzRHSpread,'b','FaceAlpha', 0.5);
			
			if strcmp(cntrLine,'Mean')
				m1 = plot(tzRHBin_mid,tzRHBinMean,'Color',[27/255 27/255 100/255],'LineWidth', 3);
			elseif strcmp(cntrLine,'Median')
				m1 = plot(tzRHBin_mid,tzRHBinMedian,'Color',[27/255 27/255 100/255],'LineWidth', 3);
			end
		end
		
		if ~isempty(esrSprls)
			esrRHBin_mid = (edgesTemp(1:end-1) + edgesTemp(2:end))/2;
			
			esrRHTempBinIx = discretize(esrTemp,edgesTemp);
			
			esrRHFlgRmvBin = [];
			esrRHBinMean = [];
			esrRHBinMedian = [];
			esrRHBinMax = [];
			esrRHBinMin = [];
			esrRHBin25p = [];
			esrRHBin75p = [];
			
			for ix=1:numBins
				esrRHBinMem = (esrRHTempBinIx == ix);
				esrBinRH = esrRH(esrRHBinMem);
				% 			esrBinRH_orig = esrRH_orig(tzRHBinMem);
				
				if ~isempty(esrBinRH)
					esrRHBinMean = horzcat(esrRHBinMean,nanmean(esrBinRH));
					esrRHBinMedian = horzcat(esrRHBinMedian,nanmedian(esrBinRH));
					esrRHBinMax = horzcat(esrRHBinMax,max(esrBinRH));
					esrRHBinMin = horzcat(esrRHBinMin,min(esrBinRH));
					esrRHBin25p = horzcat(esrRHBin25p,prctile(esrBinRH,25));
					esrRHBin75p = horzcat(esrRHBin75p,prctile(esrBinRH,75));
					
				else
					esrRHFlgRmvBin = horzcat(esrRHFlgRmvBin,ix);
				end
			end
			
			esrRHBin_mid(esrRHFlgRmvBin) = [];
			
			esrRHBin_midFlip = [esrRHBin_mid,fliplr(esrRHBin_mid)];
			% 		esrRHSpread = [esrRHBinMin, fliplr(esrRHBinMax)];
			esrRHSpread = [esrRHBin25p, fliplr(esrRHBin75p)];
			
			f2 = fill(esrRHBin_midFlip,esrRHSpread,'r','FaceAlpha', 0.5);
			
			if strcmp(cntrLine,'Mean')
				m2 = plot(esrRHBin_mid,esrRHBinMean,'Color',[100/255 27/255 27/255],'LineWidth', 3);
			elseif strcmp(cntrLine,'Median')
				m2 = plot(esrRHBin_mid,esrRHBinMedian,'Color','y','LineWidth', 3);
			end
			%
			
		end
		
		
		title([flight ' - RH vs. Temp - ' cntrLine ' and Spread']);
		set(gca,'XMinorTick','on','YMinorTick','on','XDir','reverse');
		if ~isempty(tempRangeAll)
			xlim(tempRangeAll);
		end
		if ~isempty(RHrangeAll)
			ylim(RHrangeAll)
		end
		ylabel('RH [%]')
		xlabel('Temperature (deg C)');
		set(findall(gcf,'-property','FontSize'),'FontSize',28);
		view([90 -90])
		
		legend([f1 f2],{'Transition Zone','Enhanced Stratiform'},'Location','eastoutside');
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/RH-Temp/' flight '_RH-Temp_' cntrLine '-Spread-Fill'],Ftype,Fres)
		end
	end
end
