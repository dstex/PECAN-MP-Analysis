clearvars; close all;

%% Specify various plotting/calculation parameters
flight = '20150706';

avgTime = 10;

fileIdStr = ['_Fit-CIP_' num2str(avgTime) 'secAvg_1.2cm'];
titleIdStr = '';

% Standard plots, but using hybrid (CIP obs + CIP extended) SDs
plotND				= 0;
plotMD				= 0;

plotNDall			= 0;
plotMDall			= 0;

plotNDtemp			= 0;
plotMDtemp			= 0;

plotTWCextndRatio	= 0;

getLims = 0; % If true, run without plotting and print values to assist in setting axis limits
plotNDtempBinned	= 0;
plotMDtempBinned	= 0;

% diamLim = [0.1 5];
diamLim = [0.1 inf];


saveFigs	= 0;
noDisp		= 0;
Ftype		= '-dpdf';
% Ftype		= '-dpng';

if strcmp(Ftype,'-dpdf')
	Fres = '-painters'; % Use this for fully vectorized files
else
	Fres = '-r0'; % Use this for smaller files - saves figure with same resolution/size as displayed on screen
end


savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';
dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

%% Load in struct of all original SD data
if avgTime == 1
	sDistFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.CIP.10secAvg.mat']; % 1-sec averages are saved in all averaged output files
else
	sDistFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.CIP.' num2str(avgTime) 'secAvg.mat'];
end

sDistF = load(sDistFile);


%% Create directories to save plots in if they don't already exist
if saveFigs
    saveDir = [savePath flight '/CIP-Fitting'];
	if (exist(saveDir, 'dir') ~= 7)
        mkdir(saveDir)
	end

	if ((plotND || plotNDall) && exist([saveDir '/CIP-ND_' num2str(avgTime) 's'], 'dir') ~= 7)
		mkdir([saveDir '/CIP-ND_' num2str(avgTime) 's'])
	end
	if ((plotMD || plotMDall) && exist([saveDir '/CIP-MD_' num2str(avgTime) 's'], 'dir') ~= 7)
		mkdir([saveDir '/CIP-MD_' num2str(avgTime) 's'])
	end
	if (plotNDtemp && exist([saveDir '/CIP-ND-Temp_' num2str(avgTime) 's'], 'dir') ~= 7)
		mkdir([saveDir '/CIP-ND-Temp_' num2str(avgTime) 's'])
	end
	if (plotMDtemp && exist([saveDir '/CIP-MD-Temp_' num2str(avgTime) 's'], 'dir') ~= 7)
		mkdir([saveDir '/CIP-MD-Temp_' num2str(avgTime) 's'])
	end
	if (plotTWCextndRatio && exist([saveDir '/CIP-TWCextndRatio_' num2str(avgTime) 's'], 'dir') ~= 7)
		mkdir([saveDir '/CIP-TWCextndRatio_' num2str(avgTime) 's'])
	end
	if (plotNDtempBinned && exist([saveDir '/CIP-ND-TempBinned'], 'dir') ~= 7)
		mkdir([saveDir '/CIP-ND-TempBinned'])
	end
	if (plotMDtempBinned && exist([saveDir '/CIP-MD-TempBinned'], 'dir') ~= 7)
		mkdir([saveDir '/CIP-MD-TempBinned'])
	end
end

load([dataPath 'mp-data/' flight '/sDist/' flight fileIdStr '.mat']);

loopVctr = 1:length(sprlNames);
% loopVctr = [4,7];

%% Specify flight-specific plotting parameters
tempRangeAll = [-18.5 22];

switch flight
	case '20150617'
		NDLim = [1e-10 3.5e-5];
		NDtempBinLim = [1e-8 30];
		MDtempBinLim = [];
		NDLogLim = [-5 2];
		MDLim = [];
		MDLogLim = [-10 -4];
		NtRangeAll = [];
		TWCrangeAll = [];
		DmmRangeAll = [];
		NtRangeFill = [];
		TWCRangeFill = [];
	case '20150620'
		NDLim = [];
		NDtempBinLim = [1e-10 20];
		MDtempBinLim = [1e-10 2.5e-5];
		NDLogLim = [-5 2];
		MDLim = [];
		MDLogLim = [-10 -4];
		NtRangeAll = [];
		TWCrangeAll = [];
		DmmRangeAll = [];
		NtRangeFill = [];
		TWCRangeFill = [];
	case '20150701'
		NDLim = [];
		NDtempBinLim = [1e-8 5];
		MDtempBinLim = [1e-12 8e-6];
		NDLogLim = [-5 2];
		MDLim = [];
		MDLogLim = [-10 -4];
		NtRangeAll = [];
		TWCrangeAll = [];
		DmmRangeAll = [];
		NtRangeFill = [];
		TWCRangeFill = [];
	case '20150702'
		NDLim = [];
		NDtempBinLim = [1e-8 5];
		MDtempBinLim = [1e-10 3e-6];
		NDLogLim = [-5 2];
		MDLim = [];
		MDLogLim = [-10 -4];
		NtRangeAll = [];
		TWCrangeAll = [];
		DmmRangeAll = [];
		NtRangeFill = [];
		TWCRangeFill = [];
	case '20150706'
		NDLim = [1e-10 10];
		NDtempBinLim = [1e-10 20];
		MDtempBinLim = [1e-10 2.5e-5];
		NDLogLim = [-5 2];
		MDLim = [1e-10 1e-4];
		MDLogLim = [-10 -4];
		NtRangeAll = [];
		TWCrangeAll = [];
		DmmRangeAll = [];
		NtRangeFill = [];
		TWCRangeFill = [];
	case '20150709'
		NDLim = [];
		NDtempBinLim = [1e-8 5];
		MDtempBinLim = [1e-10 8e-6];
		NDLogLim = [-5 2];
		MDLim = [];
		MDLogLim = [-10 -4];
		NtRangeAll = [];
		TWCrangeAll = [];
		DmmRangeAll = [];
		NtRangeFill = [];
		TWCRangeFill = [];
end


%% Standard plots
if plotND
    for ix = loopVctr
		
		cipConc = cipConc_hybrid_igfWhl.(sprlNames{ix});
		
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		
		stairs(cipExt_binMin(1:length(cip_binMin))*10, cipConc(1:length(cip_binMin))', 'Color',[0 0 0.3], 'LineWidth', 2);
		hold on
		stairs(cipExt_binMin(length(cip_binMin):end)*10, cipConc(length(cip_binMin):end)','b', 'LineWidth', 2);

		title([flight ' - Spiral ' num2str(ix) ' - CIP (Extended) - ' num2str(avgTime) 's Avg']);

		xlabel('D [mm]');
		ylabel('N(D) [cm^{-4}]');
		set(gca,'Yscale','log');
		set(gca,'Xscale','log');
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
			print([saveDir '/CIP-ND_' num2str(avgTime) 's/' flight '_CIP-Fit_ND_' num2str(avgTime) 's_S' num2str(ix)],Ftype,Fres)
        end
    end
end

if plotMD
    for ix = loopVctr
		
		meanMass = cipMass_hybrid_igfWhl.(sprlNames{ix});
		
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		set(gcf,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
		
		yyaxis left
		stairs(cipExt_binMin(1:length(cip_binMin))*10, meanMass(1:length(cip_binMin))', 'Color',[0 0 0.3], 'LineWidth', 2);
		hold on
		stairs(cipExt_binMin(length(cip_binMin):end)*10, meanMass(length(cip_binMin):end)','b-', 'LineWidth', 2);
		
		title([flight ' - Spiral ' num2str(ix) ' - CIP (Extended) - ' num2str(avgTime) 's Avg']);
		
		xlabel('D [mm]');
		ylabel('M_{twc}(D) [g cm^{-4}]');
		set(gca,'Yscale','log');
		set(gca,'Xscale','log');
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
		plot(cipExt_binMin(1:length(cip_binMin))*10,massCDF(1:length(cip_binMin))*100,'Color',[0.3 0 0],'LineWidth',2);
		hold on
		plot(cipExt_binMin(length(cip_binMin):end)*10,massCDF(length(cip_binMin):end)*100,'r-','LineWidth',2);
		ylabel('Mass_{twc}(D) CDF [%]');
		
		
        if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/CIP-MD_' num2str(avgTime) 's/' flight '_CIP-Fit_MD_' num2str(avgTime) 's_S' num2str(ix)],Ftype,Fres)
        end
    end
end

if plotNDall
    for ix = loopVctr
		
		cipConc = cipConc_hybrid_igf.(sprlNames{ix});
		
		colors = varycolor(size(cipConc,1));
		
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		
		hold on
		
		for icx=1:size(cipConc,1)
			plot(cipExt_binMax(1:length(cip_binMax))*10, cipConc(icx,1:length(cip_binMax))','-','Color',colors(icx,:), 'LineWidth', 2);
			
			plot(cipExt_binMax(length(cip_binMax):end)*10, cipConc(icx,length(cip_binMax):end)','--','Color',colors(icx,:), 'LineWidth', 2);
		end
		title([flight ' - Spiral ' num2str(ix) ' - CIP (Extended) - ' num2str(avgTime) 's Avg']);

		xlabel('D [mm]');
		ylabel('N(D) [cm^{-4}]');
		set(gca,'Yscale','log');
		if ~isempty(diamLim)
			xlim(diamLim);
		end
% 		if ~isempty(NDLim)
% 			ylim(NDLim);
% 		end
		set(gca,'XMinorTick','on','YMinorTick','on');
		set(findall(gcf,'-property','FontSize'),'FontSize',28)
		grid


        if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/CIP-ND_' num2str(avgTime) 's/' flight '_CIP-Fit_ND_' num2str(avgTime) 's_S' num2str(ix) '_all'],Ftype,Fres)
        end
    end
end

if plotMDall
    for ix = loopVctr
		
		cipMass = cipMass_hybrid_igf.(sprlNames{ix});
		
		colors = varycolor(size(cipMass,1));
		
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		
		hold on
		
		for icx=1:size(cipMass,1)
			plot(cipExt_binMax(1:length(cip_binMax))*10, cipMass(icx,1:length(cip_binMax))','-','Color',colors(icx,:), 'LineWidth', 2);
			
			plot(cipExt_binMax(length(cip_binMax):end)*10, cipMass(icx,length(cip_binMax):end)','--','Color',colors(icx,:), 'LineWidth', 2);
		end
		title([flight ' - Spiral ' num2str(ix) ' - CIP (Extended) - ' num2str(avgTime) 's Avg']);

		xlabel('D [mm]');
		ylabel('M(D) [g cm^{-4}]');
		set(gca,'Yscale','log');
		if ~isempty(diamLim)
			xlim(diamLim);
		end

		set(gca,'XMinorTick','on','YMinorTick','on');
		set(findall(gcf,'-property','FontSize'),'FontSize',28)
		grid


        if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/CIP-MD_' num2str(avgTime) 's/' flight '_CIP-Fit_MD_' num2str(avgTime) 's_S' num2str(ix) '_all'],Ftype,Fres)
        end
    end
end

if plotNDtemp
    for ix = loopVctr
		
		cipConc = cipConc_hybrid_igf.(sprlNames{ix});
		
		tempCsprl = sDistF.tempC_orig.(sprlNames{ix});
		time_fl = sDistF.time_secsFL_orig.(sprlNames{ix});
		
		if avgTime == 1
			time_secs = sDistF.time_secs_orig.(sprlNames{ix});	
		else
			time_secs = sDistF.time_secs_avg.(sprlNames{ix});	
		end
		

        if saveFigs && noDisp
            figure('visible','off','Position', [10,10,1500,1000]);
        else
            figure('Position', [10,10,1500,1000]);
		end
		
		set(gcf,'defaultAxesColorOrder',[0 0 0; 0 0 0]);

        yyaxis left
		ax = gca;
        contourf(cipExt_binMid*10,time_secs/24/3600,log10(cipConc),NDLogLim(1):0.1:NDLogLim(2),'LineColor','none');
        xlabel('D [mm]');
        ylabel('Time');
		set(ax,'XMinorTick','on');
		colormap(jetmod);
        c=colorbar;
		set(c,'Location','southoutside');
		ylabel(c,'log_{10}N(D) [cm^{-4}]');
		set(ax, 'CLim', NDLogLim);
		set(ax,'XLim',diamLim);
		
		% Plot ML top/bottom locations and annotate with temp
		hold on
		plot(diamLim, [1 1]*mlTopTime(ix)/24/3600,'k--')
		topStr = sprintf('%.3f %cC',mlTopTemp(ix),char(176));
		tMT = text(diamLim(2)*0.95,mlTopTime(ix)/24/3600,topStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
		plot(diamLim, [1 1]*mlBotTime(ix)/24/3600,'k--')
		botStr = sprintf('%.3f %cC',mlBotTemp(ix),char(176));
		tMB = text(diamLim(2)*0.95,mlBotTime(ix)/24/3600,botStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
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
		
		dummyX = zeros(size(time_fl));
		
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
		time_fl_s = time_fl/24/3600;
		
		if(tempCsprl(1) < tempCsprl(end))
			
			plot(dummyX,time_fl_s,'Color','w');
			yLblFinal = flip(yLbl(~isnan(index)));
			ax.YTick = flip(time_fl_s(indexFinal));
			set(gca,'YDir','reverse');
			ax.YTickLabel = yLblFinal;
		else
			plot(dummyX,time_fl_s,'Color','w');
			yLblFinal = yLbl(~isnan(index));
			ax.YTick = time_fl_s(indexFinal);
			ax.YTickLabel = yLblFinal;
		end
		
		ylabel(sprintf('T (%cC)', char(176)));
		
		title([flight ' - Spiral ' num2str(ix) ' - CIP (Extended) - ' num2str(avgTime) 's Avg']);

		set(findall(gcf,'-property','FontSize'),'FontSize',26)
		set(tMB,'FontSize',14);
		set(tMT,'FontSize',14);
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/CIP-ND-Temp_' num2str(avgTime) 's/' flight '_CIP_ND-Temp_' num2str(avgTime) 's_S' num2str(ix)],Ftype,Fres)
		end
    end
end

if plotMDtemp
    for ix = loopVctr
		
		mass_twc = cipMass_hybrid_igf.(sprlNames{ix});
		
		tempCsprl = sDistF.tempC_orig.(sprlNames{ix});
		time_fl = sDistF.time_secsFL_orig.(sprlNames{ix});
		
		if avgTime == 1
			time_secs = sDistF.time_secs_orig.(sprlNames{ix});	
		else
			time_secs = sDistF.time_secs_avg.(sprlNames{ix});	
		end
		

        if saveFigs && noDisp
            figure('visible','off','Position', [10,10,1500,1000]);
        else
            figure('Position', [10,10,1500,1000]);
		end
		
		set(gcf,'defaultAxesColorOrder',[0 0 0; 0 0 0]);

        yyaxis left
		ax = gca;
        contourf(cipExt_binMid*10,time_secs/24/3600,log10(mass_twc),MDLogLim(1):0.1:MDLogLim(2),'LineColor','none');
        xlabel('D [mm]');
        ylabel('Time');
		set(ax,'XMinorTick','on');
		colormap(jetmod);
        c=colorbar;
		set(c,'Location','southoutside');
		ylabel(c,'log_{10}Mass_{TWC}M(D) [g cm^{-4}]');
		set(ax, 'CLim', MDLogLim);
		set(ax,'XLim',diamLim);
		
		% Plot ML top/bottom locations and annotate with temp
		hold on
		plot(diamLim, [1 1]*mlTopTime(ix)/24/3600,'k--')
		topStr = sprintf('%.3f %cC',mlTopTemp(ix),char(176));
		tMT = text(diamLim(2)*0.95,mlTopTime(ix)/24/3600,topStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
		plot(diamLim, [1 1]*mlBotTime(ix)/24/3600,'k--')
		botStr = sprintf('%.3f %cC',mlBotTemp(ix),char(176));
		tMB = text(diamLim(2)*0.95,mlBotTime(ix)/24/3600,botStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
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
		
		dummyX = zeros(size(time_fl));
		
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
		time_fl_s = time_fl/24/3600;
		
		if(tempCsprl(1) < tempCsprl(end))
			
			plot(dummyX,time_fl_s,'Color','w');
			yLblFinal = flip(yLbl(~isnan(index)));
			ax.YTick = flip(time_fl_s(indexFinal));
			set(gca,'YDir','reverse');
			ax.YTickLabel = yLblFinal;
		else
			plot(dummyX,time_fl_s,'Color','w');
			yLblFinal = yLbl(~isnan(index));
			ax.YTick = time_fl_s(indexFinal);
			ax.YTickLabel = yLblFinal;
		end
		
		ylabel(sprintf('T (%cC)', char(176)));
		

		title([flight ' - Spiral ' num2str(ix) ' - CIP - ' num2str(avgTime) 's Avg']);

		set(findall(gcf,'-property','FontSize'),'FontSize',26)
		set(tMB,'FontSize',14);
		set(tMT,'FontSize',14);
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/CIP-MD-Temp_' num2str(avgTime) 's/' flight '_CIP_MD-Temp_' num2str(avgTime) 's_S' num2str(ix)],Ftype,Fres)
		end
    end
end

if plotTWCextndRatio
    for ix = loopVctr
		
		cipMass = cipMass_hybrid_igf.(sprlNames{ix});
		
		cipTWC_obsOnly = nansum(cipMass(:,1:length(cip_binMin)),2);
		cipTWC_extOnly = nansum(cipMass(:,length(cip_binMin)+1:end),2);
		
		twcRatio = cipTWC_extOnly./cipTWC_obsOnly;

		ratioExcd = find(twcRatio > 0.5);
		
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		
		nanRatio = find(isnan(twcRatio));
		plotNanRatio = NaN(size(twcRatio));
		plotNanRatio(nanRatio) = 0;
		
		
		plot(twcRatio,'b.','MarkerSize',25);
		hold on
		plot(plotNanRatio,'r*','MarkerSize',5);

		title([flight ' - Spiral ' num2str(ix) ' - CIP - ' num2str(avgTime) 's Avg']);
		
		yl = ylim;
		ylL = 0 - ((yl(2)-yl(1))/20);
		ylim([ylL yl(2)])

		xlabel(['Time Dimension [each pt = ' num2str(avgTime) ' sec]']);
		ylabel('$$\frac{TWC _{extndOnly}}{TWC _{obsOnly}}$$','Interpreter','latex');
% 		set(gca,'Yscale','log');
		set(gca,'XMinorTick','on','YMinorTick','on');
		set(findall(gcf,'-property','FontSize'),'FontSize',28)
		grid


        if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/CIP-TWCextndRatio_' num2str(avgTime) 's/' flight '_CIP-TWCextndRatio_' num2str(avgTime) 's_S' num2str(ix)],Ftype,Fres)
		end
		
		% Plot individual M(D) for all points with extended portion of TWC exceeding 50% of observed TWC
% 		for irx = 1:length(ratioExcd)
% 			if saveFigs && noDisp
% 				figure('visible','off','Position', [10,10,1200,700]);
% 			else
% 				figure('Position', [10,10,1200,700]);
% 			end
% 			
% 			tempMassObs = (cipMass(ratioExcd(irx),1:length(cip_binMin)).*cipExt_binwidth(1:length(cip_binMin))').*1e6; % Convert to g m-3
% 			tempTWCObs = nansum(tempMassObs);
% 			tempMassExt = (cipMass(ratioExcd(irx),length(cip_binMin)+1:end).*cipExt_binwidth(length(cip_binMin)+1:end)').*1e6; % Convert to g m-3
% 			tempTWCExt = nansum(tempMassExt);
% 			
% 			
% 			stairs(cipExt_binMin, cipMass(ratioExcd(irx),:)','LineWidth', 2);
% 			
% 			title({sprintf('%s - Spiral %d - CIP - SDs where TWC_{extnd} > 50%% of TWC_{obs} - #%d',flight,ix,ratioExcd(irx)),...
% 				sprintf('TWC_{obs} = %.5f g m^{-3}     TWC_{extnd} = %.5f g m^{-3}',tempTWCObs,tempTWCExt)});
% 			
% 			xlabel('D [cm]');
% 			ylabel('M_{twc}(D) [g cm^{-4}]');
% 			set(gca,'Yscale','log','XScale','log');
% 			if ~isempty(diamLim)
% 				xlim(diamLim);
% 			end
% 			if ~isempty(MDLim)
% 				ylim(MDLim);
% 			end
% 			set(gca,'XMinorTick','on','YMinorTick','on');
% 			set(findall(gcf,'-property','FontSize'),'FontSize',28)
% 			grid;
% 			
% 			if saveFigs
% 				set(gcf,'Units','Inches');
% 				pos = get(gcf,'Position');
% 				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% 				print([saveDir '/CIP-TWCextndRatio_' num2str(avgTime) 's/' flight '_CIP-MD-TWCratioExcd_S' num2str(ix) '_' num2str(ratioExcd(irx)) '_5cm'],Ftype,Fres)
% 			end
% 		end
    end
end

if plotNDtempBinned
	if getLims
		maxConc = [];
		minConc = [];
	end
	for ix = loopVctr
		
		cipConcSprl = cipConc_hybrid_igf.(sprlNames{ix});
		
		tempCsprl = sDistF.tempC_avg.(sprlNames{ix});
		
		tempCrnd = NaN(length(tempCsprl),1);
		for ii=1:length(tempCsprl)
			if tempCsprl(ii) < 0
				tempCrnd(ii) = ceil(tempCsprl(ii));
			elseif tempCsprl(ii) > 0
				tempCrnd(ii) = floor(tempCsprl(ii));
			end
		end
		
		tempBins = min(tempCrnd):max(tempCrnd);
		
		for iii=1:length(tempBins)
			cipConc = cipConcSprl(tempCrnd == tempBins(iii),:); 
			if all(isnan(cipConc))
				continue
			end
			
			if getLims
				tmpMean = nanmean(cipConc(:,1:45),1);%1:45 gives us 0-4mm bins
				maxConc = [maxConc; nanmax(tmpMean(tmpMean>0))];
				minConc = [minConc; nanmin(tmpMean(tmpMean>0))];
				% fprintf('\n%dC\tminConc = %.2e\n',tempBins(iii),nanmin(tmpMean(tmpMean>0)))
				% fprintf('%dC\tmaxConc = %.2g\n',tempBins(iii),nanmax(tmpMean(tmpMean>0)))
			else
				if saveFigs && noDisp
					figure('visible','off','Position', [10,10,1500,1000]);
				else
					figure('Position', [10,10,1500,1000]);
				end
				
				stairs(cipExt_binMin(1:length(cip_binMin))*10, nanmean(cipConc(:,1:length(cip_binMin)),1)', 'Color',[0 0 0.3], 'LineWidth', 2);
				hold on
				stairs(cipExt_binMin(length(cip_binMin):end)*10, nanmean(cipConc(:,length(cip_binMin):end),1)','b', 'LineWidth', 2);
				
				title(sprintf('%s - Spiral %d - CIP - %d%cC avg',flight,ix,tempBins(iii),char(176)));
				
				xlabel('D [mm]');
				ylabel('N(D) [cm^{-4}]');
				set(gca,'Yscale','log','XScale','log');
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
					print(sprintf('%s/CIP-ND-TempBinned/%s_CIP_ND-Temp_S%d_%ddegC',saveDir,flight,ix,tempBins(iii)),Ftype,Fres)
				end
			end
		end
	end
	if getLims
		fprintf('minConc = %.2e\n',nanmin(minConc))
		fprintf('mean_minConc = %.2e\n',nanmean(minConc))
		fprintf('maxConc = %.2g\n',nanmax(maxConc))
	end
end

if plotMDtempBinned
	if getLims
		maxMass = [];
		minMass = [];
	end
	for ix = loopVctr
		
		cipMassSprl = cipMass_hybrid_igf.(sprlNames{ix});
		
		tempCsprl = sDistF.tempC_avg.(sprlNames{ix});
		
		tempCrnd = NaN(length(tempCsprl),1);
		for ii=1:length(tempCsprl)
			if tempCsprl(ii) < 0
				tempCrnd(ii) = ceil(tempCsprl(ii));
			elseif tempCsprl(ii) > 0
				tempCrnd(ii) = floor(tempCsprl(ii));
			end
		end
		
		tempBins = min(tempCrnd):max(tempCrnd);
		
		for iii=1:length(tempBins)
			cipMass = cipMassSprl(tempCrnd == tempBins(iii),:); 
			if all(isnan(cipMass))
				continue
			end
			
			if getLims
				tmpMean = nanmean(cipMass(:,1:45),1);%1:45 gives us 0-4mm bins
				maxMass = [maxMass; nanmax(tmpMean(tmpMean>0))];
				minMass = [minMass; nanmin(tmpMean(tmpMean>0))];
				% fprintf('\n%dC\tminMass = %.2e\n',tempBins(iii),nanmin(tmpMean(tmpMean>0)))
				% fprintf('%dC\tmaxMass = %.2g\n',tempBins(iii),nanmax(tmpMean(tmpMean>0)))
			else
				if saveFigs && noDisp
					figure('visible','off','Position', [10,10,1500,1000]);
				else
					figure('Position', [10,10,1500,1000]);
				end
				
				meanMass = nanmean(cipMass,1);
				
				yyaxis left
				stairs(cipExt_binMin(1:length(cip_binMin))*10, meanMass(1:length(cip_binMin))', 'Color',[0 0 0.3], 'LineWidth', 2);
				hold on
				stairs(cipExt_binMin(length(cip_binMin):end)*10, meanMass(length(cip_binMin):end)', 'b-', 'LineWidth', 2);

				title(sprintf('%s - Spiral %d - CIP - %d%cC avg',flight,ix,tempBins(iii),char(176)));

				xlabel('D [mm]');
				ylabel('M(D) [g cm^{-4}]');
				set(gca,'Yscale','log');
				set(gca,'Xscale','log');
				if ~isempty(diamLim)
					xlim(diamLim);
				end
				if ~isempty(MDtempBinLim)
					ylim(MDtempBinLim);
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
				plot(cipExt_binMin(1:length(cip_binMin))*10,massCDF(1:length(cip_binMin))*100,'Color',[0.3 0 0],'LineWidth',2);
				hold on
				plot(cipExt_binMin(length(cip_binMin):end)*10,massCDF(length(cip_binMin):end)*100,'r-','LineWidth',2);
				ylabel('M(D) CDF [%]');
				
				
				
				if saveFigs
					set(gcf,'Units','Inches');
					pos = get(gcf,'Position');
					set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
					print(sprintf('%s/CIP-MD-TempBinned/%s_CIP_MD-Temp_S%d_%ddegC',saveDir,flight,ix,tempBins(iii)),Ftype,Fres)
				end
			end
		end
	end
	if getLims
		fprintf('minMass = %.2e\n',nanmin(minMass))
		fprintf('mean_minMass = %.2e\n',nanmean(minMass))
		fprintf('maxMass = %.2g\n',nanmax(maxMass))
	end
end