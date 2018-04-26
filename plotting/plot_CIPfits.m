clearvars; close all;

%% Specify various plotting/calculation parameters
flights = {'20150617','20150620','20150701','20150702','20150706','20150709'};
% flights = {'20150706'};

avgTime = 10;

zoom = 0;

contLevs = 64;

fileIdStr = ['_Fit-CIP_' num2str(avgTime) 'secAvg_1.2cm'];
% fileIdStr = ['_Fit-CIP_' num2str(avgTime) 'secAvg_5cm'];
titleIdStr = '';
if zoom
	fNameAppnd = '_zoom';
else
	fNameAppnd = '';
end

% Standard plots, but using hybrid (CIP obs + CIP extended) SDs
getLims = 0; % If true, run without plotting and print values to assist in setting axis limits

plotND				= 1;
plotNDall			= 1;
plotNDtemp			= 1;
plotNDtempBinned	= 1;

plotMD				= 0;
plotMDall			= 0;
plotMDtemp			= 0;
plotMDtempBinned	= 0;

plotARatTemp        = 0;

plotTWCextndRatio	= 0;
plotMDratioExcd		= 0; % Plot indvdl M(D) for periods where mass ratio between obs and extended is exceeded
plotNDratioExcd		= 0; % Plot indvdl N(D) for periods where mass ratio between obs and extended is exceeded


% Plots using all spirals over all flights
plotARsprls  = 0;
plotTWCsprls = 0;


if zoom
	diamLim = [0.1 2.5];
else
	diamLim = [0.1 6];
end
% diamLim = [0.1 inf];


saveFigs	= 1;
noDisp		= 1;
fullVec		= 0;
Ftype		= '-dpdf';
% Ftype		= '-dpng';

if fullVec && strcmp(Ftype,'-dpdf')
	Fres = '-painters'; % Use this for fully vectorized files
else
	Fres = '-r0'; % Use this for smaller files - saves figure with same resolution/size as displayed on screen
end


savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';
dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';
% load('/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/Code/PECAN-MP-Analysis/utilities/HomeyerRainbowCmap.mat');


% Initialize variables to be analyzed over all flights/spirals
cipAR_sprlAvgs = NaN(42,34); % 42 spirals, with 34 diameter bins each

iSprl = 1; % Counter for total number of spirals (used in whole-project variable concat)


%% Create figure handles for any multi-flight figs
if plotTWCsprls
	twcSfig = figure('Position', [10,10,1500,1000]);
	twcColors = varycolor(42);
end

for iFlt = 1:length(flights)
	flight = flights{iFlt};
	fprintf('\nPlotting %s...\n',flight);
	%% Load in struct of all original SD data and the CIP Fitted Data
	if avgTime == 1
		sDistFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.CIP.10secAvg.mat']; % 1-sec averages are saved in all averaged output files
	else
		sDistFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.CIP.' num2str(avgTime) 'secAvg.mat'];
	end

	sDistF = load(sDistFile);

	load([dataPath 'mp-data/' flight '/sDist/' flight fileIdStr '.mat']);

	%% Create directories to save plots in if they don't already exist
	if saveFigs
		%%% Save info and directory operations for CIP Fit plots
		saveDir = [savePath flight '/CIP-Ext'];
		if (exist(saveDir, 'dir') ~= 7)
			mkdir(saveDir)
		end

		if (plotND && exist([saveDir '/CIP-ND-avg_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-ND-avg_' num2str(avgTime) 's'])
		end
		if (plotNDall && exist([saveDir '/CIP-ND_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-ND_' num2str(avgTime) 's'])
		end
		if (plotMD && exist([saveDir '/CIP-MD-avg_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-MD-avg_' num2str(avgTime) 's'])
		end
		if (plotMDall && exist([saveDir '/CIP-MD_' num2str(avgTime) 's'], 'dir') ~= 7)
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
		if (plotMDratioExcd && exist([savePath '/CIP-MD-twcRatioExcd'], 'dir') ~= 7)
			mkdir([savePath '/CIP-MD-TWCratioExcd'])
		end
		if (plotNDratioExcd && exist([savePath '/CIP-ND-twcRatioExcd'], 'dir') ~= 7)
			mkdir([savePath '/CIP-ND-TWCratioExcd'])
		end
		if (plotNDtempBinned && exist([saveDir '/CIP-ND-TempBinned'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-ND-TempBinned'])
		end
		if (plotMDtempBinned && exist([saveDir '/CIP-MD-TempBinned'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-MD-TempBinned'])
		end
		if (plotARatTemp && exist([saveDir '/CIP-ARat-Temp_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-ARat-Temp_' num2str(avgTime) 's'])
		end
		
		%%% Save info and directory operations for multi-flight plots
		saveDirObs = savePath;
		if (exist(saveDirObs, 'dir') ~= 7)
			mkdir(saveDirObs)
		end
		
		if (plotARsprls && exist([saveDirObs '/CIP-AR-Sprls_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDirObs '/CIP-AR-Sprls_' num2str(avgTime) 's'])
		end
		if (plotTWCsprls && exist([saveDirObs '/CIP-TWC-Sprls_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([savePath '/CIP-TWC-Sprls_' num2str(avgTime) 's'])
		end
	end

	%%% Save a list of variables we want to keep across all iterations
	if iFlt == 1
		initialVars = who;
		initialVars{end+1} = 'initialVars';
	end
	
	loopVctr = 1:length(sprlNames);
% 	loopVctr = [7];

	%% Specify flight-specific plotting parameters
	tempRangeAll = [-18.5 22];

	if zoom
		NDLim = [1e-4 10];
		MDLim = [1e-8 5e-5];
	else
		NDLim = [1e-6 20];
		MDLim = [1e-10 5e-5];
	end
	NDLogLim = [-4 2];
	MDLogLim = [-8 -4];
	ARlim = [0.25 0.85];

	
	% Switch statement to allow for flight-specific plot limits
	%{
	switch flight
		case '20150617'
			if zoom
				NDLim = [1e-4 10];
				MDLim = [1e-8 5e-5]; 
			else
				NDLim = [1e-6 20];
				MDLim = [1e-10 5e-5]; 
			end
			NDLogLim = [-4 2];
			MDLogLim = [-8 -4];
			NtLim = [5.5e-6 0.5];
			TWClim = [4.4e-5 2.1];
			DmmLim = [0 2.5];
		case '20150620'
			if zoom
				NDLim = [1e-4 10];
				MDLim = [1e-8 5e-5]; 
			else
				NDLim = [1e-6 20];
				MDLim = [1e-10 5e-5]; 
			end
			NDLogLim = [-4 2];
			MDLogLim = [-8 -4];
			NtLim = [5.5e-6 0.5];
			TWClim = [4.4e-5 2.1];
			DmmLim = [0 2.5];
		case '20150701'
			if zoom
				NDLim = [1e-4 10];
				MDLim = [1e-8 5e-5]; 
			else
				NDLim = [1e-6 20];
				MDLim = [1e-10 5e-5]; 
			end
			NDLogLim = [-4 2];
			MDLogLim = [-8 -4];
			NtLim = [5.5e-6 0.5];
			TWClim = [4.4e-5 2.1];
			DmmLim = [0 2.5];
		case '20150702'
			if zoom
				NDLim = [1e-4 10];
				MDLim = [1e-8 5e-5]; 
			else
				NDLim = [1e-6 20];
				MDLim = [1e-10 5e-5]; 
			end
			NDLogLim = [-4 2];
			MDLogLim = [-8 -4];
			NtLim = [5.5e-6 0.5];
			TWClim = [4.4e-5 2.1];
			DmmLim = [0 2.5];
		case '20150706'
			if zoom
				NDLim = [1e-4 10];
				MDLim = [1e-8 5e-5]; 
			else
				NDLim = [1e-6 20];
				MDLim = [1e-10 5e-5]; 
			end
			NDLogLim = [-4 2];
			MDLogLim = [-8 -4];
			NtLim = [5.5e-6 0.5];
			TWClim = [4.4e-5 2.1];
			DmmLim = [0 2.5];
		case '20150709'
			if zoom
				NDLim = [1e-4 10];
				MDLim = [1e-8 5e-5]; 
			else
				NDLim = [1e-6 20];
				MDLim = [1e-10 5e-5]; 
			end
			NDLogLim = [-4 2];
			MDLogLim = [-8 -4];
			NtLim = [5.5e-6 0.5];
			TWClim = [4.4e-5 2.1];
			DmmLim = [0 2.5];
	end
	%}

	%% Standard plots
	if plotND
		for ix = loopVctr

			cipConc = nanmean(cipConc_hybrid_igf.(sprlNames{ix}),1);

			if saveFigs && noDisp
				figure('visible','off','Position', [10,10,1200,700]);
			else
				figure('Position', [10,10,1200,700]);
			end

			stairs(cipExt_binMin(1:length(cip_binMin))*10, cipConc(1:length(cip_binMin))', 'Color',[0 0 0.3], 'LineWidth', 2);
			hold on
			stairs(cipExt_binMin(length(cip_binMin):end)*10, cipConc(length(cip_binMin):end)','b', 'LineWidth', 2);

			title([flight ' - Spiral ' num2str(ix) ' - CIP (Extended) (using ' num2str(avgTime) 's Avg PSDs)']);

			xlabel('D (mm)');
			ylabel('N(D) (cm^{-4})');
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
				tightfig(gcf);
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				print([saveDir '/CIP-ND-avg_' num2str(avgTime) 's/' flight '_CIP_ND_' num2str(avgTime) 's_S' num2str(ix) fNameAppnd],Ftype,Fres)
			end
		end
	end
	if plotNDall
		if getLims
			maxConc = [];
			minConc = [];
		end
		for ix = loopVctr

			cipConc = cipConc_hybrid_igf.(sprlNames{ix});

			colors = varycolor(size(cipConc,1));
			if getLims
				if zoom
					tmpMean = nanmean(cipConc(:,1:37),1);%1:37 gives us 0-2.6mm bins
					maxConc = [maxConc; nanmax(tmpMean(tmpMean>0))];
					minConc = [minConc; nanmin(tmpMean(tmpMean>0))];
				else
					tmpMean = nanmean(cipConc(:,1:53),1);%1:53 gives us 0-5.8mm bins
					maxConc = [maxConc; nanmax(tmpMean(tmpMean>0))];
					minConc = [minConc; nanmin(tmpMean(tmpMean>0))];
				end
			else
				if saveFigs && noDisp
					figure('visible','off','Position', [10,10,1200,700]);
				else
					figure('Position', [10,10,1200,700]);
				end

				hold on

				for icx=1:size(cipConc,1)
					stairs(cipExt_binMin(1:length(cip_binMin))*10, cipConc(icx,1:length(cip_binMin))','Color',colors(icx,:), 'LineWidth', 2,'DisplayName',num2str(icx));
					stairs(cipExt_binMin(length(cip_binMin):end)*10, cipConc(icx,length(cip_binMin):end)','Color',colors(icx,:), 'LineWidth', 0.5,'DisplayName',[num2str(icx) 'ext']);
				end
				title([flight ' - Spiral ' num2str(ix) ' - CIP (Extended) (using ' num2str(avgTime) 's Avg PSDs)']);

				xlabel('D (mm)');
				ylabel('N(D) (cm^{-4})');
				set(gca,'Yscale','log');
				set(gca,'Xscale','log');
				if ~isempty(diamLim)
					xlim(diamLim);
				end
				if ~isempty(NDLim)
					ylim(NDLim);
				end
				
				colormap(colors)
				cbar = colorbar('Ticks',[0 1],'TickLabels',{'FirstPSD','LastPSD'});
				
				set(gca,'XMinorTick','on','YMinorTick','on');
				set(findall(gcf,'-property','FontSize'),'FontSize',28)
				set(cbar,'FontSize',14);
				grid


				if saveFigs
					set(gcf,'Units','Inches');
					pos = get(gcf,'Position');
					set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
					print([saveDir '/CIP-ND_' num2str(avgTime) 's/' flight '_CIP_ND_' num2str(avgTime) 's_S' num2str(ix) '_all' fNameAppnd],Ftype,Fres)
				end
			end
		end
		if getLims
			fprintf('minConc = %.2e\n',nanmin(minConc))
			fprintf('mean_minConc = %.2e\n',nanmean(minConc))
			fprintf('maxConc = %.2g\n',nanmax(maxConc))
		end
	end
	if plotNDtemp
		if any(isinf(diamLim)) || isempty(diamLim)
			diamLim = [min(cipExt_binMid*10) max(cipExt_binMid*10)];
		end
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
			contourf(cipExt_binMid*10,time_secs/24/3600,log10(cipConc),linspace(NDLogLim(1),NDLogLim(2),contLevs),'LineColor','none');
			xlabel('D (mm)');
			ylabel('Time');
			set(ax,'XMinorTick','on');
			colormap(HomeyerRainbow(contLevs));
			c=colorbar;
			set(c,'Location','southoutside');
			ylabel(c,'log_{10}N(D) (cm^{-4})');
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
				set(gca,'YLim',[time_fl_s(1) time_fl_s(end)]);
				set(gca,'YDir','reverse');
				ax.YTickLabel = yLblFinal;
			else
				plot(dummyX,time_fl_s,'Color','w');
				yLblFinal = yLbl(~isnan(index));
				ax.YTick = time_fl_s(indexFinal);
				set(gca,'YLim',[time_fl_s(1) time_fl_s(end)]);
				ax.YTickLabel = yLblFinal;
			end


			ylabel(sprintf('T (%cC)', char(176)));

			title([flight ' - Spiral ' num2str(ix) ' - CIP (Extended) (using ' num2str(avgTime) 's Avg PSDs)']);

			set(findall(gcf,'-property','FontSize'),'FontSize',26)
			set(tMB,'FontSize',14);
			set(tMT,'FontSize',14);

			if saveFigs
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				print([saveDir '/CIP-ND-Temp_' num2str(avgTime) 's/' flight '_CIP_ND-Temp_' num2str(avgTime) 's_S' num2str(ix) fNameAppnd],'-dpng',Fres)
			end
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
					if zoom
						tmpMean = nanmean(cipConc(:,1:37),1);%1:37 gives us 0-2.6mm bins
						maxConc = [maxConc; nanmax(tmpMean(tmpMean>0))];
						minConc = [minConc; nanmin(tmpMean(tmpMean>0))];
					else
						tmpMean = nanmean(cipConc(:,1:53),1);%1:53 gives us 0-5.8mm bins
						maxConc = [maxConc; nanmax(tmpMean(tmpMean>0))];
						minConc = [minConc; nanmin(tmpMean(tmpMean>0))];
					end
				else
					if saveFigs && noDisp
						figure('visible','off','Position', [10,10,1500,1000]);
					else
						figure('Position', [10,10,1500,1000]);
					end

					stairs(cipExt_binMin(1:length(cip_binMin))*10, nanmean(cipConc(:,1:length(cip_binMin)),1)', 'Color',[0 0 0.3], 'LineWidth', 2);
					hold on
					stairs(cipExt_binMin(length(cip_binMin):end)*10, nanmean(cipConc(:,length(cip_binMin):end),1)','b', 'LineWidth', 2);

					title(sprintf('%s - Spiral %d - CIP (Extended) - %d%cC avg',flight,ix,tempBins(iii),char(176)));

					xlabel('D (mm)');
					ylabel('N(D) (cm^{-4})');
					set(gca,'Yscale','log','XScale','log');
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
						tightfig(gcf);
						set(gcf,'Units','Inches');
						pos = get(gcf,'Position');
						set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
						print(sprintf('%s/CIP-ND-TempBinned/%s_CIP_ND-TempB_S%d_%d_%ddegC%s',saveDir,flight,ix,iii,tempBins(iii),fNameAppnd),Ftype,Fres)
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
	
	if plotMD
		for ix = loopVctr

			meanMass = nanmean(cipMass_hybrid_igf.(sprlNames{ix}),1);

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

			title([flight ' - Spiral ' num2str(ix) ' - CIP (Extended) (using ' num2str(avgTime) 's Avg PSDs)']);

			xlabel('D (mm)');
			ylabel('M(D) (g cm^{-4})');
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
			ylabel('M(D) CDF (%)');


			if saveFigs
				tightfig(gcf);
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				print([saveDir '/CIP-MD-avg_' num2str(avgTime) 's/' flight '_CIP_MD_' num2str(avgTime) 's_S' num2str(ix) fNameAppnd],Ftype,Fres)
			end
		end
	end
	if plotMDall
		if getLims
			maxMass = [];
			minMass = [];
		end
		for ix = loopVctr
			cipMass = cipMass_hybrid_igf.(sprlNames{ix});

			colors = varycolor(size(cipMass,1));
			if getLims
				if zoom
					tmpMean = nanmean(cipMass(:,1:37),1);%1:37 gives us 0-2.6mm bins
					maxMass = [maxMass; nanmax(tmpMean(tmpMean>0))];
					minMass = [minMass; nanmin(tmpMean(tmpMean>0))];
				else
					tmpMean = nanmean(cipMass(:,1:53),1);%1:53 gives us 0-5.8mm bins
					maxMass = [maxMass; nanmax(tmpMean(tmpMean>0))];
					minMass = [minMass; nanmin(tmpMean(tmpMean>0))];
				end
			else
				if saveFigs && noDisp
					figure('visible','off','Position', [10,10,1200,700]);
				else
					figure('Position', [10,10,1200,700]);
				end

				hold on

				for icx=1:size(cipMass,1)
					stairs(cipExt_binMin(1:length(cip_binMin))*10, cipMass(icx,1:length(cip_binMin))','Color',colors(icx,:), 'LineWidth', 2,'DisplayName',num2str(icx));
					stairs(cipExt_binMin(length(cip_binMin):end)*10, cipMass(icx,length(cip_binMin):end)','Color',colors(icx,:), 'LineWidth', 0.5,'DisplayName',[num2str(icx) 'ext']);
				end
				title([flight ' - Spiral ' num2str(ix) ' - CIP (Extended) (using ' num2str(avgTime) 's Avg PSDs)']);

				xlabel('D (mm)');
				ylabel('M(D) (g cm^{-4})');
				set(gca,'Yscale','log');
				set(gca,'Xscale','log');
				if ~isempty(diamLim)
					xlim(diamLim);
				end
				if ~isempty(MDLim)
					ylim(MDLim);
				end
				
				colormap(colors)
				cbar = colorbar('Ticks',[0 1],'TickLabels',{'FirstPSD','LastPSD'});
				
				set(gca,'XMinorTick','on','YMinorTick','on');
				set(findall(gcf,'-property','FontSize'),'FontSize',28)
				set(cbar,'FontSize',14);
				grid


				if saveFigs
					set(gcf,'Units','Inches');
					pos = get(gcf,'Position');
					set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
					print([saveDir '/CIP-MD_' num2str(avgTime) 's/' flight '_CIP_MD_' num2str(avgTime) 's_S' num2str(ix) '_all' fNameAppnd],Ftype,Fres)
				end
			end
		end
		if getLims
			fprintf('minMass = %.2e\n',nanmin(minMass))
			fprintf('mean_minMass = %.2e\n',nanmean(minMass))
			fprintf('maxMass = %.2g\n',nanmax(maxMass))
		end
	end
	if plotMDtemp
		if any(isinf(diamLim)) || isempty(diamLim)
			diamLim = [min(cipExt_binMid*10) max(cipExt_binMid*10)];
		end
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
			contourf(cipExt_binMid*10,time_secs/24/3600,log10(mass_twc),linspace(MDLogLim(1),MDLogLim(2),contLevs),'LineColor','none');
			xlabel('D (mm)');
			ylabel('Time');
			set(ax,'XMinorTick','on');
			colormap(HomeyerRainbow(contLevs));
			c=colorbar;
			set(c,'Location','southoutside');
			ylabel(c,'log_{10}M(D) (g cm^{-4})');
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
				set(gca,'YLim',[time_fl_s(1) time_fl_s(end)]);
				set(gca,'YDir','reverse');
				ax.YTickLabel = yLblFinal;
			else
				plot(dummyX,time_fl_s,'Color','w');
				yLblFinal = yLbl(~isnan(index));
				ax.YTick = time_fl_s(indexFinal);
				set(gca,'YLim',[time_fl_s(1) time_fl_s(end)]);
				ax.YTickLabel = yLblFinal;
			end

			ylabel(sprintf('T (%cC)', char(176)));


			title([flight ' - Spiral ' num2str(ix) ' - CIP (Extended) (using ' num2str(avgTime) 's Avg PSDs)']);

			set(findall(gcf,'-property','FontSize'),'FontSize',26)
			set(tMB,'FontSize',14);
			set(tMT,'FontSize',14);

			if saveFigs
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				print([saveDir '/CIP-MD-Temp_' num2str(avgTime) 's/' flight '_CIP_MD-Temp_' num2str(avgTime) 's_S' num2str(ix) fNameAppnd],'-dpng',Fres)
			end
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
					if zoom
						tmpMean = nanmean(cipMass(:,1:37),1);%1:37 gives us 0-2.6mm bins
						maxMass = [maxMass; nanmax(tmpMean(tmpMean>0))];
						minMass = [minMass; nanmin(tmpMean(tmpMean>0))];
					else
						tmpMean = nanmean(cipMass(:,1:53),1);%1:53 gives us 0-5.8mm bins
						maxMass = [maxMass; nanmax(tmpMean(tmpMean>0))];
						minMass = [minMass; nanmin(tmpMean(tmpMean>0))];
					end
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

					title(sprintf('%s - Spiral %d - CIP (Extended) - %d%cC avg',flight,ix,tempBins(iii),char(176)));

					xlabel('D (mm)');
					ylabel('M(D) (g cm^{-4})');
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
					ylabel('M(D) CDF (%)');



					if saveFigs
						tightfig(gcf);
						set(gcf,'Units','Inches');
						pos = get(gcf,'Position');
						set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
						print(sprintf('%s/CIP-MD-TempBinned/%s_CIP_MD-TempB_S%d_%d_%ddegC%s',saveDir,flight,ix,iii,tempBins(iii),fNameAppnd),Ftype,Fres)
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
	
	if plotARatTemp
		for ix = loopVctr
			diamLim = [0.1 2.1];
			
			cipARat = sDistF.mean_areaRatio_avg.(sprlNames{ix});

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
			contourf(sDistF.bin_mid,time_secs/24/3600,cipARat,linspace(ARlim(1),ARlim(2),contLevs),'LineColor','none');
			xlabel('D (mm)');
			ylabel('Time');
			set(ax,'XMinorTick','on');
			colormap(HomeyerRainbow(contLevs));
			c=colorbar;
			set(c,'Location','southoutside');
			ylabel(c,'log_{10}N(D) (cm^{-4})');
			set(ax, 'CLim', ARlim);
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
				set(gca,'YLim',[time_fl_s(1) time_fl_s(end)]);
				set(gca,'YDir','reverse');
				ax.YTickLabel = yLblFinal;
			else
				plot(dummyX,time_fl_s,'Color','w');
				yLblFinal = yLbl(~isnan(index));
				ax.YTick = time_fl_s(indexFinal);
				set(gca,'YLim',[time_fl_s(1) time_fl_s(end)]);
				ax.YTickLabel = yLblFinal;
			end


			ylabel(sprintf('T (%cC)', char(176)));

			title([flight ' - Spiral ' num2str(ix) ' - CIP (Extended) (using ' num2str(avgTime) 's Avg PSDs)']);

			set(findall(gcf,'-property','FontSize'),'FontSize',26)
			set(tMB,'FontSize',14);
			set(tMT,'FontSize',14);

			if saveFigs
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				print([saveDir '/CIP-ARat-Temp_' num2str(avgTime) 's/' flight '_CIP_ARat-Temp_' num2str(avgTime) 's_S' num2str(ix) fNameAppnd],'-dpng',Fres)
			end
		end
	end
	
	
	if plotTWCextndRatio
		for ix = loopVctr

			cipMass = cipMass_hybrid_igf.(sprlNames{ix});
			twcRatioSprl = twcRatio.(sprlNames{ix});
			twcRexcdSprl = twcRatioExcdIx.(sprlNames{ix});


			if saveFigs && noDisp
				figure('visible','off','Position', [10,10,1200,700]);
			else
				figure('Position', [10,10,1200,700]);
			end

			nanRatio = find(isnan(twcRatioSprl));
			nanPltRatio = NaN(size(twcRatioSprl));
			nanPltRatio(nanRatio) = 0;

			excdPltRatio = NaN(size(twcRatioSprl));
			excdPltRatio(twcRexcdSprl) = twcRatioSprl(twcRexcdSprl);

			plot(twcRatioSprl,'b.','MarkerSize',25);
			hold on
			plot(excdPltRatio,'rd','MarkerSize',15);
			plot(nanPltRatio,'k*','MarkerSize',5);

			title([flight ' - Spiral ' num2str(ix) '  - CIP (Extended) (using ' num2str(avgTime) 's Avg PSDs)']);

			yl = ylim;
			ylL = 0 - ((yl(2)-yl(1))/20);
			ylim([ylL yl(2)])

			xlabel(['Time Dimension (each pt = ' num2str(avgTime) ' sec)']);
			ylabel('$$\frac{TWC _{extndOnly}}{TWC _{obsOnly}}$$','Interpreter','latex');
	% 		set(gca,'Yscale','log');
			set(gca,'XMinorTick','on','YMinorTick','on');
			set(findall(gcf,'-property','FontSize'),'FontSize',28)
			grid
			legend('TWC Ratio','Ratio Thresh Exceeded','NaN','Location','best');


			if saveFigs
				tightfig(gcf);
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				print([saveDir '/CIP-TWCextndRatio_' num2str(avgTime) 's/' flight '_CIP_TWCextndRatio_' num2str(avgTime) 's_S' num2str(ix) fNameAppnd],Ftype,Fres)
			end
		end
	end

	if plotMDratioExcd
		% Plot individual M(D) for all points with extended portion of TWC exceeding threshold of observed TWC
		for ix = loopVctr
			cipMass = cipMass_hybrid_igf.(sprlNames{ix});
			cipMass_fitOnly = cipMass_ext_igf.(sprlNames{ix});
			tempCsprl = sDistF.tempC_orig.(sprlNames{ix});
			twcRatioSprl = twcRatio.(sprlNames{ix});
% 			twcRexcdSprl = find(twcRatioSprl >= 0.4 & twcRatioSprl <= 9.0);
			twcRexcdSprl = twcRatioExcdIx.(sprlNames{ix});
			
			for irx = 1:length(twcRexcdSprl)
				if saveFigs && noDisp
					figure('visible','off','Position', [10,10,1200,700]);
				else
					figure('Position', [10,10,1200,700]);
				end
				
				tempMassObs = (cipMass(twcRexcdSprl(irx),1:length(cip_binMin)).*cipExt_binwidth(1:length(cip_binMin))').*1e6; % Convert to g m-3
				tempTWCObs = nansum(tempMassObs);
				tempMassExt = (cipMass(twcRexcdSprl(irx),length(cip_binMin)+1:end).*cipExt_binwidth(length(cip_binMin)+1:end)').*1e6; % Convert to g m-3
				tempTWCExt = nansum(tempMassExt);
				
				
				% Determine average difference between the obs and fit
				tmpFitMass = cipMass_fitOnly(twcRexcdSprl(irx),1:length(cip_binMin));
				tmpMass = cipMass(twcRexcdSprl(irx),1:length(cip_binMin));
				tmpFitMass(tmpFitMass == 0) = NaN;
				tmpMass(tmpMass == 0) = NaN;
				massDiff = nanmean(tmpFitMass - tmpMass);
				
				
				stairs(cipExt_binMin(1:length(cip_binMin))*10, cipMass(twcRexcdSprl(irx),1:length(cip_binMin))','b-','LineWidth', 2);
				hold on
				stairs(cipExt_binMin(length(cip_binMin)+1:end)*10, cipMass(twcRexcdSprl(irx),length(cip_binMin)+1:end)','r-','LineWidth', 2);
				plot(cipExt_binMid*10,cipMass_fitOnly(twcRexcdSprl(irx),:)','k--','LineWidth', 2);
				
				t1 = sprintf('%s - Spiral %d - CIP - TWC_{extnd} = \\color{magenta}%.2f %% \\color{black}of TWC_{obs} - #%d',flight,ix,twcRatioSprl(twcRexcdSprl(irx))*100,twcRexcdSprl(irx));
				t2 = sprintf('TWC_{obs} = %.4f g m^{-3}  TWC_{extnd} = %.4f g m^{-3}  \\color{blue}\\mu = %.2f  \\lambda = %.2f',...
					tempTWCObs,tempTWCExt,cip_igf_nml.(sprlNames{ix})(twcRexcdSprl(irx),2),cip_igf_nml.(sprlNames{ix})(twcRexcdSprl(irx),3));
				% if abs(massDiff) > 1e-7
				%	t3 = sprintf('\\color{black}(M(D)_{ext}-M(D)_{obs})_{avg} = \\color{red}%.3g g cm^{-4}',massDiff);
				%else
					t3 = sprintf('\\color{black}(M(D)_{ext}-M(D)_{obs})_{avg} = %.3g g cm^{-4}',massDiff);
				%end
				
				title({t1,...
					t2,...
					t3});
				
				xlabel('D (mm)');
				ylabel('M(D) (g cm^{-4})');
				set(gca,'Yscale','log','XScale','log');
				if ~isempty(diamLim)
					xlim(diamLim);
				end
				if ~isempty(MDLim)
					ylim(MDLim);
				end
				set(gca,'XMinorTick','on','YMinorTick','on');
				set(findall(gcf,'-property','FontSize'),'FontSize',24)
				grid;
				
				fPre = sprintf('%.3f_',twcRatioSprl(twcRexcdSprl(irx)));
				
				if saveFigs
					tightfig(gcf);
					set(gcf,'Units','Inches');
					pos = get(gcf,'Position');
					set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
					print([savePath '/CIP-MD-TWCratioExcd/' fPre flight '_CIP_MD-TWCratioExcd_S' num2str(ix) '_' num2str(twcRexcdSprl(irx)) '.pdf'],Ftype,Fres)
				end
			end
		end
	end
	
	if plotNDratioExcd
		% Plot individual N(D) for all points with extended portion of TWC exceeding threshold of observed TWC
		for ix = loopVctr
			cipMass = cipMass_hybrid_igf.(sprlNames{ix});
			cipMass_fitOnly = cipMass_ext_igf.(sprlNames{ix});
			cipConc = cipConc_hybrid_igf.(sprlNames{ix});
			cipConc_fitOnly = cipConc_ext_igf.(sprlNames{ix});
			tempCsprl = sDistF.tempC_orig.(sprlNames{ix});
			twcRatioSprl = twcRatio.(sprlNames{ix});
% 			twcRexcdSprl = find(twcRatioSprl >= 0.4 & twcRatioSprl <= 9.0);
			twcRexcdSprl = twcRatioExcdIx.(sprlNames{ix});
			
			for irx = 1:length(twcRexcdSprl)
				if saveFigs && noDisp
					figure('visible','off','Position', [10,10,1200,700]);
				else
					figure('Position', [10,10,1200,700]);
				end
				
				tempMassObs = (cipMass(twcRexcdSprl(irx),1:length(cip_binMin)).*cipExt_binwidth(1:length(cip_binMin))').*1e6; % Convert to g m-3
				tempTWCObs = nansum(tempMassObs);
				tempMassExt = (cipMass(twcRexcdSprl(irx),length(cip_binMin)+1:end).*cipExt_binwidth(length(cip_binMin)+1:end)').*1e6; % Convert to g m-3
				tempTWCExt = nansum(tempMassExt);
				
				
				% Determine average difference between the obs and fit
				tmpFitConc = cipConc_fitOnly(twcRexcdSprl(irx),1:length(cip_binMin));
				tmpConc = cipConc(twcRexcdSprl(irx),1:length(cip_binMin));
				tmpFitConc(tmpFitConc == 0) = NaN;
				tmpConc(tmpConc == 0) = NaN;
				concDiff = nanmean(tmpFitConc - tmpConc);
				
				
				stairs(cipExt_binMin(1:length(cip_binMin))*10, cipConc(twcRexcdSprl(irx),1:length(cip_binMin))','b-','LineWidth', 2);
				hold on
				stairs(cipExt_binMin(length(cip_binMin)+1:end)*10, cipConc(twcRexcdSprl(irx),length(cip_binMin)+1:end)','r-','LineWidth', 2);
				plot(cipExt_binMid*10,cipConc_fitOnly(twcRexcdSprl(irx),:)','k--','LineWidth', 2);
				
				t1 = sprintf('%s - Spiral %d - CIP - TWC_{extnd} = \\color{magenta}%.2f %% \\color{black}of TWC_{obs} - #%d',flight,ix,twcRatioSprl(twcRexcdSprl(irx))*100,twcRexcdSprl(irx));
				t2 = sprintf('TWC_{obs} = %.4f g m^{-3}  TWC_{extnd} = %.4f g m^{-3}  \\color{blue}\\mu = %.2f  \\lambda = %.2f',...
					tempTWCObs,tempTWCExt,cip_igf_nml.(sprlNames{ix})(twcRexcdSprl(irx),2),cip_igf_nml.(sprlNames{ix})(twcRexcdSprl(irx),3));
				%if abs(concDiff) > 1e-7
				%	t3 = sprintf('\\color{black}(N(D)_{ext}-N(D)_{obs})_{avg} = \\color{red}%.3g cm^{-4}',concDiff);
				%else
					t3 = sprintf('\\color{black}(N(D)_{ext}-N(D)_{obs})_{avg} = %.3g cm^{-4}',concDiff);
				%end
				
				title({t1,...
					t2,...
					t3});
				
				xlabel('D (mm)');
				ylabel('N(D) (cm^{-4})');
				set(gca,'Yscale','log','XScale','log');
				if ~isempty(diamLim)
					xlim(diamLim);
				end
				if ~isempty(NDLim)
					ylim(NDLim);
				end
				set(gca,'XMinorTick','on','YMinorTick','on');
				set(findall(gcf,'-property','FontSize'),'FontSize',24)
				grid;
				
				fPre = sprintf('%.3f_',twcRatioSprl(twcRexcdSprl(irx)));
				
				if saveFigs
					tightfig(gcf);
					set(gcf,'Units','Inches');
					pos = get(gcf,'Position');
					set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
					print([savePath '/CIP-ND-TWCratioExcd/' fPre flight '_CIP_ND-TWCratioExcd_S' num2str(ix) '_' num2str(twcRexcdSprl(irx)) '.pdf'],Ftype,Fres)
				end
			end
		end
	end
	%% Concatenate/plot multi-flight variables
	for ix = loopVctr
		% Calculate area ratio averaged over all times (per bin) above ML only
		if ~isnan(mlTopTime(ix))
			[~, topIx] = min(abs(sDistF.time_secs_avg.(sprlNames{ix}) - mlTopTime(ix))); % ix of nearest time matching ML top
		else
			[~, topIx] = min(abs(sDistF.tempC_avg.(sprlNames{ix}))); % ix of nearest temp nearest 0 deg C
		end
		if sDistF.tempC_avg.(sprlNames{ix})(1) < sDistF.tempC_avg.(sprlNames{ix})(end) % Spiral down
			cipAR_sprlAvgs(iSprl,:) = nanmean(sDistF.mean_areaRatio_avg.(sprlNames{ix})(1:topIx-1,:),1);
		else % Spiral up
			cipAR_sprlAvgs(iSprl,:) = nanmean(sDistF.mean_areaRatio_avg.(sprlNames{ix})(topIx+1:end,:),1);
		end
		
		
		if plotTWCsprls
			TWCsprl = cipTWC_hybrid_igf.(sprlNames{ix});
			tempCsprl = sDistF.tempC_avg.(sprlNames{ix});
			
			set(0, 'CurrentFigure', twcSfig);
			hold on
			plot(TWCsprl,tempCsprl,'LineWidth',2,'Color',twcColors(iSprl,:));
			title('PECAN - CIP (Extended) - TWC - All Spirals');
			
			
			ylabel(sprintf('Temperature (%cC)', char(176)));
			xlabel('TWC (g m^{-3})')
			if ~isempty(tempRangeAll)
				ylim(tempRangeAll);
			end
			if ~isempty(TWClim)
				xlim(TWClim);
			end
			set(gca,'XMinorTick','on','YMinorTick','on','YDir','reverse','Xscale','log');
			set(findall(gcf,'-property','FontSize'),'FontSize',28)
			grid on
		end
		
		
		iSprl = iSprl+1;
	end
	% 	clearvars('-except',initialVars{:});
end

%% Save any multi-flight plots
if plotTWCsprls && saveFigs
	set(0, 'CurrentFigure', twcSfig);
	tightfig(gcf);
	set(gcf,'Units','Inches');
	pos = get(gcf,'Position');
	set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
	print([saveDir '/CIP-TWC-Sprls_' num2str(avgTime) 's/CIP_TWC-AllSprls_Temp' fNameAppnd],Ftype,Fres)
end


%% Make any whole-project plots
if plotARsprls
	colors = varycolor(42);
	if saveFigs && noDisp
		figure('visible','off','Position', [10,10,1000,1000]);
	else
		figure('Position', [10,10,1000,1000]);
	end
	
	hold on
		
	for ix = 1:42
		plot(cip_binMid*10, cipAR_sprlAvgs(ix,:)*100','Color',colors(ix,:), 'LineWidth', 2);
	end
	
	title('PECAN - Mean Area Ratio ($T \le 0^{\circ}C$) - All Spirals','Interpreter','latex');
	
	xlabel('D (mm)');
	ylabel('Area Ratio (%)');
	set(gca,'Xscale','log');
	xlim([0.03 2]);

	
	colormap(colors)
	cbar = colorbar('Ticks',[0 1],'TickLabels',{'FirstPSD','LastPSD'});
	
	set(gca,'XMinorTick','on','YMinorTick','on');
	set(findall(gcf,'-property','FontSize'),'FontSize',28)
	grid
	set(cbar,'FontSize',14);
	
	
	if saveFigs
		set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		print([saveDirObs '/CIP-AR-Sprls_' num2str(avgTime) 's/CIP_AR-AllSprls_Temp-lte0' fNameAppnd],Ftype,Fres)
	end

end