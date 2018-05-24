clearvars; close all;

%% Specify various plotting/calculation parameters
flights = {'20150617','20150620','20150701','20150702','20150706','20150709'};
% flights = {'20150709'};

avgTime = 10;

zoom = 0;

contLevs = 70;
logCont = 1;

rmvMLmass = 1;

fileIdStr = ['_Fit-CIP_' num2str(avgTime) 'secAvg_12mm'];
% fileIdStr = ['_Fit-CIP_' num2str(avgTime) 'secAvg_55mm'];
titleIdStr = '';
if zoom
	fNameAppnd = '_zoom';
else
	fNameAppnd = '';
end

% Standard plots, but using hybrid (CIP obs + CIP extended) SDs
getLims = 0; % If true, run without plotting and print values to assist in setting axis limits

plotND				= 1;
plotNDallOne		= 1;
plotNDall			= 0;
plotNDtemp			= 0;
plotNDtempBinned	= 0;

plotMD				= 1;
plotMDallOne		= 1;
plotMDall			= 0;
plotMDtemp			= 0;
plotMDtempBinned	= 0;

plotTempIWC_LWC		= 0;
plotTWCextndRatio	= 0;

plotMDratioExcd		= 0; % Plot indvdl M(D) for periods where mass ratio between obs and extended is exceeded
plotNDratioExcd		= 0; % Plot indvdl N(D) for periods where mass ratio between obs and extended is exceeded
plotARatTemp        = 0;

% Plots using all spirals over all flights
plotARsprls  = 0;


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


% Initialize variables to be analyzed over all flights/spirals
cipAR_sprlAvgs = NaN(42,34); % 42 spirals, with 34 diameter bins each

iSprl = 1; % Counter for total number of spirals (used in whole-project variable concat)


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
	if saveFigs && ~getLims
		%%% Save info and directory operations for CIP Fit plots
		saveDir = [savePath flight '/CIP-Ext'];
		if (exist(saveDir, 'dir') ~= 7)
			mkdir(saveDir)
		end

		if (plotND && exist([saveDir '/CIP-ND-avg_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-ND-avg_' num2str(avgTime) 's'])
		end
		if (plotNDallOne && exist([saveDir '/CIP-ND_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-ND_' num2str(avgTime) 's'])
		end
		if (plotNDall && exist([saveDir '/CIP-ND-All_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-ND-All_' num2str(avgTime) 's'])
		end
		if (plotMD && exist([saveDir '/CIP-MD-avg_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-MD-avg_' num2str(avgTime) 's'])
		end
		if (plotMDallOne && exist([saveDir '/CIP-MD_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-MD_' num2str(avgTime) 's'])
		end
		if (plotMDall && exist([saveDir '/CIP-MD-All_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-MD-All_' num2str(avgTime) 's'])
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
		if (plotTempIWC_LWC && exist([saveDir '/CIP-IWC-LWC-Temp_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-IWC-LWC-Temp_' num2str(avgTime) 's'])
		end
		
		%%% Save info and directory operations for multi-flight plots
		saveDirObs = savePath;
		if (exist(saveDirObs, 'dir') ~= 7)
			mkdir(saveDirObs)
		end
		
		if (plotARsprls && exist([saveDirObs '/CIP-AR-Sprls_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDirObs '/CIP-AR-Sprls_' num2str(avgTime) 's'])
		end
	end

	%%% Save a list of variables we want to keep across all iterations
	if iFlt == 1
		initialVars = who;
		initialVars{end+1} = 'initialVars';
	end
	
	loopVctr = 1:length(sprlNames);
% 	loopVctr = [16];

	%% Specify flight-specific plotting parameters
	tempRangeAll = [-18.5 22];

	if zoom
		NDLim = [1e-7 30];
		NDavgLim = [3e-4 5];
		MDLim = [1e-10 5e-5];
		MDavgLim = [1e-9 1e-5];
	else
		NDLim = [1e-7 30];
		NDavgLim = [1e-6 5];
		MDLim = [1e-10 5e-5];
		MDavgLim = [1e-9 1e-5];
	end
	NDLogCLim = [-5 2];
	NDLogLim = [-7 2];
	MDLogCLim = [-8 -4];
	MDLogLim = [-10 -4];
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
		if getLims
			maxConc = [];
			minConc = [];
		end
		for ix = loopVctr
			
			cipConc = nanmean(cipConc_cm4_hybrid_igf.(sprlNames{ix}),1);
			
			if getLims
				if zoom
					tmpMean = cipConc(:,1:37);%1:37 gives us 0-2.6mm bins
					maxConc = [maxConc; nanmax(tmpMean(tmpMean>0))];
					minConc = [minConc; nanmin(tmpMean(tmpMean>0))];
				else
					tmpMean = cipConc(:,1:53);%1:53 gives us 0-5.8mm bins
					maxConc = [maxConc; nanmax(tmpMean(tmpMean>0))];
					minConc = [minConc; nanmin(tmpMean(tmpMean>0))];
				end
			else
				if saveFigs && noDisp
					figure('visible','off','Position', [10,10,1200,700]);
				else
					figure('Position', [10,10,1200,700]);
				end
				
				stairs(cipExt_binMin_mm(1:numObs_bins), cipConc(1:numObs_bins)', 'Color',[0 0 0.3], 'LineWidth', 2);
				hold on
				stairs(cipExt_binMin_mm(numObs_bins:end), cipConc(numObs_bins:end)','b', 'LineWidth', 2);
				
				title([flight ' - Spiral ' num2str(ix) ' - CIP (Extended) (using ' num2str(avgTime) 's Avg PSDs)']);
				
				xlabel('D (mm)');
				ylabel('N(D) (cm^{-4})');
				set(gca,'Yscale','log');
				set(gca,'Xscale','log');
				if ~isempty(diamLim)
					xlim(diamLim);
				end
				if ~isempty(NDavgLim)
					ylim(NDavgLim);
				end
				set(gca,'XMinorTick','on','YMinorTick','on');
				set(findall(gcf,'-property','FontSize'),'FontSize',28)
				grid
				
				
				if saveFigs
					tightfig(gcf);
					set(gcf,'Units','Inches');
					pos = get(gcf,'Position');
					set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
					print(sprintf('%s/CIP-ND-avg_%ds/%s_CIP_ND_%ds_S%02d%s',saveDir,avgTime,flight,avgTime,ix,fNameAppnd),Ftype,Fres)
				end
			end
		end
		if getLims
			fprintf('minConc = %.2e\n',nanmin(minConc))
			fprintf('mean_minConc = %.2e\n',nanmean(minConc))
			fprintf('maxConc = %.2g\n',nanmax(maxConc))
		end
	end
	if plotNDallOne
		if getLims
			maxConc = [];
			minConc = [];
		end
		for ix = loopVctr

			cipConc = cipConc_cm4_hybrid_igf.(sprlNames{ix});

			colors = varycolor(size(cipConc,1));
			if getLims
				if zoom
					tmpMean = cipConc(:,1:37);%1:37 gives us 0-2.6mm bins
					maxConc = [maxConc; nanmax(tmpMean(tmpMean>0))];
					minConc = [minConc; nanmin(tmpMean(tmpMean>0))];
				else
					tmpMean = cipConc(:,1:53);%1:53 gives us 0-5.8mm bins
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
					stairs(cipExt_binMin_mm(1:numObs_bins), cipConc(icx,1:numObs_bins)','Color',colors(icx,:), 'LineWidth', 2,'DisplayName',num2str(icx));
					stairs(cipExt_binMin_mm(numObs_bins:end), cipConc(icx,numObs_bins:end)','Color',colors(icx,:), 'LineWidth', 0.5,'DisplayName',[num2str(icx) 'ext']);
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
					print(sprintf('%s/CIP-ND_%ds/%s_CIP_ND_%ds_S%02d_all%s',saveDir,avgTime,flight,avgTime,ix,fNameAppnd),Ftype,Fres)
				end
			end
		end
		if getLims
			fprintf('minConc = %.2e\n',nanmin(minConc))
			fprintf('mean_minConc = %.2e\n',nanmean(minConc))
			fprintf('maxConc = %.2g\n',nanmax(maxConc))
		end
	end
	if plotNDall
		for ix = loopVctr
			cipConc = cipConc_cm4_hybrid_igf.(sprlNames{ix});
			cipConc_fit = cipConc_cm4_ext_igf.(sprlNames{ix});
			cipTime = sDistF.time_secs_avg.(sprlNames{ix});
			tempC = sDistF.tempC_avg.(sprlNames{ix});
			
			for iPsd = 1:size(cipConc,1)
				psdTime = datestr(cipTime(iPsd)/3600/24,'HH:MM:SS');
				if saveFigs && noDisp
					figure('visible','off','Position', [10,10,1200,700]);
				else
					figure('Position', [10,10,1200,700]);
				end
				
				stairs(cipExt_binMin_mm(1:numObs_bins), cipConc(iPsd,1:numObs_bins)', 'b', 'LineWidth', 2);
				hold on
				stairs(cipExt_binMin_mm(numObs_bins:end), cipConc(iPsd,numObs_bins:end)','r', 'LineWidth', 2);
				plot(cipExt_binMid_mm, cipConc_fit(iPsd,:)','k','LineWidth', 2);

				title(sprintf('%s - Spiral %d - %s UTC |  %.1f%cC  | #%d - CIP (Extended)',flight,ix,psdTime,tempC(iPsd),char(176),iPsd));
				
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
					print(sprintf('%s/CIP-ND-All_%ds/%s_CIP_ND_S%02d_%02d%s',saveDir,avgTime,flight,ix,iPsd,fNameAppnd),Ftype,Fres)
				end
			end
		end
	end
	if plotNDtemp
		if any(isinf(diamLim)) || isempty(diamLim)
			diamLim = [min(cipExt_binMid_mm) max(cipExt_binMid_mm)];
		end
		for ix = loopVctr

			cipConc = cipConc_cm4_hybrid_igf.(sprlNames{ix});
			cipConc(cipConc == 0) = NaN;

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
			contourf(cipExt_binMid_mm,time_secs/24/3600,log10(cipConc),linspace(NDLogLim(1),NDLogLim(2),contLevs),'LineColor','none');
% 			h = pcolor(cipExt_binMid_mm,time_secs/24/3600,log10(cipConc));
% 			set(h,'EdgeColor','none')
			xlabel('D (mm)');
			ylabel('Time');
			set(ax,'XMinorTick','on');
			colormap(HomeyerRainbow(contLevs));
			c=colorbar;
			set(c,'Location','southoutside');
			ylabel(c,'log_{10}N(D) (cm^{-4})');
			set(ax, 'CLim', NDLogCLim);
			
			if logCont
				set(ax,'XScale','log');
				labLocFac = 0.83;
				scaleStr = '_log';
				diamLim = [0.15 diamLim(2)];
			else
				labLocFac = 0.95;
				scaleStr = '';
			end
			set(ax,'XLim',diamLim);

			% Plot ML top/bottom locations and annotate with temp
			topStr = sprintf('%.3f %cC',mlTopTemp(ix),char(176));
			botStr = sprintf('%.3f %cC',mlBotTemp(ix),char(176));
			hold on
			plot(diamLim, [1 1]*mlTopTime(ix)/24/3600,'k--')
			tMT = text(diamLim(2)*labLocFac,mlTopTime(ix)/24/3600,topStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
			plot(diamLim, [1 1]*mlBotTime(ix)/24/3600,'k--')
			tMB = text(diamLim(2)*labLocFac,mlBotTime(ix)/24/3600,botStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
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
			set(gca,'layer','top'); % Puts tick marks and such on top of the pcolor

			if saveFigs
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				print(sprintf('%s/CIP-ND-Temp_%ds/%s_CIP_ND-Temp_%ds_S%02d%s%s',saveDir,avgTime,flight,avgTime,ix,fNameAppnd,scaleStr),'-dpng',Fres)
			end
		end
	end
	if plotNDtempBinned
		if getLims
			maxConc = [];
			minConc = [];
		end
		for ix = loopVctr

			cipConcSprl = cipConc_cm4_hybrid_igf.(sprlNames{ix});

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

					stairs(cipExt_binMin_mm(1:numObs_bins), nanmean(cipConc(:,1:numObs_bins),1)', 'Color',[0 0 0.3], 'LineWidth', 2);
					hold on
					stairs(cipExt_binMin_mm(numObs_bins:end), nanmean(cipConc(:,numObs_bins:end),1)','b', 'LineWidth', 2);

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
						print(sprintf('%s/CIP-ND-TempBinned/%s_CIP_ND-TempB_S%02d_%d_%ddegC%s',saveDir,flight,ix,iii,tempBins(iii),fNameAppnd),Ftype,Fres)
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
		if getLims
			maxMass = [];
			minMass = [];
		end
		for ix = loopVctr
			cipMass = cipMass_gcm4_hybrid_igf.(sprlNames{ix});
			
			% Set any mass values in the ML to NaN
			if rmvMLmass
				tempC = sDistF.tempC_avg.(sprlNames{ix});
				if ~isnan(mlTopTemp(ix)) && ~isnan(mlBotTemp(ix))
					mlIXs = find(tempC <= mlBotTemp(ix) & tempC >= mlTopTemp(ix));
					cipMass(mlIXs,:) = NaN;
				end
			end
			
			meanMass = nanmean(cipMass,1);
			
			if getLims
				if zoom
					tmpMean = meanMass(:,1:37);%1:37 gives us 0-2.6mm bins
					maxMass = [maxMass; nanmax(tmpMean(tmpMean>0))];
					minMass = [minMass; nanmin(tmpMean(tmpMean>0))];
				else
					tmpMean = meanMass(:,1:53);%1:53 gives us 0-5.8mm bins
					maxMass = [maxMass; nanmax(tmpMean(tmpMean>0))];
					minMass = [minMass; nanmin(tmpMean(tmpMean>0))];
				end
			else
				if saveFigs && noDisp
					figure('visible','off','Position', [10,10,1200,700]);
				else
					figure('Position', [10,10,1200,700]);
				end
				set(gcf,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
				
				yyaxis left
				stairs(cipExt_binMin_mm(1:numObs_bins), meanMass(1:numObs_bins)', 'Color',[0 0 0.3], 'LineWidth', 2);
				hold on
				stairs(cipExt_binMin_mm(numObs_bins:end), meanMass(numObs_bins:end)','b-', 'LineWidth', 2);
				
				title([flight ' - Spiral ' num2str(ix) ' - CIP (Extended) (using ' num2str(avgTime) 's Avg PSDs)']);
				
				xlabel('D (mm)');
				ylabel('M(D) (g cm^{-4})');
				set(gca,'Yscale','log');
				set(gca,'Xscale','log');
				if ~isempty(diamLim)
					xlim(diamLim);
				end
				if ~isempty(MDavgLim)
					ylim(MDavgLim);
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
				plot(cipExt_binMin_mm(1:numObs_bins),massCDF(1:numObs_bins)*100,'Color',[0.3 0 0],'LineWidth',2);
				hold on
				plot(cipExt_binMin_mm(numObs_bins:end),massCDF(numObs_bins:end)*100,'r-','LineWidth',2);
				ylabel('M(D) CDF (%)');
				
				
				if saveFigs
					tightfig(gcf);
					set(gcf,'Units','Inches');
					pos = get(gcf,'Position');
					set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
					print(sprintf('%s/CIP-MD-avg_%ds/%s_CIP_MD_%ds_S%02d%s',saveDir,avgTime,flight,avgTime,ix,fNameAppnd),Ftype,Fres)
				end
			end
		end
		if getLims
			fprintf('minMass = %.2e\n',nanmin(minMass))
			fprintf('mean_minMass = %.2e\n',nanmean(minMass))
			fprintf('maxMass = %.2g\n',nanmax(maxMass))
		end
	end
	if plotMDallOne
		if getLims
			maxMass = [];
			minMass = [];
		end
		for ix = loopVctr
			cipMass = cipMass_gcm4_hybrid_igf.(sprlNames{ix});
			
			% Set any mass values in the ML to NaN
			if rmvMLmass
				tempC = sDistF.tempC_avg.(sprlNames{ix});
				if ~isnan(mlTopTemp(ix)) && ~isnan(mlBotTemp(ix))
					mlIXs = find(tempC <= mlBotTemp(ix) & tempC >= mlTopTemp(ix));
					cipMass(mlIXs,:) = NaN;
				end
			end

			colors = varycolor(size(cipMass,1));
			if getLims
				if zoom
					tmpMean = cipMass(:,1:37);%1:37 gives us 0-2.6mm bins
					maxMass = [maxMass; nanmax(tmpMean(tmpMean>0))];
					minMass = [minMass; nanmin(tmpMean(tmpMean>0))];
				else
					tmpMean = cipMass(:,1:53);%1:53 gives us 0-5.8mm bins
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
					stairs(cipExt_binMin_mm(1:numObs_bins), cipMass(icx,1:numObs_bins)','Color',colors(icx,:), 'LineWidth', 2,'DisplayName',num2str(icx));
					stairs(cipExt_binMin_mm(numObs_bins:end), cipMass(icx,numObs_bins:end)','Color',colors(icx,:), 'LineWidth', 0.5,'DisplayName',[num2str(icx) 'ext']);
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
					print(sprintf('%s/CIP-MD_%ds/%s_CIP_MD_%ds_S%02d_all%s',saveDir,avgTime,flight,avgTime,ix,fNameAppnd),Ftype,Fres)
				end
			end
		end
		if getLims
			fprintf('minMass = %.2e\n',nanmin(minMass))
			fprintf('mean_minMass = %.2e\n',nanmean(minMass))
			fprintf('maxMass = %.2g\n',nanmax(maxMass))
		end
	end
	if plotMDall
		for ix = loopVctr
			cipMass = cipMass_gcm4_hybrid_igf.(sprlNames{ix});
			cipMass_fit = cipMass_gcm4_ext_igf.(sprlNames{ix});
			cipTime = sDistF.time_secs_avg.(sprlNames{ix});
			tempC = sDistF.tempC_avg.(sprlNames{ix});
			
			% Set any mass values in the ML to NaN
			if rmvMLmass
				if ~isnan(mlTopTemp(ix)) && ~isnan(mlBotTemp(ix))
					mlIXs = find(tempC <= mlBotTemp(ix) & tempC >= mlTopTemp(ix));
					cipMass(mlIXs,:) = NaN;
					cipMass_fit(mlIXs,:) = NaN;
				end
			end
			
			for iPsd = 1:size(cipMass,1)
				if all(isnan(cipMass(iPsd,:)))
					continue
				end
				psdTime = datestr(cipTime(iPsd)/3600/24,'HH:MM:SS');
				if saveFigs && noDisp
					figure('visible','off','Position', [10,10,1200,700]);
				else
					figure('Position', [10,10,1200,700]);
				end
				
				stairs(cipExt_binMin_mm(1:numObs_bins), cipMass(iPsd,1:numObs_bins)', 'b', 'LineWidth', 2);
				hold on
				stairs(cipExt_binMin_mm(numObs_bins:end), cipMass(iPsd,numObs_bins:end)','r', 'LineWidth', 2);
				plot(cipExt_binMid_mm, cipMass_fit(iPsd,:)','k','LineWidth', 2);

				title(sprintf('%s - Spiral %d - %s UTC |  %.1f%cC  | #%d - CIP (Extended)',flight,ix,psdTime,tempC(iPsd),char(176),iPsd));
				
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
				grid
				
				
				if saveFigs
					tightfig(gcf);
					set(gcf,'Units','Inches');
					pos = get(gcf,'Position');
					set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
					print(sprintf('%s/CIP-MD-All_%ds/%s_CIP_MD_S%02d_%02d%s',saveDir,avgTime,flight,ix,iPsd,fNameAppnd),Ftype,Fres)
				end
			end
		end
	end
	if plotMDtemp
		if any(isinf(diamLim)) || isempty(diamLim)
			diamLim = [min(cipExt_binMid_mm) max(cipExt_binMid_mm)];
		end
		for ix = loopVctr
			cipMass = cipMass_gcm4_hybrid_igf.(sprlNames{ix});
			cipMass(cipMass == 0) = NaN;
			time_secs = sDistF.time_secs_avg.(sprlNames{ix});
			
			tempC = sDistF.tempC_orig.(sprlNames{ix});
			time_fl = sDistF.time_secsFL_orig.(sprlNames{ix});

			% Set any mass values in the ML to NaN
			if rmvMLmass
				if ~isnan(mlTopTemp(ix)) && ~isnan(mlBotTemp(ix))
					tempCavg = sDistF.tempC_avg.(sprlNames{ix});
					mlIXs = find(tempCavg <= mlBotTemp(ix) & tempCavg >= mlTopTemp(ix));
					cipMass(mlIXs,:) = NaN;
				end
			end

			if saveFigs && noDisp
				figure('visible','off','Position', [10,10,1500,1000]);
			else
				figure('Position', [10,10,1500,1000]);
			end

			set(gcf,'defaultAxesColorOrder',[0 0 0; 0 0 0]);

			yyaxis left
			ax = gca;
			contourf(cipExt_binMid_mm,time_secs/24/3600,log10(cipMass),linspace(MDLogLim(1),MDLogLim(2),contLevs),'LineColor','none');
% 			h = pcolor(cipExt_binMid_mm,time_secs/24/3600,log10(mass_twc));
% 			set(h,'EdgeColor','none')
			xlabel('D (mm)');
			ylabel('Time');
			set(ax,'XMinorTick','on');
			colormap(HomeyerRainbow(contLevs));
			c=colorbar;
			set(c,'Location','southoutside');
			ylabel(c,'log_{10}M(D) (g cm^{-4})');
			set(ax, 'CLim', MDLogCLim);
			if logCont
				set(ax,'XScale','log');
				labLocFac = 0.83;
				scaleStr = '_log';
				diamLim = [0.15 diamLim(2)];
			else
				labLocFac = 0.95;
				scaleStr = '';
			end
			set(ax,'XLim',diamLim);

			% Plot ML top/bottom locations and annotate with temp
			topStr = sprintf('%.3f %cC',mlTopTemp(ix),char(176));
			botStr = sprintf('%.3f %cC',mlBotTemp(ix),char(176));
			hold on
			plot(diamLim, [1 1]*mlTopTime(ix)/24/3600,'k--')
			tMT = text(diamLim(2)*labLocFac,mlTopTime(ix)/24/3600,topStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
			plot(diamLim, [1 1]*mlBotTime(ix)/24/3600,'k--')
			tMB = text(diamLim(2)*labLocFac,mlBotTime(ix)/24/3600,botStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
			hold off

			% Swap time direction depending on spiral direction and define top/bottom spiral temps
			if tempC(1) < tempC(end)
				set(ax,'YDir','reverse');
				topT = floor(tempC(1));
				botT = ceil(tempC(end));
			else
				topT = floor(tempC(end));
				botT = ceil(tempC(1));
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
				[tmpMin,index(iT)] = min(abs(tempC-yTemp(iT)));
				% fprintf('Deviation from desired %.1f: %.2f\n',yTemp(iT),tmpMin);
				if (index(iT) == 1 ||  index(iT) == length(tempC) || tmpMin > 0.3 )
					index(iT) = NaN;
				end
			end
			indexFinal = index(~isnan(index));
			time_fl_s = time_fl/24/3600;

			if(tempC(1) < tempC(end))

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
			set(gca,'layer','top'); % Puts tick marks and such on top of the pcolor

			if saveFigs
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				print(sprintf('%s/CIP-MD-Temp_%ds/%s_CIP_MD-Temp_%ds_S%02d%s%s',saveDir,avgTime,flight,avgTime,ix,fNameAppnd,scaleStr),'-dpng',Fres)
			end
		end
	end
	if plotMDtempBinned
		if getLims
			maxMass = [];
			minMass = [];
		end
		for ix = loopVctr
			cipMassSprl = cipMass_gcm4_hybrid_igf.(sprlNames{ix});
			tempC = sDistF.tempC_avg.(sprlNames{ix});

			% Set any mass values in the ML to NaN
			if rmvMLmass
				if ~isnan(mlTopTemp(ix)) && ~isnan(mlBotTemp(ix))
					mlIXs = find(tempC <= mlBotTemp(ix) & tempC >= mlTopTemp(ix));
					cipMassSprl(mlIXs,:) = NaN;
				end
			end
			
			tempCrnd = NaN(length(tempC),1);
			for ii=1:length(tempC)
				if tempC(ii) < 0
					tempCrnd(ii) = ceil(tempC(ii));
				elseif tempC(ii) > 0
					tempCrnd(ii) = floor(tempC(ii));
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
					stairs(cipExt_binMin_mm(1:numObs_bins), meanMass(1:numObs_bins)', 'Color',[0 0 0.3], 'LineWidth', 2);
					hold on
					stairs(cipExt_binMin_mm(numObs_bins:end), meanMass(numObs_bins:end)', 'b-', 'LineWidth', 2);

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
					plot(cipExt_binMin_mm(1:numObs_bins),massCDF(1:numObs_bins)*100,'Color',[0.3 0 0],'LineWidth',2);
					hold on
					plot(cipExt_binMin_mm(numObs_bins:end),massCDF(numObs_bins:end)*100,'r-','LineWidth',2);
					ylabel('M(D) CDF (%)');



					if saveFigs
						tightfig(gcf);
						set(gcf,'Units','Inches');
						pos = get(gcf,'Position');
						set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
						print(sprintf('%s/CIP-MD-TempBinned/%s_CIP_MD-TempB_S%02d_%d_%ddegC%s',saveDir,flight,ix,iii,tempBins(iii),fNameAppnd),Ftype,Fres)
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
			%contourf(cip_binMid_mm,time_secs/24/3600,cipARat,linspace(ARlim(1),ARlim(2),contLevs),'LineColor','none');
			h = pcolor(cip_binMid_mm',time_secs/24/3600,cipARat);
			set(h,'EdgeColor','none')
			xlabel('D (mm)');
			ylabel('Time');
			set(ax,'XMinorTick','on');
			colormap(HomeyerRainbow(contLevs));
			c=colorbar;
			set(c,'Location','southoutside');
			ylabel(c,'Area Ratio (%)');
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
				print(sprintf('%s/CIP-ARat-Temp_%ds/%s_CIP_ARat-Temp_%ds_S%02d%s',saveDir,avgTime,flight,avgTime,ix,fNameAppnd),Ftype,Fres)
			end
		end
	end
	
	if plotTempIWC_LWC
		if getLims
			maxWC = [];
			minWC = [];
		end
		for ix = loopVctr
			iceFlg = iceFlag.(sprlNames{ix});
			
			nanMassIWC_obs = find(all(isnan(cipMassIWC_gcm4_hybrid_igf.(sprlNames{ix})),2));
			cipIWC_obs = nansum((cipMassIWC_gcm4_hybrid_igf.(sprlNames{ix})(:,1:numObs_bins).*cipExt_binwidth_cm(1:numObs_bins)').*1e6,2); % g m-3
			cipIWC_obs(nanMassIWC_obs) = NaN; % Set times with NaN in all contributing mass bins to NaN
			cipIWC_obs(iceFlg == 0) = 0;
			
			nanMassIWC_all = find(all(isnan(cipMassIWC_gcm4_hybrid_igf.(sprlNames{ix})),2));
			cipIWC_all = nansum((cipMassIWC_gcm4_hybrid_igf.(sprlNames{ix}).*cipExt_binwidth_cm').*1e6,2); % g m-3
			cipIWC_all(nanMassIWC_all) = NaN; % Set times with NaN in all contributing mass bins to NaN
			cipIWC_all(iceFlg == 0) = 0;
			
			nanMassLWC_obs = find(all(isnan(cipMassLWC_gcm4_hybrid_igf.(sprlNames{ix})),2));
			cipLWC_obs = nansum((cipMassLWC_gcm4_hybrid_igf.(sprlNames{ix})(:,1:numObs_bins).*cipExt_binwidth_cm(1:numObs_bins)').*1e6,2); % g m-3
			cipLWC_obs(nanMassLWC_obs) = NaN; % Set times with NaN in all contributing mass bins to NaN
			cipLWC_obs(iceFlg == 1) = 0;
			
			nanConc = find(all(isnan(cipConc_cm4_hybrid_igf.(sprlNames{ix})),2));
			cipConc = nansum(cipConc_cm4_hybrid_igf.(sprlNames{ix}),2);
			cipConc(nanConc) = NaN; % Set times with NaN in all contributing bins to NaN
			
			time_secs = sDistF.time_secs_avg.(sprlNames{ix});
			tempC = sDistF.tempC_avg.(sprlNames{ix});
			
			
			if getLims
				maxWC = [maxWC; nanmax(cipIWC_obs(cipIWC_obs>0)); nanmax(cipIWC_all(cipIWC_all>0)); nanmax(cipLWC_obs(cipLWC_obs>0))];
				minWC = [minWC; nanmin(cipIWC_obs(cipIWC_obs>0)); nanmin(cipIWC_all(cipIWC_all>0)); nanmin(cipLWC_obs(cipLWC_obs>0))];
			else
				if saveFigs && noDisp
					figure('visible','off','Position', [10,10,1000,1200]);
				else
					figure('Position', [10,10,1000,1200]);
				end
				
				plot(cipLWC_obs,tempC,'r-','LineWidth',2);
				hold on
				plot(cipIWC_obs,tempC,'b-','LineWidth',2);
				plot(cipIWC_all,tempC,'c-','LineWidth',2);
				% plot(cipConc,tempC,'k-','LineWidth',2);
				ylabel(sprintf('Temperature (%cC)',char(176)));
				
				% Plot ML top/bottom locations and annotate with temp
				xl = [4e-5 10];
				xlim(xl);
				plot(xl, [1 1]*mlTopTemp(ix),'k--')
				topStr = sprintf('%.3f %cC',mlTopTemp(ix),char(176));
				tMT = text(xl(2)*0.6,mlTopTemp(ix),topStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
				plot(xl, [1 1]*mlBotTemp(ix),'k--')
				botStr = sprintf('%.3f %cC',mlBotTemp(ix),char(176));
				tMB = text(xl(2)*0.6,mlBotTemp(ix),botStr,'HorizontalAlignment','center','BackgroundColor','k','color','w');
				
				title([flight ' - Spiral ' num2str(ix) ' - CIP (Extended) (using ' num2str(avgTime) 's Avg PSDs)']);
				
				xlabel('Water Content (g m^{-3})');
				set(gca,'Xscale','log');
				set(gca,'YDir','reverse');
				
				legend('LWC Obs','IWC Obs','IWC Obs+Ext','Location','northwest');
				set(gca,'XMinorTick','on','YMinorTick','on');
				set(findall(gcf,'-property','FontSize'),'FontSize',28);
				set(tMB,'FontSize',14);
				set(tMT,'FontSize',14);
				grid;
				
				if saveFigs
					set(gcf,'Units','Inches');
					pos = get(gcf,'Position');
					set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
					print(sprintf('%s/CIP-IWC-LWC-Temp_%ds/%s_CIP-IWC-LWC-Temp_%ds_S%02d%s',saveDir,avgTime,flight,avgTime,ix,fNameAppnd),Ftype,Fres)
				end
			end
			
		end
		if getLims
			fprintf('minWC = %.2e\n',nanmin(minWC))
			fprintf('mean_minWC = %.2e\n',nanmean(minWC))
			fprintf('maxWC = %.2g\n',nanmax(maxWC))
		end
	end
	
	
	if plotTWCextndRatio
		for ix = loopVctr

			cipMass = cipMass_gcm4_hybrid_igf.(sprlNames{ix});
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
			cipMass = cipMass_gcm4_hybrid_igf.(sprlNames{ix});
			cipMass_fitOnly = cipMass_gcm4_ext_igf.(sprlNames{ix});
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
				
				tempMassObs = (cipMass(twcRexcdSprl(irx),1:numObs_bins).*cipExt_binwidth_cm(1:numObs_bins)').*1e6; % Convert to g m-3
				tempTWCObs = nansum(tempMassObs);
				tempMassExt = (cipMass(twcRexcdSprl(irx),numObs_bins:end).*cipExt_binwidth_cm(numObs_bins:end)').*1e6; % Convert to g m-3
				tempTWCExt = nansum(tempMassExt);
				
				
				% Determine average difference between the obs and fit
				tmpFitMass = cipMass_fitOnly(twcRexcdSprl(irx),1:numObs_bins);
				tmpMass = cipMass(twcRexcdSprl(irx),1:numObs_bins);
				tmpFitMass(tmpFitMass == 0) = NaN;
				tmpMass(tmpMass == 0) = NaN;
				massDiff = nanmean(tmpFitMass - tmpMass);
				
				
				stairs(cipExt_binMin_mm(1:numObs_bins), cipMass(twcRexcdSprl(irx),1:numObs_bins)','b-','LineWidth', 2);
				hold on
				stairs(cipExt_binMin_mm(numObs_bins:end), cipMass(twcRexcdSprl(irx),numObs_bins:end)','r-','LineWidth', 2);
				plot(cipExt_binMid_mm,cipMass_fitOnly(twcRexcdSprl(irx),:)','k--','LineWidth', 2);
				
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
			cipMass = cipMass_gcm4_hybrid_igf.(sprlNames{ix});
			cipMass_fitOnly = cipMass_gcm4_ext_igf.(sprlNames{ix});
			cipConc = cipConc_cm4_hybrid_igf.(sprlNames{ix});
			cipConc_fitOnly = cipConc_cm4_ext_igf.(sprlNames{ix});
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
				
				tempMassObs = (cipMass(twcRexcdSprl(irx),1:numObs_bins).*cipExt_binwidth_cm(1:numObs_bins)').*1e6; % Convert to g m-3
				tempTWCObs = nansum(tempMassObs);
				tempMassExt = (cipMass(twcRexcdSprl(irx),numObs_bins:end).*cipExt_binwidth_cm(numObs_bins:end)').*1e6; % Convert to g m-3
				tempTWCExt = nansum(tempMassExt);
				
				
				% Determine average difference between the obs and fit
				tmpFitConc = cipConc_fitOnly(twcRexcdSprl(irx),1:numObs_bins);
				tmpConc = cipConc(twcRexcdSprl(irx),1:numObs_bins);
				tmpFitConc(tmpFitConc == 0) = NaN;
				tmpConc(tmpConc == 0) = NaN;
				concDiff = nanmean(tmpFitConc - tmpConc);
				
				
				stairs(cipExt_binMin_mm(1:numObs_bins), cipConc(twcRexcdSprl(irx),1:numObs_bins)','b-','LineWidth', 2);
				hold on
				stairs(cipExt_binMin_mm(numObs_bins:end), cipConc(twcRexcdSprl(irx),numObs_bins:end)','r-','LineWidth', 2);
				plot(cipExt_binMid_mm,cipConc_fitOnly(twcRexcdSprl(irx),:)','k--','LineWidth', 2);
				
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
		
		iSprl = iSprl+1;
	end
	% 	clearvars('-except',initialVars{:});
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
		plot(cip_binMid_mm, cipAR_sprlAvgs(ix,:)*100','Color',colors(ix,:), 'LineWidth', 2);
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