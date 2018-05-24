clearvars; close all;

%% Specify various plotting/calculation parameters
% flights = {'20150617','20150620','20150701','20150702','20150706','20150709'};
flights = {'20150706'};

avgTime = 10;

rmvMLmass = 1;
if rmvMLmass
	mExcldTstr = ' (ML Data Excld)';
else
	mExcldTstr = '';
end

contLevs = 64;

titleIdStr = '';

fNameAppnd = '';


% Standard plots, but using hybrid (CIP obs + CIP extended) SDs
getLims = 0; % If true, run without plotting and print values to assist in setting axis limits

plotND				= 0;
plotNDallOne		= 0;
plotNDtemp			= 0;
plotNDtempBinned	= 0;

plotMD				= 1;
plotMDallOne		= 1;
plotMDtemp			= 1;
plotMDtempBinned	= 0;


diamLim = [0.1 2.1];

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


for iFlt = 1:length(flights)
	flight = flights{iFlt};
	fprintf('\nPlotting %s...\n',flight);
	%% Load in original SD data
	if avgTime == 1
		sDistFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.CIP.10secAvg.mat']; % 1-sec averages are saved in all averaged output files
	else
		sDistFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.CIP.' num2str(avgTime) 'secAvg.mat'];
	end

	load(sDistFile);
	
	sprlNames = fieldnames(time_secs_orig); % Variable used unimportant - just needs to be one of the structs
	
	
	mlTopTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlTopTime');
	mlTopTemp = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlTopTemp');
	mlBotTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlBotTime');
	mlBotTemp = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlBotTemp');

	%% Create directories to save plots in if they don't already exist
	if saveFigs
		%%% Save info and directory operations for CIP Obs plots
		saveDir = [savePath flight '/CIP-Obs'];
		if (exist(saveDir, 'dir') ~= 7)
			mkdir(saveDir)
		end

		if (plotND && exist([saveDir '/CIP-ND-avg_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-ND-avg_' num2str(avgTime) 's'])
		end
		if (plotNDallOne && exist([saveDir '/CIP-ND_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-ND_' num2str(avgTime) 's'])
		end
		if (plotMD && exist([saveDir '/CIP-MD-avg_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-MD-avg_' num2str(avgTime) 's'])
		end
		if (plotMDallOne && exist([saveDir '/CIP-MD_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-MD_' num2str(avgTime) 's'])
		end
		if (plotNDtemp && exist([saveDir '/CIP-ND-Temp_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-ND-Temp_' num2str(avgTime) 's'])
		end
		if (plotMDtemp && exist([saveDir '/CIP-MD-Temp_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-MD-Temp_' num2str(avgTime) 's'])
		end
		if (plotNDtempBinned && exist([saveDir '/CIP-ND-TempBinned'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-ND-TempBinned'])
		end
		if (plotMDtempBinned && exist([saveDir '/CIP-MD-TempBinned'], 'dir') ~= 7)
			mkdir([saveDir '/CIP-MD-TempBinned'])
		end
	end

	%%% Save a list of variables we want to keep across all iterations
	if iFlt == 1
		initialVars = who;
		initialVars{end+1} = 'initialVars';
	end
	
	loopVctr = 1:length(sprlNames);
% 	loopVctr = [2];

	%% Specify flight-specific plotting parameters
	tempRangeAll = [-18.5 22];

	NDavgLim = [3e-4 5];
	NDLim = [1e-4 20]; %1e-5 for tempB
	MDavgLim = [5e-8 5e-6];
	MDLim = [1e-9 5e-5]; %2e-10 for tempB
	NDLogLim = [-4 2];
	MDLogLim = [-8 -4];
	
	% Switch statement to allow for flight-specific plot limits
	%{
	switch flight
		case '20150617'
			NDLim = [1e-4 10];
			MDLim = [1e-8 5e-5]; 
			NDLogLim = [-4 2];
			MDLogLim = [-8 -4];
		case '20150620'
			NDLim = [1e-4 10];
			MDLim = [1e-8 5e-5];  
			NDLogLim = [-4 2];
			MDLogLim = [-8 -4];
		case '20150701'
			NDLim = [1e-4 10];
			MDLim = [1e-8 5e-5]; 
			NDLogLim = [-4 2];
			MDLogLim = [-8 -4];
		case '20150702'
			NDLim = [1e-4 10];
			MDLim = [1e-8 5e-5]; 
			NDLogLim = [-4 2];
			MDLogLim = [-8 -4];
		case '20150706'
			NDLim = [1e-4 10];
			MDLim = [1e-8 5e-5]; 
			NDLogLim = [-4 2];
			MDLogLim = [-8 -4];
		case '20150709'
			NDLim = [1e-4 10];
			MDLim = [1e-8 5e-5]; 
			NDLogLim = [-4 2];
			MDLogLim = [-8 -4];
	end
	%}

	%% Standard plots
	if plotND
		for ix = loopVctr

			cipConc = nanmean(conc_minR_cm4_avg.(sprlNames{ix}),1);

			if saveFigs && noDisp
				figure('visible','off','Position', [10,10,1200,700]);
			else
				figure('Position', [10,10,1200,700]);
			end

			stairs(bin_min_mm, cipConc', 'Color','b', 'LineWidth', 2);

			title([flight ' - Spiral ' num2str(ix) ' - CIP (using ' num2str(avgTime) 's Avg PSDs)']);

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
				print([saveDir '/CIP-ND-avg_' num2str(avgTime) 's/' flight '_CIP_ND-avg_' num2str(avgTime) 's_S' num2str(ix) fNameAppnd],Ftype,Fres)
			end
		end
	end
	if plotNDallOne
		if getLims
			maxConc = [];
			minConc = [];
		end
		for ix = loopVctr

			cipConc = conc_minR_cm4_avg.(sprlNames{ix});

			colors = varycolor(size(cipConc,1));
			if getLims
% 				tmpMean = nanmean(cipConc,1);
				tmpMean = cipConc;
				maxConc = [maxConc; nanmax(tmpMean(tmpMean>0))];
				minConc = [minConc; nanmin(tmpMean(tmpMean>0))];
			else
				if saveFigs && noDisp
					figure('visible','off','Position', [10,10,1200,700]);
				else
					figure('Position', [10,10,1200,700]);
				end

				hold on

				for icx=1:size(cipConc,1)
					stairs(bin_min_mm, cipConc(icx,:)','Color',colors(icx,:), 'LineWidth', 2);
				end
				title([flight ' - Spiral ' num2str(ix) ' - CIP (using ' num2str(avgTime) 's Avg PSDs)']);

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
% 					tightfig(gcf);
					set(gcf,'Units','Inches');
					pos = get(gcf,'Position');
					set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
					print([saveDir '/CIP-ND_' num2str(avgTime) 's/' flight '_CIP_ND_' num2str(avgTime) 's_S' num2str(ix) fNameAppnd],Ftype,Fres)
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
			diamLim = [min(bin_mid) max(bin_mid)];
		end
		for ix = loopVctr

			tempCsprl = tempC_orig.(sprlNames{ix});
			time_fl = time_secsFL_orig.(sprlNames{ix});

			if avgTime == 1
				time_secs = time_secs_orig.(sprlNames{ix});
				cipConc = conc_minR_cm4_orig.(sprlNames{ix});
			else
				time_secs = time_secs_avg.(sprlNames{ix});
				cipConc = conc_minR_cm4_avg.(sprlNames{ix});
			end


			if saveFigs && noDisp
				figure('visible','off','Position', [10,10,1500,1000]);
			else
				figure('Position', [10,10,1500,1000]);
			end

			set(gcf,'defaultAxesColorOrder',[0 0 0; 0 0 0]);

			yyaxis left
			ax = gca;
			contourf(bin_mid_mm,time_secs/24/3600,log10(cipConc),linspace(NDLogLim(1),NDLogLim(2),contLevs),'LineColor','none');
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

			title([flight ' - Spiral ' num2str(ix) ' - CIP (using ' num2str(avgTime) 's Avg PSDs)']);

			set(findall(gcf,'-property','FontSize'),'FontSize',26)
			set(tMB,'FontSize',14);
			set(tMT,'FontSize',14);

			if saveFigs
% 				tightfig(gcf);
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

			cipConcSprl = conc_minR_cm4_avg.(sprlNames{ix});

			tempCsprl = tempC_avg.(sprlNames{ix});

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
					tmpMean = nanmean(cipConc,1);
					maxConc = [maxConc; nanmax(tmpMean(tmpMean>0))];
					minConc = [minConc; nanmin(tmpMean(tmpMean>0))];
				else
					if saveFigs && noDisp
						figure('visible','off','Position', [10,10,1500,1000]);
					else
						figure('Position', [10,10,1500,1000]);
					end

					stairs(bin_min_mm, nanmean(cipConc,1)', 'Color','b', 'LineWidth', 2);

					title(sprintf('%s - Spiral %d - CIP (%d%cC avg)',flight,ix,tempBins(iii),char(176)));

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
			cipMass = mass_twc_gcm4_avg.(sprlNames{ix});
			% Set any mass values in the ML to NaN
			if rmvMLmass
				tempC = tempC_avg.(sprlNames{ix});
				if ~isnan(mlTopTemp(ix)) && ~isnan(mlBotTemp(ix))
					mlIXs = find(tempC <= mlBotTemp(ix) & tempC >= mlTopTemp(ix));
					cipMass(mlIXs,:) = NaN;
				end
			end
			
			meanMass = nanmean(cipMass,1);

			if saveFigs && noDisp
				figure('visible','off','Position', [10,10,1200,700]);
			else
				figure('Position', [10,10,1200,700]);
			end
			set(gcf,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);

			yyaxis left
			stairs(bin_min_mm, meanMass', 'Color','b', 'LineWidth', 2);

			title([flight ' - Spiral ' num2str(ix) ' - CIP' mExcldTstr ' (using ' num2str(avgTime) 's Avg PSDs)']);

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
			plot(bin_min_mm,massCDF*100,'Color','r','LineWidth',2);
			ylabel('M(D) CDF (%)');


			if saveFigs
				tightfig(gcf);
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				print([saveDir '/CIP-MD-avg_' num2str(avgTime) 's/' flight '_CIP_MD-avg_' num2str(avgTime) 's_S' num2str(ix) fNameAppnd],Ftype,Fres)
			end
		end
	end
	if plotMDallOne
		if getLims
			maxMass = [];
			minMass = [];
		end
		for ix = loopVctr
			cipMass = mass_twc_gcm4_avg.(sprlNames{ix});

			% Set any mass values in the ML to NaN
			if rmvMLmass
				tempC = tempC_avg.(sprlNames{ix});
				if ~isnan(mlTopTemp(ix)) && ~isnan(mlBotTemp(ix))
					mlIXs = find(tempC <= mlBotTemp(ix) & tempC >= mlTopTemp(ix));
					cipMass(mlIXs,:) = NaN;
				end
			end
			
			colors = varycolor(size(cipMass,1));
			if getLims
% 				tmpMean = nanmean(cipMass,1);
				tmpMean = cipMass;
				maxMass = [maxMass; nanmax(tmpMean(tmpMean>0))];
				minMass = [minMass; nanmin(tmpMean(tmpMean>0))];
			else
				if saveFigs && noDisp
					figure('visible','off','Position', [10,10,1200,700]);
				else
					figure('Position', [10,10,1200,700]);
				end

				hold on

				for icx=1:size(cipMass,1)
					stairs(bin_min_mm, cipMass(icx,:)','Color',colors(icx,:), 'LineWidth', 2);
				end
				title([flight ' - Spiral ' num2str(ix) ' - CIP' mExcldTstr ' (using ' num2str(avgTime) 's Avg PSDs)']);

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
% 					tightfig(gcf);
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
			diamLim = [min(bin_mid) max(bin_mid)];
		end
		for ix = loopVctr
			tempCsprl = tempC_orig.(sprlNames{ix});
			time_fl = time_secsFL_orig.(sprlNames{ix});

			if avgTime == 1
				time_secs = time_secs_orig.(sprlNames{ix});	
				cipMass = mass_twc_gcm4_orig.(sprlNames{ix});
			else
				time_secs = time_secs_avg.(sprlNames{ix});
				cipMass = mass_twc_gcm4_avg.(sprlNames{ix});
			end
			
			% Set any mass values in the ML to NaN
			if rmvMLmass
				tempC = tempC_avg.(sprlNames{ix});
				if ~isnan(mlTopTemp(ix)) && ~isnan(mlBotTemp(ix))
					mlIXs = find(tempC <= mlBotTemp(ix) & tempC >= mlTopTemp(ix));
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
			contourf(bin_mid_mm,time_secs/24/3600,log10(cipMass),linspace(MDLogLim(1),MDLogLim(2),contLevs),'LineColor','none');
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


			title([flight ' - Spiral ' num2str(ix) ' - CIP' mExcldTstr ' (using ' num2str(avgTime) 's Avg PSDs)']);

			set(findall(gcf,'-property','FontSize'),'FontSize',26)
			set(tMB,'FontSize',14);
			set(tMT,'FontSize',14);

			if saveFigs
% 				tightfig(gcf);
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

			cipMassSprl = mass_twc_gcm4_avg.(sprlNames{ix});

			tempCsprl = tempC_avg.(sprlNames{ix});
			
			% Set any mass values in the ML to NaN
			if rmvMLmass
				if ~isnan(mlTopTemp(ix)) && ~isnan(mlBotTemp(ix))
					mlIXs = find(tempCsprl <= mlBotTemp(ix) & tempCsprl >= mlTopTemp(ix));
					cipMassSprl(mlIXs,:) = NaN;
				end
			end

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
					tmpMean = nanmean(cipMass,1);
					maxMass = [maxMass; nanmax(tmpMean(tmpMean>0))];
					minMass = [minMass; nanmin(tmpMean(tmpMean>0))];
				else
					if saveFigs && noDisp
						figure('visible','off','Position', [10,10,1500,1000]);
					else
						figure('Position', [10,10,1500,1000]);
					end

					meanMass = nanmean(cipMass,1);

					yyaxis left
					stairs(bin_min_mm, meanMass', 'Color','b', 'LineWidth', 2);

					title(sprintf('%s - Spiral %d - CIP%s (%d%cC avg)',flight,ix,mExcldTstr,tempBins(iii),char(176)));

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
					plot(bin_min_mm,massCDF*100,'Color','r','LineWidth',2);
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
	
% 	clearvars('-except',initialVars{:});
end