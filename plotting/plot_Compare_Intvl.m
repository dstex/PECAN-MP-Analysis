close all;clearvars;

flight = '20150706';

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';
savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';

timeStep = 10;

plotMD		= 0;
plotTWC		= 0;
plotND		= 1;

pipRjct		= 1;

saveFigs	= 1;
noDisp		= 1;

%% Make any directories that are needed
if saveFigs
	saveDir = [savePath flight];
	if (exist(saveDir, 'dir') ~= 7)
		mkdir(saveDir)
	end
	if (exist([saveDir '/Comparisons'], 'dir') ~= 7)
		mkdir([saveDir '/Comparisons'])
	end
	if (plotMD && exist([saveDir '/Comparisons/MD_CIP-PIP'], 'dir') ~= 7)
		mkdir([saveDir '/Comparisons/MD_CIP-PIP'])
	end
	if (plotTWC && exist([saveDir '/Comparisons/TWC_CIP-PIP-SEA'], 'dir') ~= 7)
		mkdir([saveDir '/Comparisons/TWC_CIP-PIP-SEA'])
	end
	if (plotND && exist([saveDir '/Comparisons/ND_CIP-PIP'], 'dir') ~= 7)
		mkdir([saveDir '/Comparisons/ND_CIP-PIP'])
	end
end

%% Import data from all relevant data files
startT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'startT');
endT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'endT');
PIP_rjctStartT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'PIP_rjctStartT');
PIP_rjctEndT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'PIP_rjctEndT');

cipDataF = load([dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.CIP.10secAvg.mat']);
pipDataF = load([dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.PIP.10secAvg.mat']);
SEAfile = [dataPath 'mp-data/' flight '/' flight '_SEA_TWC.txt'];

loadSEAcsv; % Imports SEA data (SEA_Time, SEA_TWC)

CIP_bin_min = cipDataF.bin_min;
CIP_bin_size = cipDataF.bin_size;
CIP_conc_minR_orig = cipDataF.conc_minR_orig;
CIP_mass_twc_orig = cipDataF.mass_twc_orig;
CIP_time_secs_orig = cipDataF.time_secs_orig;
CIP_twc_orig = cipDataF.twc_orig;

PIP_bin_min = pipDataF.bin_min;
PIP_bin_size = pipDataF.bin_size;
PIP_conc_minR_orig = pipDataF.conc_minR_orig;
PIP_mass_twc_orig = pipDataF.mass_twc_orig;
PIP_time_secs_orig = pipDataF.time_secs_orig;
PIP_twc_orig = pipDataF.twc_orig;

time_secsFL_orig = cipDataF.time_secsFL_orig;
tempC_orig = cipDataF.tempC_orig;


sprlNames = fieldnames(CIP_time_secs_orig);
loopVctr = 1:length(sprlNames);



%% Loop through each spiral and time interval therein
for ix = loopVctr
	sprlIx_SEA = find(SEA_Time >= startT(ix) & SEA_Time < endT(ix));
	twc_SEA = SEA_TWC(sprlIx_SEA);
	sea_TimeSecs = SEA_Time(sprlIx_SEA);
	
	fl_TimeSecs = time_secsFL_orig.(sprlNames{ix});
	fl_tempC = tempC_orig.(sprlNames{ix});
	
	pip_TimeSecs = PIP_time_secs_orig.(sprlNames{ix});
	cip_TimeSecs = CIP_time_secs_orig.(sprlNames{ix});
	
	cipMassSprl = CIP_mass_twc_orig.(sprlNames{ix});
	pipMassSprl = PIP_mass_twc_orig.(sprlNames{ix});
	
	cipTWCSprl = CIP_twc_orig.(sprlNames{ix});
	pipTWCSprl = PIP_twc_orig.(sprlNames{ix});
	
	cipConcSprl = CIP_conc_minR_orig.(sprlNames{ix});
	pipConcSprl = PIP_conc_minR_orig.(sprlNames{ix});
	
	if pipRjct
		rjctIx = find(pip_TimeSecs >= PIP_rjctStartT(ix) & pip_TimeSecs <= PIP_rjctEndT(ix));
		pipMassSprl(rjctIx,:) = NaN;
		pipTWCSprl(rjctIx,:) = NaN;
		pipConcSprl(rjctIx,:) = NaN;
	end
	
	% Create some initial conditions for our while loop
	pipTimeBeg = floor(pip_TimeSecs(1));
	pipTimeEnd = pipTimeBeg + timeStep;
	cipTimeBeg = floor(cip_TimeSecs(1));
	cipTimeEnd = cipTimeBeg + timeStep;
	seaTimeBeg = floor(sea_TimeSecs(1));
	seaTimeEnd = seaTimeBeg + timeStep;
	flTimeBeg = floor(fl_TimeSecs(1));
	flTimeEnd = flTimeBeg + timeStep;
	
	while ( (pipTimeEnd <= PIP_time_secs_orig.(sprlNames{ix})(end)) && (cipTimeEnd <= CIP_time_secs_orig.(sprlNames{ix})(end)) && (seaTimeEnd <= sea_TimeSecs(end)))
		pipPlotIx = find(pip_TimeSecs >= pipTimeBeg & pip_TimeSecs < pipTimeEnd);
		cipPlotIx = find(cip_TimeSecs >= cipTimeBeg & cip_TimeSecs < cipTimeEnd);
		seaPlotIx = find(sea_TimeSecs >= seaTimeBeg & sea_TimeSecs < seaTimeEnd);
		flPlotIx = find(fl_TimeSecs >= flTimeBeg & fl_TimeSecs < flTimeEnd);
		
		disp(['Now plotting ' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(1)))) '-' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(end))))]);
		
		flTempC = fl_tempC(flPlotIx);
		
		cipMass = cipMassSprl(cipPlotIx,:);
		pipMass = pipMassSprl(pipPlotIx,:);
		
		cipTWC = cipTWCSprl(cipPlotIx);
		pipTWC = pipTWCSprl(pipPlotIx);
		seaTWC = twc_SEA(seaPlotIx);
		
		seaTWC_mean = nanmean(seaTWC);
		seaMeanPlot = ones(length(PIP_bin_min),1).*seaTWC_mean;
		
		cipConc = cipConcSprl(cipPlotIx,:);
		pipConc = pipConcSprl(pipPlotIx,:);
		
		cipMeanMass = nanmean(cipMass,1);
		pipMeanMass = nanmean(pipMass,1);
		
		
		CIPmassCDF = zeros(size(cipMeanMass));
		for ii=1:length(cipMeanMass)
			CIPmassCDF(ii) = nansum(cipMeanMass(1:ii));
		end
		
		PIPmassCDF = zeros(size(pipMeanMass));
		for ii=1:length(pipMeanMass)
			PIPmassCDF(ii) = nansum(pipMeanMass(1:ii));
		end
		
		CIPmassCDF = (CIPmassCDF.*(CIP_bin_size'./10))./1e-6; %([g/cm4]*[mm/10])/[1e-6 m3/cm3] -> [g/m3]
		PIPmassCDF = (PIPmassCDF.*(PIP_bin_size'./10))./1e-6; %([g/cm4]*[mm/10])/[1e-6 m3/cm3] -> [g/m3]
		
		if plotMD
			if saveFigs && noDisp
				figure('visible','off','Position', [10,10,1200,700]);
			else
				figure('Position', [10,10,1200,700]);
			end
			
			yyaxis left
			stairs(CIP_bin_min,cipMeanMass,'b-','LineWidth',2);
			if pipRjct
				title([flight ' - Spiral ' num2str(ix) ' - ' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(1)))) '-' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(end))))]);
			else
				title([flight ' - Spiral ' num2str(ix) ' - w/All PIP Data - ' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(1)))) '-' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(end))))]);
			end
			hold on
			stairs(PIP_bin_min,pipMeanMass,'r-','LineWidth',2);
			xlabel('D [mm]');
			ylabel('M_{twc}(D) [g cm^{-4}]');
			set(gca,'Yscale','log','Xscale','log','XMinorTick','on','YMinorTick','on');
			grid
			
			yyaxis right
			plot(CIP_bin_min,CIPmassCDF,'c-','LineWidth',2);
			ylabel('Cumulative Mass_{twc} [g m^{-3}]');
			hold on
			plot(PIP_bin_min,PIPmassCDF,'m-','LineWidth',2);
			plot(PIP_bin_min,seaMeanPlot,'k--','LineWidth',2);
			
			legend('CIP M(D)','PIP M(D)','CIP CDF','PIP CDF','SEA Avg TWC','Location','southwest');
			
			set(findall(gcf,'-property','FontSize'),'FontSize',23);
			
			if saveFigs
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				if pipRjct
					print([saveDir '/Comparisons/MD_CIP-PIP/' flight '_MD_CIP-PIP_' num2str(timeStep) 's_S' num2str(ix) '_' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(end))))],'-dpdf','-r0')
				else
					print([saveDir '/Comparisons/MD_CIP-PIP/' flight '_MD_CIP-PIPall_' num2str(timeStep) 's_S' num2str(ix) '_' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(end))))],'-dpdf','-r0')
				end
			end
		end
		
		if plotTWC
			if saveFigs && noDisp
				figure('visible','off','Position', [10,10,1200,700]);
			else
				figure('Position', [10,10,1200,700]);
			end
			
			yyaxis left
			plot(cip_TimeSecs(cipPlotIx)/3600/24,cipTWC,'b-','LineWidth',1);
			hold on
			plot(pip_TimeSecs(pipPlotIx)/3600/24,pipTWC,'r-','LineWidth',1);
			plot(sea_TimeSecs(seaPlotIx)/3600/24,seaTWC,'m-','LineWidth',1);
			ylabel('TWC [g m^{-3}]');
			
			yyaxis right
			plot(fl_TimeSecs(flPlotIx)/3600/24,flTempC,'g-','LineWidth',1);
			datetickzoom('x','HH:MM:SS');
			set(gca,'XMinorTick','on','YMinorTick','on','XTickLabelRotation',45,'YDir','reverse');
			grid
			
			xlabel('Time [UTC]');
			ylabel(sprintf('Temperature (%cC)',char(176)));
			legend('CIP','PIP','SEA','Temp','Location','best');
			if pipRjct
				title([flight ' - Spiral ' num2str(ix) ' - ' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(1)))) '-' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(end))))]);
			else
				title([flight ' - Spiral ' num2str(ix) ' - w/All PIP Data - ' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(1)))) '-' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(end))))]);
			end
			set(findall(gcf,'-property','FontSize'),'FontSize',23);
			
			if saveFigs
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				if pipRjct
					print([saveDir '/Comparisons/TWC_CIP-PIP-SEA/' flight '_TWC_CIP-PIP-SEA_' num2str(timeStep) 's_S' num2str(ix) '_' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(end))))],'-dpdf','-r0')
				else
					print([saveDir '/Comparisons/TWC_CIP-PIP-SEA/' flight '_TWC_CIP-PIPall-SEA_' num2str(timeStep) 's_S' num2str(ix) '_' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(end))))],'-dpdf','-r0')
				end
			end
		end
		
		if plotND
			if saveFigs && noDisp
				figure('visible','off','Position', [10,10,1200,700]);
			else
				figure('Position', [10,10,1200,700]);
			end
			
			stairs(CIP_bin_min,nanmean(cipConc,1),'b-','LineWidth',2);
			if pipRjct
				title([flight ' - Spiral ' num2str(ix) ' - ' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(1)))) '-' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(end))))]);
			else
				title([flight ' - Spiral ' num2str(ix) ' - w/All PIP Data - ' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(1)))) '-' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(end))))]);
			end
			hold on
			stairs(PIP_bin_min,nanmean(pipConc,1),'r-','LineWidth',2);
			xlabel('D [mm]');
			ylabel('N(D) [cm^{-4}]');
			set(gca,'Yscale','log','Xscale','log','XMinorTick','on','YMinorTick','on');
			grid
			
			legend('CIP N(D)','PIP N(D)');
			
			set(findall(gcf,'-property','FontSize'),'FontSize',23);
			
			if saveFigs
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				if pipRjct
					print([saveDir '/Comparisons/ND_CIP-PIP/' flight '_ND_CIP-PIP_' num2str(timeStep) 's_S' num2str(ix) '_' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(end))))],'-dpdf','-r0')
				else
					print([saveDir '/Comparisons/ND_CIP-PIP/' flight '_ND_CIP-PIPall_' num2str(timeStep) 's_S' num2str(ix) '_' num2str(insec2hhmmss(pip_TimeSecs(pipPlotIx(end))))],'-dpdf','-r0')
				end
			end
		end
		
		pipTimeBeg = pipTimeEnd;
		pipTimeEnd = pipTimeBeg + timeStep;
		cipTimeBeg = cipTimeEnd;
		cipTimeEnd = cipTimeBeg + timeStep;
		seaTimeBeg = seaTimeEnd;
		seaTimeEnd = seaTimeBeg + timeStep;
		flTimeBeg = flTimeEnd;
		flTimeEnd = flTimeBeg + timeStep;
	end
end