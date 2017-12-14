close all;clearvars;

flight = '20150706';

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';
savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';

plotMD		= 1;
plotTWC		= 0;
plotND		= 1;

pipRjct		= 0;

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

% cipDataF = load([dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.CIP.10secAvg.mat']);
% pipDataF = load([dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.PIP.10secAvg.mat']);
cipDataF = load([dataPath 'mp-data/' flight '/sDist_matchBins/sdistCI.' flight '.CIP.10secAvg.mat']);
pipDataF = load([dataPath 'mp-data/' flight '/sDist_matchBins/sdistCI.' flight '.PIP.10secAvg.mat']);
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

%% Loop through each spiral and plot
for ix = loopVctr
	sprlIx_SEA = find(SEA_Time >= startT(ix) & SEA_Time < endT(ix));
	twc_SEA = SEA_TWC(sprlIx_SEA);
	sea_TimeSecs = SEA_Time(sprlIx_SEA);
	
	seaTWC_mean = nanmean(twc_SEA);
	seaMeanPlot = ones(length(PIP_bin_min),1).*seaTWC_mean;
	
	fl_TimeSecs = time_secsFL_orig.(sprlNames{ix});
	fl_tempC = tempC_orig.(sprlNames{ix});
	
	pip_TimeSecs = PIP_time_secs_orig.(sprlNames{ix});
	cip_TimeSecs = CIP_time_secs_orig.(sprlNames{ix});
	
	cipMass = CIP_mass_twc_orig.(sprlNames{ix});
	pipMass = PIP_mass_twc_orig.(sprlNames{ix});
	
	cipTWC = CIP_twc_orig.(sprlNames{ix});
	pipTWC = PIP_twc_orig.(sprlNames{ix});
	
	cipConc = CIP_conc_minR_orig.(sprlNames{ix});
	pipConc = PIP_conc_minR_orig.(sprlNames{ix});
	
	if pipRjct
		rjctIx = find(pip_TimeSecs >= PIP_rjctStartT(ix) & pip_TimeSecs <= PIP_rjctEndT(ix));
		pipMass(rjctIx,:) = NaN;
		pipTWC(rjctIx,:) = NaN;
		pipConc(rjctIx,:) = NaN;
	end
	
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
			title([flight ' - Spiral ' num2str(ix)]);
		else
			title([flight ' - Spiral ' num2str(ix) ' - w/All PIP Data']);
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
		
		legend('CIP M(D)','PIP M(D)','CIP CDF','PIP CDF','SEA Avg TWC');
		
		set(findall(gcf,'-property','FontSize'),'FontSize',23);
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			if pipRjct
				print([saveDir '/Comparisons/MD_CIP-PIP/' flight '_MD_CIP-PIP_S' num2str(ix)],'-dpdf','-r0')
			else
				print([saveDir '/Comparisons/MD_CIP-PIP/' flight '_MD_CIP-PIPall_S' num2str(ix)],'-dpdf','-r0')
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
		plot(cip_TimeSecs/3600/24,cipTWC,'b-','LineWidth',1);
		hold on
		plot(pip_TimeSecs/3600/24,pipTWC,'r-','LineWidth',1);
		plot(sea_TimeSecs/3600/24,twc_SEA,'m-','LineWidth',1);
		ylabel('TWC [g m^{-3}]');
		
		yyaxis right
		plot(fl_TimeSecs/3600/24,fl_tempC,'g-','LineWidth',1);
		datetickzoom('x','HH:MM:SS');
		set(gca,'XMinorTick','on','YMinorTick','on','XTickLabelRotation',45,'YDir','reverse');
		grid
		
		xlabel('Time [UTC]');
		ylabel(sprintf('Temperature (%cC)',char(176)));
		legend('CIP','PIP','SEA','Temp');
		if pipRjct
			title([flight ' - Spiral ' num2str(ix)]);
		else
			title([flight ' - Spiral ' num2str(ix) ' - w/All PIP Data']);
		end
		set(findall(gcf,'-property','FontSize'),'FontSize',23);
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			if pipRjct
				print([saveDir '/Comparisons/TWC_CIP-PIP-SEA/' flight '_TWC_CIP-PIP-SEA_S' num2str(ix)],'-dpdf','-r0')
			else
				print([saveDir '/Comparisons/TWC_CIP-PIP-SEA/' flight '_TWC_CIP-PIPall-SEA_S' num2str(ix)],'-dpdf','-r0')
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
			title([flight ' - Spiral ' num2str(ix)]);
		else
			title([flight ' - Spiral ' num2str(ix) ' - w/All PIP Data']);
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
				print([saveDir '/Comparisons/ND_CIP-PIP/' flight '_ND_CIP-PIP_S' num2str(ix)],'-dpdf','-r0')
			else
				print([saveDir '/Comparisons/ND_CIP-PIP/' flight '_ND_CIP-PIPall_S' num2str(ix)],'-dpdf','-r0')
			end
		end
	end
end