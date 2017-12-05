close all;clearvars;

flight = '20150706';

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';
savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';

plotWholeSprl	= 0;

plotND			= 1;

pipRjct			= 1;

saveFigs		= 1;
noDisp			= 1;

%% Make any directories that are needed
if saveFigs
	saveDir = [savePath flight];
	if (exist(saveDir, 'dir') ~= 7)
		mkdir(saveDir)
	end
	if (exist([saveDir '/Comparisons'], 'dir') ~= 7)
		mkdir([saveDir '/Comparisons'])
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

cipDataF = load([dataPath 'mp-data/' flight '/sDist-match/sdistCI.' flight '.CIP.10secAvg.mat']);
pipDataF = load([dataPath 'mp-data/' flight '/sDist-match/sdistCI.' flight '.PIP.10secAvg.mat']);

CIP_bin_min = cipDataF.bin_min;
CIP_bin_size = cipDataF.bin_size;
CIP_conc_minR_avg = cipDataF.conc_minR_avg;
CIP_time_secs_avg = cipDataF.time_secs_avg;


PIP_bin_min = pipDataF.bin_min;
PIP_bin_size = pipDataF.bin_size;
PIP_conc_minR_avg = pipDataF.conc_minR_avg;
PIP_time_secs_avg = pipDataF.time_secs_avg;

time_secsFL_avg = cipDataF.time_secsFL_avg;
tempC_avg = cipDataF.tempC_avg;


sprlNames = fieldnames(CIP_time_secs_avg);
% loopVctr = 1:length(sprlNames);
loopVctr = 1:2;



%% Loop through each spiral and time interval therein
for ix = loopVctr
	if strcmp(flight,'20150620') && ix == 4
		continue
	end
	fl_TimeSecs = time_secsFL_avg.(sprlNames{ix});
	fl_tempC = tempC_avg.(sprlNames{ix});
	
	pip_TimeSecs = PIP_time_secs_avg.(sprlNames{ix});
	cip_TimeSecs = CIP_time_secs_avg.(sprlNames{ix});
	
	cipConcSprl = CIP_conc_minR_avg.(sprlNames{ix});
	pipConcSprl = PIP_conc_minR_avg.(sprlNames{ix});
	
	cipConcSprl_ol = cipConcSprl(:,15:21);
	pipConcSprl_ol = pipConcSprl(:,14:20);
	cip_binMin_ol = CIP_bin_min(15:21);
	pip_binMin_ol = PIP_bin_min(14:20);
	
	if pipRjct
		rjctIx = find(pip_TimeSecs >= PIP_rjctStartT(ix) & pip_TimeSecs <= PIP_rjctEndT(ix));
		pipConcSprl_ol(rjctIx,:) = NaN;
	end
	
	if plotWholeSprl
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		
		plot(nanmean(cipConcSprl_ol,1),nanmean(cipConcSprl_ol,1)./nanmean(pipConcSprl_ol,1),'b.','MarkerSize',40);
		
		title(sprintf('%s - Spiral %d',flight,ix));
		hold on
		
		yL = ylim;
		ylim([0 yL(2)])
		
		ylabel('Ratio of CIP N(D_{800-1250\mum}) to PIP N(D_{800-1250\mum})');
		xlabel('CIP N(D_{800-1250\mum}) [cm^{-4}]');
		set(gca,'Xscale','log','XMinorTick','on','YMinorTick','on');
		plot(xlim,[1 1],'r--','LineWidth',4);
		grid
		
		set(findall(gcf,'-property','FontSize'),'FontSize',23);
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/Comparisons/ND_CIP-PIP/' flight '_ND_CIP-PIP_S' num2str(ix) '_Whole'],'-dpng','-r0')
		end
	end

	for ii = 1:length(cip_TimeSecs)
		
		if plotND
			concRatio = cipConcSprl_ol(ii,:)./pipConcSprl_ol(ii,:);
			if all(isnan(concRatio))
				continue
			end
			if saveFigs && noDisp
				figure('visible','off','Position', [10,10,1200,700]);
			else
				figure('Position', [10,10,1200,700]);
			end
			
			
			plot(cipConcSprl_ol(ii,:),concRatio,'b.','MarkerSize',40);
			title(sprintf('%s - Spiral %d - %d',flight,ix,insec2hhmmss(cip_TimeSecs(ii))));
			hold on
			
			yL = ylim;
			ylim([0 yL(2)])
			
			ylabel('Ratio of CIP N(D_{800-1250\mum}) to PIP N(D_{800-1250\mum})');
			xlabel('CIP N(D_{800-1250\mum}) [cm^{-4}]');
			set(gca,'Xscale','log','XMinorTick','on','YMinorTick','on');
			plot(xlim,[1 1],'r--','LineWidth',4);
			grid
			
			set(findall(gcf,'-property','FontSize'),'FontSize',23);
			
			if saveFigs
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				print([saveDir '/Comparisons/ND_CIP-PIP/' flight '_ND_CIP-PIP_S' num2str(ix) '_' num2str(insec2hhmmss(cip_TimeSecs(ii)))],'-dpng','-r0')
			end
		end
	end
end