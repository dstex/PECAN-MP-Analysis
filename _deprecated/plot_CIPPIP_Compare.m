close all;clearvars;

flight = '20150709';

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';
savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';

avgTime = 1;

plotN	= 0;
plotND	= 0;

plotN_allSprls = 1;

pipRjct		= 0;

saveFigs	= 1;
noDisp		= 0;

%% Make any directories that are needed
if saveFigs
	saveDir = [savePath flight];
	if (exist(saveDir, 'dir') ~= 7)
		mkdir(saveDir)
	end
	if (exist([saveDir '/Comparisons'], 'dir') ~= 7)
		mkdir([saveDir '/Comparisons'])
	end
	if avgTime ~= 1
		if (plotN && exist([saveDir '/Comparisons/N_CIP-PIP_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/Comparisons/N_CIP-PIP_' num2str(avgTime) 's'])
		end
		if (plotND && exist([saveDir '/Comparisons/ND_CIP-PIP_' num2str(avgTime) 's'], 'dir') ~= 7)
			mkdir([saveDir '/Comparisons/ND_CIP-PIP_' num2str(avgTime) 's'])
		end
	else
		if (plotN && exist([saveDir '/Comparisons/N_CIP-PIP_1s'], 'dir') ~= 7)
			mkdir([saveDir '/Comparisons/N_CIP-PIP_1s'])
		end
		if (plotND && exist([saveDir '/Comparisons/ND_CIP-PIP_1s'], 'dir') ~= 7)
			mkdir([saveDir '/Comparisons/ND_CIP-PIP_1s'])
		end
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
PIP_bin_min = pipDataF.bin_min;
PIP_bin_size = pipDataF.bin_size;

if avgTime == 1
	CIP_conc_minR = cipDataF.conc_minR_orig;
	CIP_time_secs = cipDataF.time_secs_orig;
	PIP_conc_minR = pipDataF.conc_minR_orig;
	PIP_time_secs = pipDataF.time_secs_orig;

	time_secsFL = cipDataF.time_secsFL_orig;
	tempC = cipDataF.tempC_orig;
else
	CIP_conc_minR = cipDataF.conc_minR_avg;
	CIP_time_secs = cipDataF.time_secs_avg;
	PIP_conc_minR = pipDataF.conc_minR_avg;
	PIP_time_secs = pipDataF.time_secs_avg;

	time_secsFL = cipDataF.time_secsFL_avg;
	tempC = cipDataF.tempC_avg;
end

% Overlap region bin subset
cip_binMin_ol = CIP_bin_min(15:21);
pip_binMin_ol = PIP_bin_min(14:20);
cip_binWidth_ol = CIP_bin_size(15:21);
pip_binWidth_ol = PIP_bin_size(14:20);


sprlNames = fieldnames(CIP_time_secs);
loopVctr = 1:length(sprlNames);
% loopVctr = 1:2;

if plotN_allSprls
	cip_n_ol = [];
	pip_n_ol = [];
	fAll = figure('Position', [10,10,1200,700]);
	colors = varycolor(length(loopVctr));
	hold on
end

%% Loop through each spiral and time interval therein
for ix = loopVctr
	if strcmp(flight,'20150620') && ix == 4
		continue
	end
	% Spiral-specific subset
	fl_TimeSecs = time_secsFL.(sprlNames{ix});
	fl_tempC = tempC.(sprlNames{ix});
	pip_TimeSecs = PIP_time_secs.(sprlNames{ix});
	cip_TimeSecs = CIP_time_secs.(sprlNames{ix});
	cipConcSprl = CIP_conc_minR.(sprlNames{ix});
	pipConcSprl = PIP_conc_minR.(sprlNames{ix});
	
	if pipRjct
		rjctIx = find(pip_TimeSecs >= PIP_rjctStartT(ix) & pip_TimeSecs <= PIP_rjctEndT(ix));
		pipConcSprl(rjctIx,:) = NaN;
	end
	
	
	
	% Overlap region subset
	cipConcSprl_ol = cipConcSprl(:,15:21);
	pipConcSprl_ol = pipConcSprl(:,14:20);
	cip_n_ol = nansum(cipConcSprl_ol.*cip_binWidth_ol'./10,2); % #/L (#/cm-3)
	pip_n_ol = nansum(pipConcSprl_ol.*cip_binWidth_ol'./10,2);
	
	
	if plotN_allSprls
		figure(fAll);
		plot(cip_n_ol,cip_n_ol./pip_n_ol,'.','Color',colors(ix,:),'MarkerSize',20);
	end
	
	
	if plotN
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		
		semilogy(cip_n_ol,cip_n_ol./pip_n_ol,'b.','MarkerSize',20);
		
		
		if avgTime == 1
			title(sprintf('%s - Spiral %d',flight,ix));
		else
			title(sprintf('%s - Spiral %d - %d sec avg',flight,ix,avgTime));
		end
		hold on
		
		yL = ylim;
		ylim([0.3 yL(2)*1.1])
		
		ylabel('Ratio of CIP N_{800-1250\mum} to PIP N_{800-1250\mum}');
		xlabel('CIP N_{800-1250\mum} [cm^{-3}]');
		set(gca,'Xscale','log','XMinorTick','on','YMinorTick','on');
		plot(xlim,[1 1],'r--','LineWidth',4);
		grid
		
		set(findall(gcf,'-property','FontSize'),'FontSize',23);
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			if avgTime == 1
				print([saveDir '/Comparisons/N_CIP-PIP_1s/' flight '_N_CIP-PIP_S' num2str(ix) '_Whole'],'-dpng','-r0')
			else
				print([saveDir '/Comparisons/N_CIP-PIP_' num2str(avgTime) 's/' flight '_N_CIP-PIP_S' num2str(ix) '_' num2str(avgTime) 's_Whole'],'-dpng','-r0')
			end
		end
	end

	if plotND
		for ii = 1:length(cip_TimeSecs)
		
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
				if avgTime == 1
					print([saveDir '/Comparisons/ND_CIP-PIP_1s/' flight '_ND_CIP-PIP_S' num2str(ix) '_' num2str(insec2hhmmss(cip_TimeSecs(ii)))],'-dpng','-r0')
				else
					print([saveDir '/Comparisons/ND_CIP-PIP_' avgTime 's/' flight '_ND_CIP-PIP_S' num2str(ix) '_' num2str(insec2hhmmss(cip_TimeSecs(ii)))],'-dpng','-r0')
				end
				
			end
		end
	end
end

if plotN_allSprls
	% Get concentration of single particle for comparison
	cip_bin_mid = cipDataF.bin_mid;
	cip_sa2 = calc_sa_randombins(cip_bin_mid,25,100,64,0, 1);
	pip_bin_mid = pipDataF.bin_mid;
	pip_sa2 = calc_sa_randombins(pip_bin_mid,100,260,64,0, 1);
	cip_sa_m2 = 160*(1e-6);
	pip_sa_m2 = 1664*(1e-6);
	cip_sv_m3 = cip_sa_m2*130;
	pip_sv_m3 = pip_sa_m2*130;
	cip_sv_cm3 = cip_sv_m3*100^3;
	pip_sv_cm3 = pip_sv_m3*100^3;
	cip_snglPrt_conc = 1/cip_sv_cm3;
	pip_snglPrt_conc = 1/pip_sv_cm3;
	
	pip_snglPrt_conc_all = ones(size(pipConcSprl_ol)).*pip_snglPrt_conc;
	
	figure(fAll);
	
	if avgTime == 1
		title(sprintf('%s - All Spirals',flight));
	else
		title(sprintf('%s - All Spirals - %d sec avg',flight,avgTime));
	end
	
	yL = ylim;
	ylim([0.3 yL(2)*1.1])
	
	ylabel('Ratio of CIP N_{800-1250\mum} to PIP N_{800-1250\mum}');
	xlabel('CIP N_{800-1250\mum} [cm^{-3}]');
	set(gca,'Xscale','log','YScale','log','XMinorTick','on','YMinorTick','on');
	plot(xlim,[1 1],'r--','LineWidth',4);
	plot(cip_snglPrt_conc*[1 1],ylim,'b--','LineWidth',4);
	plot(cip_n_ol,cip_n_ol./pip_snglPrt_conc_all,'k--','LineWidth',4);
	grid
	
	
	% Create legends depending on which flight it is
	if strcmp(flight,'20150617')
		legend('S1','S2','S3','S4','S5','S6','S7','Location','eastoutside');
	elseif strcmp(flight,'20150620')
		legend('S1','S2','S3','S5','S6','S7','Location','eastoutside');
	elseif strcmp(flight,'20150702')
		legend('S1','S2','S3','Location','eastoutside');
	elseif strcmp(flight,'20150706')
		legend('S1','S2','S3','S4','S5','S6','S7','S8','Location','eastoutside');
	elseif strcmp(flight,'20150709')
		legend('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10',...
			'S11','S12','S13','S14','S15','S16','Location','eastoutside');
	end
	
	set(findall(gcf,'-property','FontSize'),'FontSize',23);
	
	if saveFigs
		set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		if avgTime == 1
			print([saveDir '/Comparisons/N_CIP-PIP_1s/' flight '_N_CIP-PIP_AllSprls_Whole'],'-dpng','-r0')
		else
			print([saveDir '/Comparisons/N_CIP-PIP_' num2str(avgTime) 's/' flight '_N_CIP-PIP_AllSprls_' num2str(avgTime) 's_Whole'],'-dpng','-r0')
		end
	end
end