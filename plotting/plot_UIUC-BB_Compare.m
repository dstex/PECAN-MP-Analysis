% Reads in data from Bob Black's PECAN output files
% Also reads in variables from our own processed SD datasets and does any necessary conversions
% (such as unnormalizing by sample volume or binwidth) to line up with Black's data

close all;clearvars;

flight = '20150706';

plotBB_CipPip		= 0;
plotUIUC_CipPip		= 0;
plotBB_UIUC_Compare = 0;
plotCommonT			= 0;
plotSingleT			= 0;

saveFigs	= 1;
noDisp		= 1;
Ftype	= '-dpdf';

savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';

%% Read in required data
readBobBlack_SDs;

%% Times common to both probes and datasets given 6-sec averaging
cT.sprl1 = [11995;12001;12007;12331;12337;12343;12349;12355];
cT.sprl3 = 16150;
cT.sprl5 = [20469;20979;20985;20991;21081;21087;21237;21321;21327];
cT.sprl7 = [23408;23534;23540;23546];
cT.sprl8 = [23983;24187;24193;24199;24205];

cTnames = fieldnames(cT);

%% Make any directories that are needed
if saveFigs
	saveDir = [savePath flight];
	if (exist(saveDir, 'dir') ~= 7)
		mkdir(saveDir)
	end
	if (exist([saveDir '/Black-Comparisons'], 'dir') ~= 7)
		mkdir([saveDir '/Black-Comparisons'])
	end
end


%% Plot comparisons between datasets/probes
if plotBB_CipPip
	% 	for ix=1:length(sprlNames)
	for ix=7
		bb_cipConc = bb_cip_concAvg.(sprlNames{ix});
		bb_pipConc = bb_pip_concAvg.(sprlNames{ix});
		bb_time = insec2hhmmss(bb_cip_timeSecsAvg.(sprlNames{ix}));
		
		for ii = 1:size(bb_cipConc,1)
% 		for ii = 1:15
			if saveFigs && noDisp
				figure('visible','off','Position', [10,10,1200,700]);
			else
				figure('Position', [10,10,1200,700]);
			end
			
			stairs(bb_cip_binMin,bb_cipConc(ii,:),'r','LineWidth', 1.5);
			hold on
			if (ii <= size(bb_pipConc,1))
				stairs(bb_pip_binMin,bb_pipConc(ii,:),'b','LineWidth', 1.5);
				lgnd = legend('CIP','PIP');
			else
				lgnd = legend('CIP');
			end
			lgnd.FontSize = 18;
			title([flight ' - BB CIP & PIP - Spiral ' num2str(ix) ' - ' num2str(bb_time(ii)) ' - ' num2str(ii)],'FontSize',24);
			xlabel('D [mm]');
			ylabel('N(D) [cm^{-4}]');
			set(gca,'yscale','log');
% 			set(gca,'YTick',[10e-10 10e-9 10e-8 10e-7 10e-6... 
% 				10e-5 10e-4 10e-3 10e-2 10e-1 10e0 10e1 10e2 10e3 10e4 10e5 10e6 10e7]);
			grid on
			set(gca,'YMinorGrid','on','YMinorTick','on');
			ax = ancestor(gca, 'axes');
			xRule = ax.XAxis;
			yRule = ax.YAxis;
			xRule.FontSize = 24;
			yRule.FontSize = 24;
			
			if saveFigs
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				print([saveDir '/Black-Comparisons/' flight '_bb-CIP-PIP_S' num2str(ix) '-' sprintf('%2.2d',ii) '-' num2str(bb_time(ii))],Ftype,'-r0')
			end
		end
	end
	
end


if plotUIUC_CipPip
% 	for ix=1:length(sprlNames)
	for ix=1
		uiuc_cipConc = uiuc_cip_concAvg.(sprlNames{ix});
		uiuc_pipConc = uiuc_pip_concAvg.(sprlNames{ix});
		uiuc_time = insec2hhmmss(uiuc_cip_timeSecsAvg.(sprlNames{ix}));
		
% 		for ii = 1:size(uiuc_cipConc,1)
		for ii = 1:15
			if saveFigs && noDisp
				figure('visible','off','Position', [10,10,1200,700]);
			else
				figure('Position', [10,10,1200,700]);
			end
			
			stairs(uiuc_cip_binMin,uiuc_cipConc(ii,:),'r','LineWidth', 1.5);
			hold on
			stairs(uiuc_pip_binMin,uiuc_pipConc(ii,:),'b','LineWidth', 1.5);
			title([flight ' - UIUC CIP & PIP - Spiral ' num2str(ix) ' - ' num2str(uiuc_time(ii)) ' - ' num2str(ii)],'FontSize',24);
			lgnd = legend('CIP','PIP');
			lgnd.FontSize = 18;
			xlabel('D [mm]');
			ylabel('N(D) [cm^{-4}]');
			set(gca,'yscale','log');
			set(gca,'YTick',[10e-10 10e-9 10e-8 10e-7 10e-6... 
				10e-5 10e-4 10e-3 10e-2 10e-1 10e0 10e1 10e2 10e3 10e4 10e5 10e6 10e7]);
			grid on
			set(gca,'YMinorGrid','on','YMinorTick','on');
			ax = ancestor(gca, 'axes');
			xRule = ax.XAxis;
			yRule = ax.YAxis;
			xRule.FontSize = 24;
			yRule.FontSize = 24;
			
			if saveFigs
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				print([saveDir '/Black-Comparisons/' flight '_uiuc-CIP-PIP_S' num2str(ix) '-' sprintf('%2.2d',ii) '-' num2str(uiuc_time(ii))],Ftype,'-r0')
			end
		end
	end
	
	
end

if plotBB_UIUC_Compare
	for ix=1:length(sprlNames)
% 	for ix=1
		uiuc_cipConc = nansum(uiuc_cip_concAvg.(sprlNames{ix}),1);
		uiuc_pipConc = nansum(uiuc_pip_concAvg.(sprlNames{ix}),1);
		bb_cipConc = nansum(bb_cip_concAvg.(sprlNames{ix}),1);
		bb_pipConc = nansum(bb_pip_concAvg.(sprlNames{ix}),1);
		
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		
		stairs(uiuc_cip_binMin,uiuc_cipConc,'r','LineWidth', 1.5);
		hold on
		stairs(uiuc_pip_binMin,uiuc_pipConc,'b','LineWidth', 1.5);
		stairs(bb_cip_binMin,bb_cipConc,'m','LineWidth', 1.5);
		stairs(bb_pip_binMin,bb_pipConc,'c','LineWidth', 1.5);
		title([flight ' - CIP & PIP BB/UIUC Comparison - Spiral ' num2str(ix)],'FontSize',24);
		lgnd = legend('UIUC CIP','UIUC PIP','BB CIP','BB PIP');
		lgnd.FontSize = 18;
		xlabel('D [mm]');
		ylabel('N(D) [cm^{-4}]');
		set(gca,'yscale','log');
		set(gca,'YTick',[10e-10 10e-9 10e-8 10e-7 10e-6...
			10e-5 10e-4 10e-3 10e-2 10e-1 10e0 10e1 10e2 10e3 10e4 10e5 10e6 10e7 10e8 10e9 10e10]);
		grid on
		set(gca,'YMinorGrid','on','YMinorTick','on');
		ax = ancestor(gca, 'axes');
		xRule = ax.XAxis;
		yRule = ax.YAxis;
		xRule.FontSize = 24;
		yRule.FontSize = 24;
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/Black-Comparisons/' flight '_compare-uiuc-bb_S' num2str(ix)],Ftype,'-r0')
		end
	end
end

if plotCommonT
	for ix=1:length(cTnames)
		commonT = cT.(cTnames{ix});
		uiucCIPtime = floor(uiuc_cip_timeSecsAvg.(cTnames{ix}));
		uiucPIPtime = floor(uiuc_pip_timeSecsAvg.(cTnames{ix}));
		bbCIPtime = bb_cip_timeSecsAvg.(cTnames{ix});
		bbPIPtime = bb_pip_timeSecsAvg.(cTnames{ix});
		
		for iz=1:length(commonT)
			uiucCIPix = find(uiucCIPtime == commonT(iz));
			uiucPIPix = find(uiucPIPtime == commonT(iz));
			bbCIPix = find(bbCIPtime == commonT(iz));
			bbPIPix = find(bbPIPtime == commonT(iz));
			
			% If there are more than one set of SDs for a given time, just use the first
			if length(bbCIPix) > 1
				bbCIPix = bbCIPix(1);
			end
			if length(bbPIPix) > 1
				bbPIPix = bbPIPix(1);
			end
			
			uiucCIPconc = uiuc_cip_concAvg.(cTnames{ix})(uiucCIPix,:);
			uiucPIPconc = uiuc_pip_concAvg.(cTnames{ix})(uiucPIPix,:);
			bbCIPconc = bb_cip_concAvg.(cTnames{ix})(bbCIPix,:);
			bbPIPconc = bb_pip_concAvg.(cTnames{ix})(bbPIPix,:);
			
			if saveFigs && noDisp
				figure('visible','off','Position', [10,10,1200,700]);
			else
				figure('Position', [10,10,1200,700]);
			end
			
			stairs(uiuc_cip_binMin,uiucCIPconc,'r','LineWidth', 1.5);
			hold on
			stairs(uiuc_pip_binMin,uiucPIPconc,'b','LineWidth', 1.5);
			stairs(bb_cip_binMin,bbCIPconc,'m','LineWidth', 1.5);
			stairs(bb_pip_binMin,bbPIPconc,'Color',[0.,0.808,0.82],'LineWidth', 1.5);
			title([flight ' - CIP & PIP BB/UIUC Comparison - Spiral ' num2str(cTnames{ix}(end)) ' - ' num2str(insec2hhmmss(commonT(iz)))],'FontSize',24);
			lgnd = legend('UIUC CIP','UIUC PIP','BB CIP','BB PIP');
			lgnd.FontSize = 18;
			xlabel('D [mm]');
			ylabel('N(D) [cm^{-4}]');
			set(gca,'yscale','log');
			set(gca,'YTick',[10e-10 10e-9 10e-8 10e-7 10e-6...
				10e-5 10e-4 10e-3 10e-2 10e-1 10e0 10e1 10e2 10e3 10e4 10e5 10e6 10e7 10e8 10e9 10e10]);
			grid on
			set(gca,'YMinorGrid','on','YMinorTick','on');
			ax = ancestor(gca, 'axes');
			xRule = ax.XAxis;
			yRule = ax.YAxis;
			xRule.FontSize = 24;
			yRule.FontSize = 24;
			
			if saveFigs
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				print([saveDir '/Black-Comparisons/' flight '_compare-uiuc-bb_S' num2str(cTnames{ix}(end)) '_' num2str(insec2hhmmss(commonT(iz)))],Ftype,'-r0')
			end
		end
	end
end

if plotSingleT
	pltTime = 23983;
	sprlNm = 'sprl8';
	
	doCompare		= 0;
	doBBfracTotal	= 1;
	
	uiucCIPtime = floor(uiuc_cip_timeSecsAvg.(sprlNm));
	uiucPIPtime = floor(uiuc_pip_timeSecsAvg.(sprlNm));
	bbCIPtime = bb_cip_timeSecsAvg.(sprlNm);
	bbPIPtime = bb_pip_timeSecsAvg.(sprlNm);
	
	uiucCIPix = find(uiucCIPtime == pltTime);
	uiucPIPix = find(uiucPIPtime == pltTime);
	bbCIPix = find(bbCIPtime == pltTime);
	bbPIPix = find(bbPIPtime == pltTime);
	
	% If there are more than one set of SDs for a given time, just use the first
	if length(bbCIPix) > 1
		bbCIPix = bbCIPix(1);
	end
	if length(bbPIPix) > 1
		bbPIPix = bbPIPix(1);
	end
	
	uiucCIPconc = uiuc_cip_concAvg.(sprlNm)(uiucCIPix,:);
	uiucPIPconc = uiuc_pip_concAvg.(sprlNm)(uiucPIPix,:);
	bbCIPconc = bb_cip_concAvg.(sprlNm)(bbCIPix,:);
	bbPIPconc = bb_pip_concAvg.(sprlNm)(bbPIPix,:); % Center-in
	bbCIPconcT = bb_cip_concTAvg.(sprlNm)(bbCIPix,:);
	bbCIPconcF = bb_cip_concFAvg.(sprlNm)(bbCIPix,:);
	bbPIPconcT = bb_pip_concTAvg.(sprlNm)(bbPIPix,:);
	bbPIPconcF = bb_pip_concFAvg.(sprlNm)(bbPIPix,:);
	
	if doCompare
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		
		stairs(uiuc_cip_binMin,uiucCIPconc,'r','LineWidth', 1.5);
		hold on
		stairs(uiuc_pip_binMin,uiucPIPconc,'b','LineWidth', 1.5);
		stairs(bb_cip_binMin,bbCIPconc,'m','LineWidth', 1.5);
		stairs(bb_pip_binMin,bbPIPconc,'Color',[0.,0.808,0.82],'LineWidth', 1.5);
		title([flight ' - CIP & PIP BB/UIUC Comparison - Spiral ' num2str(sprlNm(end)) ' - ' num2str(insec2hhmmss(pltTime))],'FontSize',24);
		lgnd = legend('UIUC CIP','UIUC PIP','BB CIP','BB PIP');
		lgnd.FontSize = 18;
		xlabel('D [mm]');
		ylabel('N(D) [cm^{-4}]');
		set(gca,'yscale','log');
		set(gca,'YTick',[10e-10 10e-9 10e-8 10e-7 10e-6...
			10e-5 10e-4 10e-3 10e-2 10e-1 10e0 10e1 10e2 10e3 10e4 10e5 10e6 10e7 10e8 10e9 10e10]);
		grid on
		set(gca,'YMinorGrid','on','YMinorTick','on');
		ax = ancestor(gca, 'axes');
		xRule = ax.XAxis;
		yRule = ax.YAxis;
		xRule.FontSize = 24;
		yRule.FontSize = 24;
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/Black-Comparisons/' flight '_compare-uiuc-bb_S' num2str(sprlNm(end)) '_' num2str(insec2hhmmss(pltTime))],Ftype,'-r0')
		end
	end
	
	if doBBfracTotal
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		
		stairs(bb_cip_binMin,bbCIPconcT,'b','LineWidth', 1.5);
		hold on
		stairs(bb_cip_binMin,bbCIPconcF,'Color',[0.,0.808,0.82],'LineWidth', 1.5);
		stairs(bb_pip_binMin,bbPIPconcT,'r','LineWidth', 1.5);
		stairs(bb_pip_binMin,bbPIPconcF,'m','LineWidth', 1.5);
		title([flight ' - BB Total & Fractional N(D) Comparison - Spiral ' num2str(sprlNm(end)) ' - ' num2str(insec2hhmmss(pltTime))],'FontSize',24);
		lgnd = legend('BB CIP Total','BB CIP Partial','BB PIP Total','BB PIP Partial');
		lgnd.FontSize = 18;
		xlabel('D [mm]');
		ylabel('N(D) [cm^{-4}]');
		set(gca,'yscale','log');
		set(gca,'YTick',[10e-10 10e-9 10e-8 10e-7 10e-6...
			10e-5 10e-4 10e-3 10e-2 10e-1 10e0 10e1 10e2 10e3 10e4 10e5 10e6 10e7 10e8 10e9 10e10]);
		grid on
		set(gca,'YMinorGrid','on','YMinorTick','on');
		ax = ancestor(gca, 'axes');
		xRule = ax.XAxis;
		yRule = ax.YAxis;
		xRule.FontSize = 24;
		yRule.FontSize = 24;
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/Black-Comparisons/' flight '_compare-bb-Total-Partial_S' num2str(sprlNm(end)) '_' num2str(insec2hhmmss(pltTime))],Ftype,'-r0')
		end
	end
end