% This code was used to produce some fit sensitivity and comparison plots
% The script as-is (12/20/17) will not run without some additions as it
% was meant to be used with a premade mat file containing comparison data

% Plot comparisons between IGF and Exp
plotCmpSprl		= 0;
plotCmpTempBins	= 0;
plotCmpEvryTStp	= 0;

% Plots for sensitivity tests
plotAllGood		= 0;
plotAllSkip		= 0;
plotAllNegLmda	= 0;


%% Plots for IGF and exponential fit comparison
if plotCmpSprl
	for ix=loopVctr
		cipConcSprl = nanmean(cip_concMinR.(sprlNames{ix}),1);
		
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		
		stairs(cip_binMin*10,cipConcSprl,'Color',[0.75 0.75 0.75],'LineWidth',2,'DisplayName','Obs - All');
		hold on
		stairs(cip_binMin(cipIncld)*10,cipConcSprl(cipIncld),'k','LineWidth',2,'DisplayName','Obs - Incld');
		plot(cipExt_binMid(cipIncld(1):end)*10,cipConc_ext_igfWhl.(sprlNames{ix})(cipIncld(1):end),'b','LineWidth',2,'DisplayName','IGF');
		plot(cipExt_binMid(cipIncld(1):end)*10,cipConc_ext_expWhl.(sprlNames{ix})(cipIncld(1):end),'m','LineWidth',2,'DisplayName','L-M Exp');

		title(sprintf('%s - CIP N(D) Fits - Spiral %d - %d sec avg%s',flight,ix,avgTime,titleIdStr),'FontSize',24);
		
		lgnd = legend('show');
		lgnd.FontSize = 18;

		xlabel('D [mm]');
		ylabel('N(D) [cm^{-4}]');
		set(gca,'yscale','log','XScale','log');
		grid on
		set(gca,'YMinorGrid','on','YMinorTick','on');
		set(gca,'XLim',[0.1 3])
		set(gca,'YLimMode','auto')
		
		ax = ancestor(gca, 'axes');
		xRule = ax.XAxis;
		yRule = ax.YAxis;
		xRule.FontSize = 24;
		yRule.FontSize = 24;
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir flight '_FitComparison-CIP_S' num2str(ix) '_Whole' fileIdStr],Ftype,'-r0')
		end
		
	end
end

if plotCmpTempBins
	for ix=loopVctr
		fprintf('Plotting Spiral %d\n',ix)
		cipConc = cip_concMinR.(sprlNames{ix});
		tempBins = tempBinsAll.(sprlNames{ix});
		tmpC = tempCsprl.(sprlNames{ix});
		for tmp = 1:length(tempBins)
			if all(isnan(cip_igf_nmlTb.(sprlNames{ix})(tmp)))
				fprintf('\tSkipping plots of %d - no valid fits\n',tempBins(tmp))
				continue
			else
				fprintf('\tPlotting %d\n',tempBins(tmp))
			end
			
			if saveFigs && noDisp
				figure('visible','off','Position', [10,10,1200,700]);
			else
				figure('Position', [10,10,1200,700]);
			end
			
			stairs(cip_binMin*10,nanmean(cipConc(tmpC == tempBins(tmp),:),1),'Color',[0.75 0.75 0.75],'LineWidth',2,'DisplayName','Obs - All');
			hold on
			stairs(cip_binMin(cipIncld)*10,nanmean(cipConc(tmpC == tempBins(tmp),cipIncld),1),'k','LineWidth',2,'DisplayName','Obs - Incld');
			plot(cipExt_binMid(cipIncld(1):end)*10,cipConc_ext_igfTb.(sprlNames{ix})(tmp,cipIncld(1):end),'b','LineWidth',2,'DisplayName','IGF');
			plot(cipExt_binMid(cipIncld(1):end)*10,cipConc_ext_expTb.(sprlNames{ix})(tmp,cipIncld(1):end),'m','LineWidth',2,'DisplayName','L-M Exp');

			title(sprintf('%s - CIP N(D) Fits - Spiral %d - %s',flight,ix,sprintf('%d%cC avg Temp',tempBins(tmp),char(176))),'FontSize',24);
			
			lgnd = legend('show');
			lgnd.FontSize = 18;

			xlabel('D [mm]');
			ylabel('N(D) [cm^{-4}]');
			set(gca,'yscale','log','XScale','log');
			grid on
			set(gca,'YMinorGrid','on','YMinorTick','on');
			set(gca,'XLim',[0.1 3])
			set(gca,'YLimMode','auto')
			
			ax = ancestor(gca, 'axes');
			xRule = ax.XAxis;
			yRule = ax.YAxis;
			xRule.FontSize = 24;
			yRule.FontSize = 24;
			
			if saveFigs
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				print([saveDir flight '_FitComparison-CIP_S' num2str(ix) '_' num2str(tempBins(tmp)) 'degC' fileIdStr],Ftype,'-r0')
			end
		end
		
	end
end

if plotCmpEvryTStp
	for ix=loopVctr
		fprintf('Plotting Spiral %d\n',ix)
		cipConc = cip_concMinR.(sprlNames{ix});
		cipTsec = cip_timeSecs.(sprlNames{ix});
		for time = 1:size(cipConc,1)
			if all(isnan(cip_igf_nml.(sprlNames{ix})(time)))
				fprintf('\tSkipping plots of time %d/%d - no valid fits\n',time,size(cipConc,1))
				continue
			else
				fprintf('\tPlotting time %d/%d\n',time,size(cipConc,1))
			end
			
			crntTempC = tempC.(sprlNames{ix})(time);
			
			pltT.allSec = cipTsec(time);
			pltT.hour = floor(pltT.allSec/60^2);
			pltT.minutes = floor(mod((pltT.allSec/60), 60));
			pltT.seconds = floor(mod(pltT.allSec,60));
			pltT.str = sprintf('%02d:%02d:%02d', pltT.hour, pltT.minutes, pltT.seconds);
			
			if saveFigs && noDisp
				figure('visible','off','Position', [10,10,1200,700]);
			else
				figure('Position', [10,10,1200,700]);
			end
			
			stairs(cip_binMin*10,cipConc(time,:),'Color',[0.75 0.75 0.75],'LineWidth',2,'DisplayName','Obs - All');
			hold on
			stairs(cip_binMin(cipIncld)*10,cipConc(time,cipIncld),'k','LineWidth',2,'DisplayName','Obs - Incld');
			plot(cipExt_binMid(cipIncld(1):end)*10,cipConc_ext_igf.(sprlNames{ix})(time,cipIncld(1):end),'b','LineWidth',2,'DisplayName','IGF');
			plot(cipExt_binMid(cipIncld(1):end)*10,cipConc_ext_exp.(sprlNames{ix})(time,cipIncld(1):end),'m','LineWidth',2,'DisplayName','L-M Exp');

			title({sprintf('%s - CIP N(D) Fits - Spiral %d - #%d - %s - %s',flight,ix,time,pltT.str,titleIdStr), sprintf('%.2f%cC avg Temp',crntTempC,char(176))},'FontSize',24);
			
			lgnd = legend('show');
			lgnd.FontSize = 18;

			xlabel('D [mm]');
			ylabel('N(D) [cm^{-4}]');
			set(gca,'yscale','log','XScale','log');
			grid on
			set(gca,'YMinorGrid','on','YMinorTick','on');
			set(gca,'XLim',[0.1 3])
			set(gca,'YLimMode','auto')
			
			ax = ancestor(gca, 'axes');
			xRule = ax.XAxis;
			yRule = ax.YAxis;
			xRule.FontSize = 24;
			yRule.FontSize = 24;
			
			if saveFigs
				set(gcf,'Units','Inches');
				pos = get(gcf,'Position');
				set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
				print([saveDir flight '_FitComparison-CIP_S' num2str(ix) '_' sprintf('%2.2d',time) fileIdStr],Ftype,'-r0')
			end
		end
		
	end
end

%% Plots for sensitivity tests
if plotAllGood
	for ix=loopVctr
		% Get N(D) and M(D) for all non-skipped points
		cipConcGood = cip_concMinR.(sprlNames{ix})(setdiff(1:size(cip_concMinR.(sprlNames{ix}),1),fitSkipIx.(sprlNames{ix})),:);
		cipMassGood = cip_massTWC.(sprlNames{ix})(setdiff(1:size(cip_massTWC.(sprlNames{ix}),1),fitSkipIx.(sprlNames{ix})),:);
		
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		hold on
		
		for iGd = 1:size(cipConcGood,1)
			plot(cip_binMin*10,cipConcGood(iGd,:),'LineWidth',1.5);
		end
		
		title(sprintf('%s - CIP N(D) Good - Spiral %d - %d sec avg%s',flight,ix,avgTime,titleIdStr),'FontSize',24);
		
		xlabel('D [mm]');
		ylabel('N(D) [cm^{-4}]');
		xlim([0.1 2])
		set(gca,'yscale','log','XScale','log');
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
			print([saveDir flight '_ND-GoodAll-CIP_S' num2str(ix) fileIdStr],Ftype,'-r0')
		end
		
		
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		hold on
		
		for iGd = 1:size(cipMassGood,1)
			plot(cip_binMin*10,cipMassGood(iGd,:),'LineWidth',1.5);
		end
		
		title(sprintf('%s - CIP M(D) Good - Spiral %d - %d sec avg%s',flight,ix,avgTime,titleIdStr),'FontSize',24);
		
		xlabel('D [mm]');
		ylabel('M(D) [g cm^{-4}]');
		xlim([0.1 2])
		set(gca,'yscale','log','XScale','log');
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
			print([saveDir flight '_MD-GoodAll-CIP_S' num2str(ix) fileIdStr],Ftype,'-r0')
		end
		
	end
end

if plotAllSkip
	for ix=loopVctr
		if isempty(fitSkipIx.(sprlNames{ix}))
			fprintf('\nNo fits were skipped for Spiral %d',ix)
			continue
		end
		cipConcSkip = cip_concMinR.(sprlNames{ix})(fitSkipIx.(sprlNames{ix}),:);
		cipMassSkip = cip_massTWC.(sprlNames{ix})(fitSkipIx.(sprlNames{ix}),:);
		
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		hold on
		
		for iSkp = 1:size(cipConcSkip,1)
			plot(cip_binMin*10,cipConcSkip(iSkp,:),'LineWidth',1.5);
		end
		
		title({sprintf('%s - CIP N(D) Fit Skips - Spiral %d - %d sec avg%s',flight,ix,avgTime,titleIdStr),...
			sprintf('%d/%d SDs Skipped',size(cipConcSkip,1),size(cip_concMinR.(sprlNames{ix}),1))},'FontSize',24);
		
		xlabel('D [mm]');
		ylabel('N(D) [cm^{-4}]');
		xlim([0.1 2])
		set(gca,'yscale','log','XScale','log');
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
			print([saveDir flight '_ND-FitSkipsAll-CIP_S' num2str(ix) fileIdStr],Ftype,'-r0')
		end
		
		
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		hold on
		
		for iSkp = 1:size(cipMassSkip,1)
			plot(cip_binMin*10,cipMassSkip(iSkp,:),'LineWidth',1.5);
		end
		
		title({sprintf('%s - CIP M(D) Fit Skips - Spiral %d - %d sec avg%s',flight,ix,avgTime,titleIdStr),...
			sprintf('%d/%d SDs Skipped',size(cipMassSkip,1),size(cip_concMinR.(sprlNames{ix}),1))},'FontSize',24);
		
		xlabel('D [mm]');
		ylabel('M(D) [g cm^{-4}]');
		xlim([0.1 2])
		set(gca,'yscale','log','XScale','log');
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
			print([saveDir flight '_MD-FitSkipsAll-CIP_S' num2str(ix) fileIdStr],Ftype,'-r0')
		end
		
	end
end

if plotAllNegLmda
	for ix=loopVctr
		if isempty(negLmdaIx.(sprlNames{ix}))
			fprintf('No fits with lambda < 0 for Spiral %d\n',ix)
			continue
		end
		cipConcNL = cip_concMinR.(sprlNames{ix})(negLmdaIx.(sprlNames{ix}),:);
		cipMassNL = cip_massTWC.(sprlNames{ix})(negLmdaIx.(sprlNames{ix}),:);
		
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		hold on
		
		for iNL = 1:size(cipConcNL,1)
			plot(cip_binMin*10,cipConcNL(iNL,:),'LineWidth',1.5);
		end
		lmdaStr = '\lambda';
		title({sprintf('%s - CIP N(D) w/%s<0 - Spiral %d - %d sec avg%s',flight,lmdaStr,ix,avgTime,titleIdStr),...
			sprintf('%d/%d SDs w/%s<0',size(cipConcNL,1),size(cip_concMinR.(sprlNames{ix}),1),lmdaStr)},'FontSize',24);
		
		xlabel('D [mm]');
		ylabel('N(D) [cm^{-4}]');
		xlim([0.1 2])
		set(gca,'yscale','log','XScale','log');
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
			print([saveDir flight '_ND-NegLmdaAll-CIP_S' num2str(ix) fileIdStr],Ftype,'-r0')
		end
		
		
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		hold on
		
		for iNL = 1:size(cipMassNL,1)
			plot(cip_binMin*10,cipMassNL(iNL,:),'LineWidth',1.5);
		end
		
		title({sprintf('%s - CIP M(D) w/%s<0 - Spiral %d - %d sec avg%s',flight,lmdaStr,ix,avgTime,titleIdStr),...
			sprintf('%d/%d SDs w/%s<0',size(cipMassNL,1),size(cip_concMinR.(sprlNames{ix}),1),lmdaStr)},'FontSize',24);
		
		xlabel('D [mm]');
		ylabel('M(D) [g cm^{-4}]');
		xlim([0.1 2])
		set(gca,'yscale','log','XScale','log');
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
			print([saveDir flight '_MD-NegLmdaAll-CIP_S' num2str(ix) fileIdStr],Ftype,'-r0')
		end
		
	end
end