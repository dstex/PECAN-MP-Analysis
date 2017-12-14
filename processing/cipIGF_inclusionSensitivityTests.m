% load 20150706_cipIGF_gte7-17-25-29.mat
close all;


saveFigs = 1;

saveDir = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/20150706/CIP-Fitting';
if saveFigs
	if (exist([saveDir '/Sensitivity'], 'dir') ~= 7)
		mkdir([saveDir '/Sensitivity'])
	end
end

sprlNames = fieldnames(cipConc_hybrid_igf_7);

for ix = 1:8
	%{ 	
	negLmbda_7 = find(cip_igf_nml_7.(sprlNames{ix})(:,3) < 0);
	cipConc_ext_igf_7.(sprlNames{ix})(negLmbda_7,:) = NaN;
	cipConc_hybrid_igf_7.(sprlNames{ix})(negLmbda_7,length(cip_binMid)+1:end) = NaN;
	cipMass_ext_igf_7.(sprlNames{ix})(negLmbda_7,:) = NaN;
	cipMass_hybrid_igf_7.(sprlNames{ix})(negLmbda_7,length(cip_binMid)+1:end) = NaN;
	
	negLmbda_17 = find(cip_igf_nml_17.(sprlNames{ix})(:,3) < 0);
	cipConc_ext_igf_17.(sprlNames{ix})(negLmbda_17,:) = NaN;
	cipConc_hybrid_igf_17.(sprlNames{ix})(negLmbda_17,length(cip_binMid)+1:end) = NaN;
	cipMass_ext_igf_17.(sprlNames{ix})(negLmbda_17,:) = NaN;
	cipMass_hybrid_igf_17.(sprlNames{ix})(negLmbda_17,length(cip_binMid)+1:end) = NaN;
	
	negLmbda_25 = find(cip_igf_nml_25.(sprlNames{ix})(:,3) < 0);
	cipConc_ext_igf_25.(sprlNames{ix})(negLmbda_25,:) = NaN;
	cipConc_hybrid_igf_25.(sprlNames{ix})(negLmbda_25,length(cip_binMid)+1:end) = NaN;
	cipMass_ext_igf_25.(sprlNames{ix})(negLmbda_25,:) = NaN;
	cipMass_hybrid_igf_25.(sprlNames{ix})(negLmbda_25,length(cip_binMid)+1:end) = NaN;
	
	negLmbda_29 = find(cip_igf_nml_29.(sprlNames{ix})(:,3) < 0);
	cipConc_ext_igf_29.(sprlNames{ix})(negLmbda_29,:) = NaN;
	cipConc_hybrid_igf_29.(sprlNames{ix})(negLmbda_29,length(cip_binMid)+1:end) = NaN;
	cipMass_ext_igf_29.(sprlNames{ix})(negLmbda_29,:) = NaN;
	cipMass_hybrid_igf_29.(sprlNames{ix})(negLmbda_29,length(cip_binMid)+1:end) = NaN;

	
	cipMass_hybrid_7tmp = cipMass_hybrid_igf_7.(sprlNames{ix}).*(cipExt_binwidth'); %convert to g/cm3
	cipTWC_hybrid_7.(sprlNames{ix}) = nansum(cipMass_hybrid_7tmp,2)*1e6; % Sum masses over all bins and convert to g/m3
	cipMass_hybrid_17tmp = cipMass_hybrid_igf_17.(sprlNames{ix}).*(cipExt_binwidth'); 
	cipTWC_hybrid_17.(sprlNames{ix}) = nansum(cipMass_hybrid_17tmp,2)*1e6;
	cipMass_hybrid_25tmp = cipMass_hybrid_igf_25.(sprlNames{ix}).*(cipExt_binwidth'); 
	cipTWC_hybrid_25.(sprlNames{ix}) = nansum(cipMass_hybrid_25tmp,2)*1e6;
	cipMass_hybrid_29tmp = cipMass_hybrid_igf_29.(sprlNames{ix}).*(cipExt_binwidth'); 
	cipTWC_hybrid_29.(sprlNames{ix}) = nansum(cipMass_hybrid_29tmp,2)*1e6;
	%}
	
% 	cipTWC_hybrid_total_7.(sprlNames{ix}) = nansum(cipTWC_hybrid_7.(sprlNames{ix}));
% 	cipTWC_hybrid_total_17.(sprlNames{ix}) = nansum(cipTWC_hybrid_17.(sprlNames{ix}));
% 	cipTWC_hybrid_total_25.(sprlNames{ix}) = nansum(cipTWC_hybrid_25.(sprlNames{ix}));
% 	cipTWC_hybrid_total_29.(sprlNames{ix}) = nansum(cipTWC_hybrid_29.(sprlNames{ix}));
	
% 	fprintf('\nSpiral %d',ix)
% 	fprintf('\n\tTWC_total_7 - TWC_total_17 = %.6f ',cipTWC_hybrid_total_7.(sprlNames{ix})-cipTWC_hybrid_total_17.(sprlNames{ix}))
% 	fprintf('\n\tTWC_total_7 - TWC_total_25 = %.6f',cipTWC_hybrid_total_7.(sprlNames{ix})-cipTWC_hybrid_total_25.(sprlNames{ix}))
% 	fprintf('\n\tTWC_total_7 - TWC_total_29 = %.6f',cipTWC_hybrid_total_7.(sprlNames{ix})-cipTWC_hybrid_total_29.(sprlNames{ix}))
	

	figure('Position', [10,10,1200,700]);
	plot(abs(cipTWC_hybrid_7.(sprlNames{ix})-cipTWC_hybrid_17.(sprlNames{ix})),'d','MarkerSize',15,'DisplayName','17TWC','LineWidth',2);
	hold on
	plot(abs(cipTWC_hybrid_7.(sprlNames{ix})-cipTWC_hybrid_25.(sprlNames{ix})),'s','MarkerSize',15,'DisplayName','25TWC','LineWidth',2);
	plot(abs(cipTWC_hybrid_7.(sprlNames{ix})-cipTWC_hybrid_29.(sprlNames{ix})),'k.','MarkerSize',15,'DisplayName','29TWC','LineWidth',2);
	
	
	
	title({sprintf('20150706 - Spiral %d - CIP IGF-Extended TWC - Diff. in # req. PSD values',ix),...
		sprintf('FSR-7=%.2f%%\tFSR-17=%.2f%%\tFSR-25=%.2f%%\tFSR-29=%.2f%%',fitSkipRatio_7.(sprlNames{ix})*100,...
		fitSkipRatio_17.(sprlNames{ix})*100,fitSkipRatio_25.(sprlNames{ix})*100,fitSkipRatio_29.(sprlNames{ix})*100),...
		sprintf('TWCWhl7-17 = %.6f\tTWCWhl7-25 = %.6f\tTWCWhl7-29 = %.6f g m^{-3}',cipTWC_hybrid_total_7.(sprlNames{ix})-cipTWC_hybrid_total_17.(sprlNames{ix}),...
		cipTWC_hybrid_total_7.(sprlNames{ix})-cipTWC_hybrid_total_25.(sprlNames{ix}),...
		cipTWC_hybrid_total_7.(sprlNames{ix})-cipTWC_hybrid_total_29.(sprlNames{ix}))},'FontSize',20);
	
	lgnd = legend('show');
	lgnd.FontSize = 18;
	xlabel('Time Dimension');
	ylabel('TWC Loss (relative to 7TWC) [g m^{-3}]');
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
		print([saveDir '/Sensitivity/20150706_FitSens-CIP_S' num2str(ix)],'-dpng','-r0')
		save([saveDir '/Sensitivity/20150706_FitSens-CIP_S' num2str(ix) '.fig']);
	end

end