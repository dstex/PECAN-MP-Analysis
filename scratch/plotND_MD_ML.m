for ix = 16
	cipMass_fit = cipMass_ext_igf.(sprlNames{ix});
	cipConc_fit = cipConc_ext_igf.(sprlNames{ix});
	tempC = tempC_avg.(sprlNames{ix});
	
	figure('Position', [10,10,1200,700]);
	hold on
	colors = varycolor(length(25:41));
	icx=1;
	for iPsd = 25:41
		plot(cipExt_binMid*10, cipConc_fit(iPsd,:)','Color',colors(icx,:),'LineWidth', 2,'DisplayName',sprintf('%.2f%cC',tempC(iPsd),char(176)));
		icx = icx+1;
	end
	title(sprintf('%s - Spiral %d - CIP (Extended)',flight,ix));
	
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
	
	figure('Position', [10,10,1200,700]);
	hold on
	colors = varycolor(length(25:41));
	icx = 1;
	for iPsd = 25:41
		plot(cipExt_binMid*10, cipMass_fit(iPsd,:)','Color',colors(icx,:),'LineWidth', 2,'DisplayName',sprintf('%.2f%cC',tempC(iPsd),char(176)));
		icx = icx+1;
	end
	title(sprintf('%s - Spiral %d - CIP (Extended)',flight,ix));
	
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
end
