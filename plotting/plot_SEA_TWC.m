close all;clearvars;

flight = '20150617';

probe = 'CIP'; 
avgTime = 10; % 5 or 10 will work here since we're using 1-sec data in this script


plotTWCtempTime	= 1;
plotNDtwc		= 1;


saveFigs    = 1;
noDisp      = 0;
% Ftype		= '-dpdf';
Ftype		= '-dpng';


% Range of plotted variables over subset/all flights/spirals
% Useful if it is desired to directly compare plots
% TWCrangeAll = [-0.1 3.2];
TWCrangeAll = [-3 3.2];
tempRangeAll = [-20 21];
NDLogLim = [-4 2];

savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';
dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

startT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'startT');
endT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'endT');

sDistFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.' probe '.' num2str(avgTime) 'secAvg.mat'];

SEAfile = [dataPath 'mp-data/' flight '/' flight '_SEA_TWC.txt'];

%% Load data
loadSEAcsv % Only requires that SEAfile variable be set above -- Outputs SEA_Time (secs), and SEA_TWC

load(sDistFile,'time_secsFL_orig','tempC_orig','conc_minR_orig','time_secs_orig','bin_mid');

sprlNames = fieldnames(time_secsFL_orig);


%% Create directories to save plots in if they don't already exist
if saveFigs
    saveDir = [savePath flight];
	if (exist(saveDir, 'dir') ~= 7)
        mkdir(saveDir)
	end

	if (plotTWCtempTime && exist([saveDir '/SEA-TWC-Temp-Time'], 'dir') ~= 7)
		mkdir([saveDir '/SEA-TWC-Temp-Time'])
	end
	if (plotNDtwc && exist([saveDir '/' probe '-ND-Time-TWC_1s'], 'dir') ~= 7)
		mkdir([saveDir '/' probe '-ND-Time-TWC_1s'])
	end

end

%% Plotting

if plotTWCtempTime

	for ix = 1:length(startT)
		if saveFigs && noDisp
            figure('visible','off','Position', [10,10,1500,1000]);
		else
            figure('Position', [10,10,1500,1000]);
		end

		sprlIx_SEA = find(SEA_Time >= startT(ix) & SEA_Time < endT(ix));
		timeSecs_SEA = SEA_Time(sprlIx_SEA);
		twc_SEA = SEA_TWC(sprlIx_SEA);
		
		yyaxis left
		plot(timeSecs_SEA/3600/24,twc_SEA,'b-')
		xlabel('Time (UTC)');
		datetick('x','HH:MM:SS');
        ylabel('TWC (g m^{-3})');
%         ylim(TWCrangeAll);
		
		yyaxis right
		plot(time_secsFL_orig.(sprlNames{ix})/3600/24,tempC_orig.(sprlNames{ix}),'r-')
		yLab = sprintf('Temperature (%cC)',char(176));
		ylabel(yLab);
        ylim(tempRangeAll);
		set(gca,'YDir','reverse');

        title([flight ' - Spiral ' num2str(ix) ' - SEA TWC & Temp vs. Time']);
		grid
        set(gca,'XTickLabelRotation',45)
        set(findall(gcf,'-property','FontSize'),'FontSize',28)
		if saveFigs
            set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/SEA-TWC-Temp-Time/' flight '_SEAtwcTempTime_S' num2str(ix)],Ftype,'-r0')
		end
	end
end

if plotNDtwc
    for ix = 1:length(startT)

		conc_minR = conc_minR_orig.(sprlNames{ix});
		time_secs = time_secs_orig.(sprlNames{ix});
		
		sprlIx_SEA = find(SEA_Time >= startT(ix) & SEA_Time < endT(ix));
		timeSecs_SEA = SEA_Time(sprlIx_SEA);
		twc_SEA = SEA_TWC(sprlIx_SEA);
		goodTWC = twc_SEA(~isnan(twc_SEA));
		goodTWCtimeSec = timeSecs_SEA(~isnan(twc_SEA));
		
        if saveFigs && noDisp
            figure('visible','off','Position', [10,10,1920,700]);
        else
            figure('Position', [10,10,1920,700]);
		end
		
		set(gcf,'defaultAxesColorOrder',[0 0 1; 0 0 0]); % Blue and black y-axes will be black
		
		yyaxis left
		contourf(time_secs/24/3600,bin_mid,log10(conc_minR'),NDLogLim(1):0.1:NDLogLim(2),'LineColor','none');

		title([flight ' - Spiral ' num2str(ix) ' - ' probe ' - 1s Avg']);

		ylabel('D [mm]');
		datetick('x','HH:MM:SS');
		set(gca,'XMinorTick','on','YMinorTick','on','XTickLabelRotation',45);
		colormap(jetmod); %Uses modified 'jet' colormap
		c=colorbar('southoutside');
		ylabel(c,'log_{10}N(D) [cm^{-4}]');
		set(findall(gcf,'-property','FontSize'),'FontSize',28)
		set(gca, 'CLim', NDLogLim);
		hold on
		
		yyaxis right
		plot(goodTWCtimeSec/3600/24,goodTWC,'k-','LineWidth',2)
		plot(xlim,[0 0],'b--');
        ylabel('TWC (g m^{-3})');
        ylim(TWCrangeAll);

        
        if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/' probe '-ND-Time-TWC_1s/' flight '_' probe '_NDtimeTWC_1s_S' num2str(ix)],Ftype,'-r0')
        end
    end
end