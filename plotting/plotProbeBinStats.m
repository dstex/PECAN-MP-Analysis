close all; clearvars;

flight = '20150709';
probe = 'PIP';

normalize = 1;

timeStep = -999; % set to -999 to use whole spiral

saveFigs = 1;
noDisp = 1;
ftype = '-dpdf';
% ftype = '-dpng';

savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';
dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

proc2File = [dataPath 'mp-data/' flight '/proc2.' flight '.' probe '.cdf'];

if saveFigs
	saveDir = [savePath flight];
	if (exist(saveDir, 'dir') ~= 7)
		mkdir(saveDir)
	end
	if (exist([saveDir '/' probe '-DiodeShadows'], 'dir') ~= 7 )
		mkdir([saveDir '/' probe '-DiodeShadows'])
	end
end

startT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'startT');
endT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'endT');


binStats = nc_varget(proc2File,'image_bin_stats');
timeSecs = nc_varget(proc2File,'Time_in_seconds');
timehhmmss = insec2hhmmss(timeSecs);

edges = 1:64;



for ix = 1:length(startT)
% for ix = 2
	sprlIx = find(timeSecs >= startT(ix) & timeSecs <= endT(ix));
	binStats_sprl = binStats(:,sprlIx);
	timeSecs_sprl = timeSecs(sprlIx);
	timehhmmss_sprl = timehhmmss(sprlIx);
	
	if timeStep == -999
		disp(['Now plotting stats for ' num2str(timehhmmss_sprl(1)) '-' num2str(timehhmmss_sprl(end))]);

		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,800]);
		else
			figure('Position', [10,10,1200,800]);
		end

		if normalize
			h = bar(edges,normc(sum(binStats_sprl,2)),0.95,'FaceColor','b');
			ylim([0,0.25])
		else
			h = bar(edges,sum(binStats_sprl,2),0.95,'FaceColor','b');
		end
		
		xlabel('Diode Number','FontSize',24);
		ylabel('Shadow Counts','FontSize',24);
		ax = ancestor(h, 'axes');
		yRule = ax.YAxis;
		xRule = ax.XAxis;
		title([flight ' - ' probe ' Diode Shadow Counts Spiral ' num2str(ix) ' - ' num2str(timehhmmss_sprl(1)) '-' num2str(timehhmmss_sprl(end))],'FontSize',28);
		set(gca,'XTick',1:64,'XTickLabelRotation',90);
		yRule.FontSize = 24;


		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/' probe '-DiodeShadows/' flight '_' probe '_diodeShadows_S' num2str(ix)],ftype,'-r0')
		end
		
	else
		timeBeg = floor(timeSecs_sprl(1));
		timeEnd = timeBeg + timeStep;
		
		while timeEnd <= timeSecs_sprl(end)
			
			plotIx = find(timeSecs_sprl >= timeBeg & timeSecs_sprl < timeEnd);
			
			if length(plotIx) >= 2
				disp(['Now plotting stats for ' num2str(timehhmmss_sprl(plotIx(1))) '-' num2str(timehhmmss_sprl(plotIx(end)))]);
				
				if saveFigs && noDisp
					figure('visible','off','Position', [10,10,1200,800]);
				else
					figure('Position', [10,10,1200,800]);
				end
				
				if normalize
					h = bar(edges,normc(sum(binStats_sprl(:,plotIx),2)),0.95,'FaceColor','b');
					ylim([0,0.25])
				else
					h = bar(edges,sum(binStats_sprl(:,plotIx),2),0.95,'FaceColor','b');
				end
				xlabel('Diode Number','FontSize',24);
				ylabel('Shadow Counts','FontSize',24);
				ax = ancestor(h, 'axes');
				yRule = ax.YAxis;
				xRule = ax.XAxis;
				title([flight ' - ' probe ' Diode Shadow Counts Spiral ' num2str(ix) ' - ' num2str(timehhmmss_sprl(plotIx(1))) '-' num2str(timehhmmss_sprl(plotIx(end)))],'FontSize',28);
				set(gca,'XTick',1:64,'XTickLabelRotation',90);
				yRule.FontSize = 24;
				
				
				if saveFigs
					set(gcf,'Units','Inches');
					pos = get(gcf,'Position');
					set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
					print([saveDir '/' probe '-DiodeShadows/' flight '_' probe '_diodeShadows_' num2str(timeStep) 's_S' num2str(ix) '_' num2str(timehhmmss_sprl(plotIx(end)))],ftype,'-r0')
				end
			end
			
			timeBeg = timeEnd;
			timeEnd = timeBeg + timeStep;
		end
	end
	
	
end