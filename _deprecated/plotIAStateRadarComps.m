% Plots geotiff files of NEXRAD composites obtained from Iowa State data archive
% Download tiffs using the following address, replacing the date and time appropriately, where
% time must be in 5-minute intervals
%	http://mesonet.agron.iastate.edu/request/gis/n0r2gtiff.php?dstr=YYYYMMDDHHMM

close all; clearvars;

saveFigs	= 1;
noDisp		= 1;

doICCP		= 1;
doStatic	= 0;
doAll		= 0;

doLarge = 0;
doZoom	= 0;

flight = '20150706-ICCP-6';

% Get the times we want radar composites for (must be in 5-min increments as
% this is all Iowa State has available).
% In general, I aimed for the time at the middle of each spiral
switch flight
	case '20150706'
		times = {'0335', '0430', '0445', '0550', '0605', '0630', '0645'};
		
	case '20150706-ICCP-1'
		flight = '20150706';
		times = {'0330', '0335', '0340', '0345'};
		startTiccp = [12595, 12600, 12900, 13200];
		endTiccp = [12600, 12900, 13200, 13500];
		staticT = 12930;
		staticImgT = '0335';
		
	case '20150706-ICCP-6'
		flight = '20150706';
		times = {'0620', '0625', '0630', '0635'};
		startTiccp = [22795, 22800, 23100, 23400];
		endTiccp = [22800, 23100, 23400, 23700];
		staticT = 22977;
		staticImgT = '0625';
		
	case '20150709'
		times = {'0235', '0250', '0305', '0315', '0340', '0355', '0410',...
			'0425', '0445', '0500', '0525', '0540', '0600', '0610', '0630', '0640'};
		
	case '20150709-ICCP-2'
		flight = '20150709';
		times = {'0240', '0245', '0250', '0255'};
		startTiccp = [9595, 9600, 9900, 10200];
		endTiccp = [9600, 9900, 10200, 10500];
		staticT = 9812;
		staticImgT = '0245';
		
	case '20150709-ICCP-5'
		flight = '20150709';
		times = {'0330', '0335', '0340', '0345'};
		startTiccp = [12547, 12600, 12900, 13200];
		endTiccp = [12600, 12900, 13200, 13500];
		staticT = 12565;
		staticImgT = '0330';
		
	case '20150709-ICCP-9'
		flight = '20150709';
		times = {'0435', '0440', '0445', '0450'};
		startTiccp = [16495, 16500, 16800, 17100];
		endTiccp = [16500, 16800, 17100, 17400];
		staticT = 17303;
		staticImgT = '0450';
		
	
		
end

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';
savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';

startT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'startT');
endT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'endT');

if saveFigs
    saveDir = [savePath flight];
	if (exist(saveDir, 'dir') ~= 7)
        mkdir(saveDir)
	end
	if (exist([saveDir '/RadComps-SprlLocs'], 'dir') ~= 7)
        mkdir([saveDir '/RadComps-SprlLocs'])
	end
end

FLfile = [dataPath 'FlightLevelData/Processed/' flight '_FltLvl_Processed.mat'];

importVars = {'time_secs_FL','lat','lon'};

tempLoad = load(FLfile,importVars{:});

time_secs = tempLoad.(importVars{1});
lat = tempLoad.(importVars{2});
lon = tempLoad.(importVars{3});


if doICCP
	latSprlSmall = [];
	lonSprlSmall = [];
	for ix=1:length(times)
		radComp = [dataPath 'IowaStateRadarComps/n0r_' flight times{ix} '.tif'];
		
		spiralIXsmall = find(time_secs >= startTiccp(ix) & time_secs <= endTiccp(ix));
		
		[~, sprlNum] = min(abs(endTiccp(2) - endT));
		spiralIX = find(time_secs >= startT(sprlNum) & time_secs <= endT(sprlNum)); % For map limits

		latSprlSmall = vertcat(latSprlSmall,lat(spiralIXsmall));
		lonSprlSmall = vertcat(lonSprlSmall,lon(spiralIXsmall));
		
		% Determine map limits
		latSprl = lat(spiralIX);
		lonSprl = lon(spiralIX);
		maxLat = max(latSprl);
		minLat = min(latSprl);
		maxLon = max(lonSprl);
		minLon = min(lonSprl);
		if doLarge
			mapLatLim = [minLat-2, maxLat+2];
			mapLonLim = [minLon-2, maxLon+2];
		elseif doZoom
			mapLatLim = [minLat-0.3, maxLat+0.3];
			mapLonLim = [minLon-0.3, maxLon+0.3];
		else
			mapLatLim = [minLat-1, maxLat+1];
			mapLonLim = [minLon-1, maxLon+1];
		end
		
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1100,1200]);
		else
			figure('Position', [10,10,1100,1200]);
		end
		ax = usamap(mapLatLim,mapLonLim);
		states = shaperead('usastatehi','UseGeoCoords', true, 'BoundingBox', [mapLonLim', mapLatLim']);
		geoshow(radComp);
		geoshow(ax, states, 'FaceAlpha', 0,'EdgeColor','white');
		
		if doZoom
			geoshow(latSprlSmall,lonSprlSmall,'Color','black','LineWidth',20);
			geoshow(latSprlSmall,lonSprlSmall,'Color','white','LineWidth',10);
			geoshow(latSprlSmall(end),lonSprlSmall(end),'DisplayType','Point','Marker','o','MarkerFaceColor','black',...
				'MarkerEdgeColor','black','Markersize',50);
			geoshow(latSprlSmall(end),lonSprlSmall(end),'DisplayType','Point','Marker','o','MarkerFaceColor','white',...
				'MarkerEdgeColor','white','Markersize',40);
		else
			geoshow(latSprlSmall,lonSprlSmall,'Color','black','LineWidth',8);
			geoshow(latSprlSmall,lonSprlSmall,'Color','white','LineWidth',4);
% 			geoshow(latSprlSmall(end),lonSprlSmall(end),'DisplayType','Point','Marker','o','MarkerFaceColor','black',...
% 				'MarkerEdgeColor','black','Markersize',30);
% 			geoshow(latSprlSmall(end),lonSprlSmall(end),'DisplayType','Point','Marker','o','MarkerFaceColor','white',...
% 				'MarkerEdgeColor','white','Markersize',20);
		end
		
		if ~doZoom
			title([flight ' - ' times{ix} ' Radar Composite']);
		end

		set(findall(gcf,'-property','FontSize'),'FontSize',28);
		set(gcf,'Color','w');


		if saveFigs
			set(gcf,'InvertHardCopy','off');
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			if doLarge
				print([saveDir '/RadComps-SprlLocs/' flight '_ICCP_RadCompLg-' times{ix}],'-dpdf','-r0')
			elseif doZoom
				print([saveDir '/RadComps-SprlLocs/' flight '_ICCP_RadCompZm-' times{ix}],'-dpdf','-r0')
			else
				print([saveDir '/RadComps-SprlLocs/' flight '_ICCP_RadComp-' times{ix}],'-dpdf','-r0')
			end
		end
		
	end
	
end

if doStatic
	
	radComp = [dataPath 'IowaStateRadarComps/n0r_' flight staticImgT '.tif'];
	
	spiralIXsmall = find(time_secs >= startTiccp(1) & time_secs <= staticT);
	
	[~, sprlNum] = min(abs(endTiccp(2) - endT));
	spiralIX = find(time_secs >= startT(sprlNum) & time_secs <= endT(sprlNum)); % For map limits
	
	latSprlSmall = lat(spiralIXsmall);
	lonSprlSmall = lon(spiralIXsmall);
	
	% Determine map limits
	latSprl = lat(spiralIX);
	lonSprl = lon(spiralIX);
	maxLat = max(latSprl);
	minLat = min(latSprl);
	maxLon = max(lonSprl);
	minLon = min(lonSprl);
	if doLarge
		mapLatLim = [minLat-2, maxLat+2];
		mapLonLim = [minLon-2, maxLon+2];
	else
		mapLatLim = [minLat-1, maxLat+1];
		mapLonLim = [minLon-1, maxLon+1];
	end
	
	if saveFigs && noDisp
		figure('visible','off','Position', [10,10,1100,1200]);
	else
		figure('Position', [10,10,1100,1200]);
	end
	ax = usamap(mapLatLim,mapLonLim);
	states = shaperead('usastatehi','UseGeoCoords', true, 'BoundingBox', [mapLonLim', mapLatLim']);
	geoshow(radComp);
	geoshow(ax, states, 'FaceAlpha', 0,'EdgeColor','white');
	geoshow(latSprlSmall,lonSprlSmall,'Color','black','LineWidth',12);
	geoshow(latSprlSmall,lonSprlSmall,'Color','white','LineWidth',6);
	geoshow(latSprlSmall(end),lonSprlSmall(end),'DisplayType','Point','Marker','o','MarkerFaceColor','black',...
		'MarkerEdgeColor','black','Markersize',30);
	geoshow(latSprlSmall(end),lonSprlSmall(end),'DisplayType','Point','Marker','o','MarkerFaceColor','white',...
		'MarkerEdgeColor','white','Markersize',20);
	title([flight ' - ' staticImgT ' Radar Composite']);
	
	set(findall(gcf,'-property','FontSize'),'FontSize',28);
	set(gcf,'Color','w');
	
	
	if saveFigs
		set(gcf,'InvertHardCopy','off');
		set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		if doLarge
			print([saveDir '/RadComps-SprlLocs/' flight '_Static-ICCP_RadCompLg-' staticImgT],'-dpdf','-r0')
		else
			print([saveDir '/RadComps-SprlLocs/' flight '_Static-ICCP_RadComp-' staticImgT],'-dpdf','-r0')
		end
	end
		
	
end

if doAll
	for ix=1:length(startT)
	% for ix=1:1
		radComp = [dataPath 'IowaStateRadarComps/n0r_' flight times{ix} '.tif'];

		spiralIX = find(time_secs >= startT(ix) & time_secs <= endT(ix));

		latSprl = lat(spiralIX);
		lonSprl = lon(spiralIX);

		maxLat = max(latSprl);
		minLat = min(latSprl);
		maxLon = max(lonSprl);
		minLon = min(lonSprl);


		if doLarge
			mapLatLim = [minLat-2, maxLat+2];
			mapLonLim = [minLon-2, maxLon+2];
		else
			mapLatLim = [minLat-1, maxLat+1];
			mapLonLim = [minLon-1, maxLon+1];
		end

		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1100,1200]);
		else
			figure('Position', [10,10,1100,1200]);
		end
		ax = usamap(mapLatLim,mapLonLim);
		states = shaperead('usastatehi','UseGeoCoords', true, 'BoundingBox', [mapLonLim', mapLatLim']);
		geoshow(radComp);
		geoshow(ax, states, 'FaceAlpha', 0,'EdgeColor','white');
		geoshow(latSprl,lonSprl,'Color','black','LineWidth',6);
		geoshow(latSprl,lonSprl,'Color','white','LineWidth',3);
		geoshow(latSprl(end),lonSprl(end),'DisplayType','Point','Marker','o','MarkerFaceColor','black',...
			'MarkerEdgeColor','black','Markersize',20);
		geoshow(latSprl(end),lonSprl(end),'DisplayType','Point','Marker','o','MarkerFaceColor','white',...
			'MarkerEdgeColor','white','Markersize',15);
		geoshow(latSprl,lonSprl,'Color','white','LineWidth',4);
		geoshow(latSprl(end),lonSprl(end),'DisplayType','Point','Marker','o','MarkerFaceColor','black','Markersize',20);
		geoshow(latSprl(end),lonSprl(end),'DisplayType','Point','Marker','o','MarkerFaceColor','white','Markersize',15);
		title([flight ' - ' times{ix} ' Radar Composite - Spiral ' num2str(ix)]);

		set(findall(gcf,'-property','FontSize'),'FontSize',28);
		set(gcf,'Color','w');


		if saveFigs
			set(gcf,'InvertHardCopy','off');
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			if doLarge
				print([saveDir '/RadComps-SprlLocs/' flight '_RadCompLg-' times{ix} '-Spiral' num2str(ix)],'-dpdf','-r0')
			else
				print([saveDir '/RadComps-SprlLocs/' flight '_RadComp-' times{ix} '-Spiral' num2str(ix)],'-dpdf','-r0')
			end
		end
	end
end