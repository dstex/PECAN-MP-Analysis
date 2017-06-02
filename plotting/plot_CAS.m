close all;clearvars;

flight = '20150706';

probe = 'CIP'; % CIP/PIP irrelevant - Only used for pulling in FL data from sDist average files
avgTime = 10; % 5 or 10 will work here since we're using 1-sec data in this script

CAS_headLines = 37;
CAS_sampleArea = (0.25)*0.01; % mm2 -> cm2 %0.25 mm^2 suggested by Matt Freer


plotCASdist			= 0;
plotCASNt			= 0;
plotCASLWC			= 1;


saveFigs    = 1;
noDisp      = 0;
Ftype		= '-dpdf';
% Ftype		= '-dpng';

savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';
dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

startT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'startT');
endT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'endT');

sDistFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.' probe '.' num2str(avgTime) 'secAvg.mat'];

CASfile = [dataPath 'mp-data/' flight '/' 'CAS_' flight '.csv'];


%% Load and separate CAS probe data

% CAS_binMins = [0.61,0.68,0.75,0.82,0.89,0.96,1.03,1.1,1.17,1.25,1.5,...
% 	2,2.5,3,3.5,4,5,6.5,7.2,7.9,10.2,12.5,15,20,25,30,35,40,45,50];
CAS_binMins = [2,2.5,3,3.5,4,5,6.5,7.2,7.9,10.2,12.5,15,20,25,30,35,40,45,50];
binwidths = [diff(CAS_binMins),5]*1e-4; %Convert um to cm

CAS_all = csvread(CASfile,CAS_headLines);
timeSecs_CAS_all = CAS_all(:,1);
sumTransit_CAS_all = CAS_all(:,2);
sumParts_CAS_all = CAS_all(:,3);
% numConc_CAS_all = CAS_all(:,103:132); %All data
numConc_CAS_all = CAS_all(:,114:132); %Ignoring any bins less than 2 um
TAS_CAS_all = CAS_all(:,163);
Nt_CAS_RAW_all = CAS_all(:,164);
LWC_CAS_all = CAS_all(:,165);
MVD_CAS_all = CAS_all(:,166); % Median Volume Diameter
ED_CAS_all = CAS_all(:,167); % Effective Diameter

timeStep = [1;diff(timeSecs_CAS_all)];

for jj=1:length(CAS_binMins)
	sampleVol_CAS(:,jj) = ((TAS_CAS_all*100).*timeStep).*CAS_sampleArea; % (m/s -> cm/s)*sec*cm2 -> cm3
end

for ix=1:length(timeSecs_CAS_all)
	conc_CAS_all_1(ix,:) = numConc_CAS_all(ix,:)./binwidths; % cm-1
end

conc_CAS_all = conc_CAS_all_1./sampleVol_CAS; % cm-4

for ix=1:length(timeSecs_CAS_all)
	Nt_clcd_1(ix,:) = conc_CAS_all(ix,:).*binwidths; % cm-3
end

Nt_clcd = sum(Nt_clcd_1,2);

% display(['Sample Area (cm^2): ' num2str(CAS_sampleArea)])
% display(['TAS (ix=5078) (cm/s): ' num2str(TAS_CAS_all(5078)*100)])
% display(['Time step (ix=5078) (sec): ' num2str(timeStep(5078))])
% display(['Bin counts (ix=5078): ' num2str(numConc_CAS_all(5078,:))])
% display(['Sample Volume (ix=5078) (cm^3): ' num2str(sampleVol_CAS(5078,:))])
% display(['Counts normalized by binwidth (ix=5078) (cm^-1): ' num2str(conc_CAS_all_1(5078,:))])
% display(['N(D) (ix=5078) (cm^-4): ' num2str(conc_CAS_all(5078,:))])
% display(['Nt (ix=5078) (cm^-3): ' num2str(Nt_clcd(5078))])


load(sDistFile,'time_secsFL_orig','tempC_orig');

sprlNames = fieldnames(time_secs_orig);


%% Create directories to save plots in if they don't already exist
if saveFigs
    saveDir = [savePath flight];
	if (exist(saveDir, 'dir') ~= 7)
        mkdir(saveDir)
	end

	if (plotCASdist && exist([saveDir '/CAS-ND'], 'dir') ~= 7)
		mkdir([saveDir '/CAS-ND'])
	end
	if (plotCASNt && exist([saveDir '/CAS-Nt'], 'dir') ~= 7)
		mkdir([saveDir '/CAS-Nt'])
	end
	if (plotCASLWC && exist([saveDir '/CAS-LWC'], 'dir') ~= 7)
		mkdir([saveDir '/CAS-LWC'])
	end
end

%% Plotting

if plotCASdist
	for ix = 1:length(startT)
		if saveFigs && noDisp
            figure('visible','off','Position', [10,10,1500,1000]);
		else
            figure('Position', [10,10,1500,1000]);
		end
		
		sprlIx_CAS = find(timeSecs_CAS_all >= startT(ix) & timeSecs_CAS_all < endT(ix));
		timeSecs_CAS = timeSecs_CAS_all(sprlIx_CAS);
		conc_CAS = conc_CAS_all(sprlIx_CAS,:);
		
		stairs(CAS_binMins, mean(conc_CAS,1), 'b', 'LineWidth', 2);
		title([flight ' - Spiral ' num2str(ix) ' - CAS N(D)']);
		
		xlabel('D [\mum]');
		ylabel('N(D) [cm^{-4}]');
		set(gca,'Yscale','log');
 		ylim([1e0 1.5e5]);
		set(gca,'XMinorTick','on','YMinorTick','on');
		set(findall(gcf,'-property','FontSize'),'FontSize',28)
		grid
		
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/CAS-ND/' flight '_CAS_ND_S' num2str(ix)],Ftype,'-r0')
		end
	end
	
end

if plotCASNt
	for ix = 1:length(startT)
		if saveFigs && noDisp
            figure('visible','off','Position', [10,10,1500,1000]);
		else
            figure('Position', [10,10,1500,1000]);
		end
		
		sprlIx_CAS = find(timeSecs_CAS_all >= startT(ix) & timeSecs_CAS_all < endT(ix));
		timeSecs_CAS = timeSecs_CAS_all(sprlIx_CAS);
		Nt_CAS_RAW_all_sprl = Nt_CAS_RAW_all(sprlIx_CAS);
		Nt_CAS = Nt_clcd(sprlIx_CAS);
		
		plot(timeSecs_CAS/24/3600, Nt_CAS,'b-');
		title([flight ' - Spiral ' num2str(ix) ' - CAS N_t']);
		hold on
		plot(timeSecs_CAS/24/3600, Nt_CAS_RAW_all_sprl,'r-');
		legend('UIUC','RAW');
		
		datetick('x','HH:MM:SS');
		ylabel('N_t [cm^{-3}]')
% 		ylim([1e0 1e5]);
		set(gca,'XMinorTick','on','YMinorTick','on','XTickLabelRotation',45);
		set(findall(gcf,'-property','FontSize'),'FontSize',28)
		grid
		
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/CAS-Nt/' flight '_CAS_Nt_S' num2str(ix)],Ftype,'-r0')
		end
	end
	
end

if plotCASLWC
	for ix = 1:length(startT)
		if saveFigs && noDisp
            figure('visible','off','Position', [10,10,1500,1000]);
		else
            figure('Position', [10,10,1500,1000]);
		end
		
		sprlIx_CAS = find(timeSecs_CAS_all >= startT(ix) & timeSecs_CAS_all < endT(ix));
		timeSecs_CAS = timeSecs_CAS_all(sprlIx_CAS);
		LWC_CAS = LWC_CAS_all(sprlIx_CAS);
		
		plot(timeSecs_CAS/24/3600, LWC_CAS,'b-');
		title([flight ' - Spiral ' num2str(ix) ' - CAS LWC']);
		
		datetick('x','HH:MM:SS');
		ylabel('LWC [g m^{-3}]')
% 		ylim([1e-3 1e3]);
		set(gca,'XMinorTick','on','YMinorTick','on','XTickLabelRotation',45);
		set(findall(gcf,'-property','FontSize'),'FontSize',28)
		grid
		
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/CAS-LWC/' flight '_CAS_LWC_S' num2str(ix)],Ftype,'-r0')
		end
	end
	
end