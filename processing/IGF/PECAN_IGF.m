clearvars; close all;

%% Specify various plotting/calculation parameters
flight = '20150706';

allSprls	= 1;

extendCIP	= 1;

plotAllPIP	= 1; % Consider all PIP data regardless of good/bad designation
plotGoodPIP	= 1; % Plot only the PIP periods considered good
plotBadPIP	= 0; % Plot only the PIP periods considered bad

stairFit	= 0; % Plot fit using stairs? Line plot otherwise

saveMat		= 0; % Save mat files for each spiral

saveFigs	= 1;
noDisp		= 1;
Ftype = '-dpdf';

savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';
dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

% Create save directories if they don't exist
if saveFigs
    saveDir = [savePath flight];
	if (exist(saveDir, 'dir') ~= 7)
        mkdir(saveDir)
	end
	if (exist([saveDir '/IGF'], 'dir') ~= 7)
		mkdir([saveDir '/IGF'])
	end
	if (saveMat && exist([saveDir '/IGF/matFiles'], 'dir') ~= 7)
		mkdir([saveDir '/IGF/matFiles'])
	end
end

%% Load in various PECAN parameters
PIP_acptStartT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'PIP_acptStartT');
PIP_acptEndT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'PIP_acptEndT');
PIP_rjctStartT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'PIP_rjctStartT');
PIP_rjctEndT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'PIP_rjctEndT');

%% Load in CIP and PIP SD files and then extract only the variables we need from them
cipSDall = load([dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.CIP.10secAvg.mat']);
pipSDall = load([dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.PIP.10secAvg.mat']);

cip_concMinR = cipSDall.conc_minR_orig;
cip_binMin = (cipSDall.bin_min)*1000; % convert mm to um
cip_binMax = (cipSDall.bin_max)*1000;
cip_binMid = (cipSDall.bin_mid)*1000;
cip_binEdges = [cip_binMin; cip_binMax(end)]; 
cip_timSecs = cipSDall.time_secs_orig;

pip_concMinR = pipSDall.conc_minR_orig;
pip_binMin = (pipSDall.bin_min)*1000;
pip_binMax = (pipSDall.bin_max)*1000;
pip_binMid = (pipSDall.bin_mid)*1000;
pip_binEdges = [pip_binMin; pip_binMax(end)];
pip_timSecs = pipSDall.time_secs_orig;

sprlNames = fieldnames(cip_concMinR); % Variable used unimportant - just needs to be one of the structs

%% Create loop vector to hold spiral numbers we intend to plot
if ~allSprls
	switch flight
		case '20150706'
			loopVctr = [1,3,5,7]; % Spiral numbers to plot
	end
else
	loopVctr = 1:length(sprlNames);
end

%% Create array of extended CIP bins
% Get increment between bin mids for CIP and then create an extended bin mids array
% spanning from the beginning of the CIP bins to the end of the PIP bins
if extendCIP
	% Get increment between bin mids for CIP and then create an extended bin mids array
	% spanning from the beginning of the CIP bins to the end of the PIP bins
	cip_binDiff = diff(cip_binMid);
	cipExt_incrmnt = ceil(cip_binDiff(end));
	cipExt_binMidHalf = cip_binMid(end)+cipExt_incrmnt:cipExt_incrmnt:pip_binMid(end);
	cipExt_binMid = [cip_binMid; cipExt_binMidHalf'];
end


%% Plotting

% Consider all PIP data when plotting the CIP/PIP PSDs and fits
if plotAllPIP
	for ix = loopVctr
		cipPSD = nanmean(cip_concMinR.(sprlNames{ix}),1);
		pipPSD = nanmean(pip_concMinR.(sprlNames{ix}),1);

		[cip_nml, cip_mObs, cip_mFit, cip_chisquare, cipPSD_fit, cip_fitFlag, cip_fitDetails] = f_one_mode(cipPSD, cip_binEdges', [0 2 3]);
		[pip_nml, pip_mObs, pip_mFit, pip_chisquare, pipPSD_fit, pip_fitFlag, pip_fitDetails] = f_one_mode(pipPSD, pip_binEdges', [0 2 3]);

		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end

		stairs(cip_binMid/1000,cipPSD,'r','LineWidth', 1.5);
		hold on
		if extendCIP
			cipExtPSD_fit = 10^cip_nml(1).*cipExt_binMid.^cip_nml(2).*exp(-cip_nml(3).*cipExt_binMid);
			if stairFit
				stairs(cipExt_binMid/1000, cipExtPSD_fit,'b--','LineWidth', 1.5);
				stairs(cip_binMid/1000, cipPSD_fit,'b','LineWidth', 1.5);
			else
				plot(cipExt_binMid/1000,cipExtPSD_fit,'b--','LineWidth', 1.5);
				plot(cip_binMid/1000,cipPSD_fit,'b','LineWidth', 1.5);
			end
		else
			if stairFit
				stairs(cip_binMid/1000, cipPSD_fit,'b','LineWidth', 1.5);
			else
				plot(cip_binMid/1000,cipPSD_fit,'b','LineWidth', 1.5);
			end
		end

		stairs(pip_binMid/1000,pipPSD,'m','LineWidth', 1.5);
		if stairFit
			stairs(pip_binMid/1000, pipPSD_fit,'k','LineWidth', 1.5);
		else
			plot(pip_binMid/1000,pipPSD_fit,'k','LineWidth', 1.5);
		end

		hold off

		title([flight ' - Spiral ' num2str(ix) ' - CIP and PIP N(D) Fits'],'FontSize',24);

		if extendCIP
			lgnd = legend('CIP Obs','CIP Extnd IGF','CIP IGF','PIP Obs','PIP IGF');
		else
			lgnd = legend('CIP Obs','CIP IGF','PIP Obs','PIP IGF');
		end
		lgnd.FontSize = 18;
		xlabel('D [mm]');
		ylabel('N(D) [cm^{-4}]');
		set(gca,'yscale','log','xscale','log');
		grid minor

		ax = ancestor(gca, 'axes');
		xRule = ax.XAxis;
		yRule = ax.YAxis;
		xRule.FontSize = 24;
		yRule.FontSize = 24;

		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/IGF/' flight '_IGF-CIP-allPIP_S' num2str(ix)],Ftype,'-r0')
		end

		if saveMat
			save([saveDir '/IGF/matFiles/' flight '_IGF-CIP-allPIP_S' num2str(ix)],'-regexp',...
				'^(?!(cipSDall|pipSDall|cip_concMinR|pip_concMinR|Ftype|ax|allSprls|dataPath|lgnd|noDisp|pos|saveDir|saveFigs|saveMat|savePath|sprlNames|stairFit|xRule|yRule)$).');
		end

	end
end

% Consider only good PIP data when plotting the CIP/PIP PSDs and fits
if plotGoodPIP
	for ix = loopVctr
		pipTime = pip_timSecs.(sprlNames{ix});
		pipGdIx = find(pipTime >= PIP_acptStartT(ix) & pipTime <= PIP_acptEndT(ix));
		
		pipConc = pip_concMinR.(sprlNames{ix});
		
		cipPSD = nanmean(cip_concMinR.(sprlNames{ix}),1);
		pipPSD = nanmean(pipConc(pipGdIx,:),1);

		[cip_nml, cip_mObs, cip_mFit, cip_chisquare, cipPSD_fit, cip_fitFlag, cip_fitDetails] = f_one_mode(cipPSD, cip_binEdges', [0 2 3]);
		[pip_nml, pip_mObs, pip_mFit, pip_chisquare, pipPSD_fit, pip_fitFlag, pip_fitDetails] = f_one_mode(pipPSD, pip_binEdges', [0 2 3]);

		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end

		stairs(cip_binMid/1000,cipPSD,'r','LineWidth', 1.5);
		hold on
		if extendCIP
			cipExtPSD_fit = 10^cip_nml(1).*cipExt_binMid.^cip_nml(2).*exp(-cip_nml(3).*cipExt_binMid);
			if stairFit
				stairs(cipExt_binMid/1000, cipExtPSD_fit,'b--','LineWidth', 1.5);
				stairs(cip_binMid/1000, cipPSD_fit,'b','LineWidth', 1.5);
			else
				plot(cipExt_binMid/1000,cipExtPSD_fit,'b--','LineWidth', 1.5);
				plot(cip_binMid/1000,cipPSD_fit,'b','LineWidth', 1.5);
			end
		else
			if stairFit
				stairs(cip_binMid/1000, cipPSD_fit,'b','LineWidth', 1.5);
			else
				plot(cip_binMid/1000,cipPSD_fit,'b','LineWidth', 1.5);
			end
		end

		stairs(pip_binMid/1000,pipPSD,'m','LineWidth', 1.5);
		if stairFit
			stairs(pip_binMid/1000, pipPSD_fit,'k','LineWidth', 1.5);
		else
			plot(pip_binMid/1000,pipPSD_fit,'k','LineWidth', 1.5);
		end

		hold off

		title([flight ' - Spiral ' num2str(ix) ' - CIP and PIP (Good) N(D) Fits'],'FontSize',24);

		if extendCIP
			lgnd = legend('CIP Obs','CIP Extnd IGF','CIP IGF','PIP Obs','PIP IGF');
		else
			lgnd = legend('CIP Obs','CIP IGF','PIP Obs','PIP IGF');
		end
		lgnd.FontSize = 18;
		xlabel('D [mm]');
		ylabel('N(D) [cm^{-4}]');
		set(gca,'yscale','log','xscale','log');
		grid minor

		ax = ancestor(gca, 'axes');
		xRule = ax.XAxis;
		yRule = ax.YAxis;
		xRule.FontSize = 24;
		yRule.FontSize = 24;

		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/IGF/' flight '_IGF-CIP-goodPIP_S' num2str(ix)],Ftype,'-r0')
		end

		if saveMat
			save([saveDir '/IGF/matFiles/' flight '_IGF-CIP-goodPIP_S' num2str(ix)],'-regexp',...
				'^(?!(cipSDall|pipSDall|cip_concMinR|pip_concMinR|Ftype|ax|allSprls|dataPath|lgnd|noDisp|pos|saveDir|saveFigs|saveMat|savePath|sprlNames|stairFit|xRule|yRule)$).');
		end

	end
end