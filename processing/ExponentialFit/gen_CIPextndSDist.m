clearvars; close all;

%% Specify various plotting/calculation parameters
flight = '20150706';

trimCipBeg	= 1;
trimCipEnd	= 0;
trimCipBoth	= 0;
trimIncld = 5:28; % Encompasses diams. of 125-1000 um
 

plotNDintvl = 1; % Plot N(D) for every time step
plotND		= 0;

allSprls	= 1;

use10secAvg = 1; % If false, will use original 1-sec SDs

% Brown and Francis 1995 coefficients
a = 7.38e-11;
b = 1.9;

saveMat = 1;

saveFigs	= 1;
noDisp		= 1;
Ftype = '-dpdf';

savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';
dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';




%% Load in various PECAN parameters
mlBotTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlBotTime');
mlTopTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlTopTime');
PIP_acptStartT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'PIP_acptStartT');
PIP_acptEndT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'PIP_acptEndT');
PIP_rjctStartT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'PIP_rjctStartT');
PIP_rjctEndT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'PIP_rjctEndT');



%% Load in CIP and PIP SD files and then extract only the variables we need from them
cipSDall = load([dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.CIP.10secAvg.mat']);
pipSDall = load([dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.PIP.10secAvg.mat']);

if use10secAvg
	cip_concMinR = cipSDall.conc_minR_avg; %cm-4
	pip_concMinR = pipSDall.conc_minR_avg;
	cip_timeSecs = cipSDall.time_secs_avg;
	pip_timeSecs = pipSDall.time_secs_avg;
	cip_smplVol = cipSDall.sampleVol_avg; %cm3
	pip_smplVol = pipSDall.sampleVol_avg;
	cip_massTWC = cipSDall.mass_twc_avg; %g/cm4
	pip_massTWC = pipSDall.mass_twc_avg;
	cip_TWC = cipSDall.twc_avg; %g/m3
	pip_TWC = pipSDall.twc_avg;
else
	cip_concMinR = cipSDall.conc_minR_orig;
	pip_concMinR = pipSDall.conc_minR_orig;
	cip_timeSecs = cipSDall.time_secs_orig;
	pip_timeSecs = pipSDall.time_secs_orig;
	cip_smplVol = cipSDall.sampleVol_orig;
	pip_smplVol = pipSDall.sampleVol_orig;
	cip_massTWC = cipSDall.mass_twc_orig;
	pip_massTWC = pipSDall.mass_twc_orig;
	cip_TWC = cipSDall.twc_orig;
	pip_TWC = pipSDall.twc_orig;
end

% Set bins to be included in CIP variables
if trimCipBoth
	cipIncld = trimIncld;
	titlStr = ' - Obs125-1000\mum';
	fileStr = '_Obs125-1000um';
elseif trimCipBeg
	cipIncld = trimIncld(1):size(cipSDall.conc_minR_orig.sprl1,2);
	titlStr = ' - Obs125-2000\mum';
	fileStr = '_Obs125-2000um';
elseif trimCipEnd
	cipIncld = 1:trimIncld(end);
	titlStr = ' - Obs25-1000\mum';
	fileStr = '_Obs25-1000um';
else
	cipIncld = 1:size(cipSDall.conc_minR_orig.sprl1,2);
	titlStr = '';
	fileStr = '';
end

% Create save directories if they don't exist
if saveFigs
    saveDir = [savePath flight];
	if (exist(saveDir, 'dir') ~= 7)
        mkdir(saveDir)
	end
	if (exist([saveDir '/expFit/10sFits' fileStr], 'dir') ~= 7)
		mkdir([saveDir '/expFit/10sFits' fileStr])
	end
end

cip_binMin = (cipSDall.bin_min(cipIncld))/10; % cm
cip_binMinAll = (cipSDall.bin_min)/10;
cip_binMax = (cipSDall.bin_max(cipIncld))/10;
cip_binMaxAll = (cipSDall.bin_max)/10;
cip_binMid = (cipSDall.bin_mid(cipIncld))/10;
cip_binMidAll = (cipSDall.bin_mid)/10; % Used when constructing extended cip bins
cip_binEdges = [cip_binMin; cip_binMax(end)];% cm
cip_binEdgesAll = [cip_binMinAll; cip_binMaxAll(end)]; % Used when constructing extended cip bins
cip_binwidth = (cipSDall.bin_size(cipIncld))/10; % cm

pip_binMin = (pipSDall.bin_min)/10;
pip_binMax = (pipSDall.bin_max)/10;
pip_binMid = (pipSDall.bin_mid)/10;
pip_binEdges = [pip_binMin; pip_binMax(end)];
pip_binwidth = (pipSDall.bin_size)/10;

sprlNames = fieldnames(cip_concMinR); % Variable used unimportant - just needs to be one of the structs




%% Find indices corresponding to ML top/bottom for each CIP spiral
mlBotIx = NaN(size(mlBotTime)); % Indices will be relative to the spiral, NOT the whole flight
mlTopIx = NaN(size(mlTopTime));
for iz = 1:length(sprlNames)
	if isnan(mlBotTime(iz))
		mlBotIx(iz) = NaN;
	else
		[~, mlBotIx(iz)] = min(abs(cip_timeSecs.(sprlNames{iz}) - mlBotTime(iz))); % Find the closest match
	end
	if isnan(mlTopTime(iz))
		mlTopIx(iz) = NaN;
	else
		[~, mlTopIx(iz)] = min(abs(cip_timeSecs.(sprlNames{iz}) - mlTopTime(iz)));
	end
end

%% Create loop vector to hold spiral numbers we intend to plot
if ~allSprls
	switch flight
		case '20150706'
			loopVctr = [2:8]; % Spiral numbers to plot
	end
else
	loopVctr = 1:length(sprlNames);
end

%% Create array of extended CIP bins
% Get increment between bin mids for CIP and then create an extended bin mids array
% spanning from the beginning of the CIP bins to the end of the PIP bins

% cip_binDiff = diff(cip_binMid); % cm
% cipExt_incrmnt = cip_binDiff(end); % cm
% cipExt_binMidHalf = cip_binMid(end)+cipExt_incrmnt:cipExt_incrmnt:pip_binMid(end);
% cipExt_binMid = [cip_binMid; cipExt_binMidHalf'];

cipExt_binMid = [cip_binMidAll; pip_binMid(12:end)]; % cm
cipExt_binEdges = [cip_binEdgesAll; pip_binEdges(13:end)];
cipExt_binMin = cipExt_binEdges(1:end-1);
cipExt_binMax = cipExt_binEdges(2:end);
cipExt_binwidth = diff(cipExt_binEdges); % cm
numExt_bins = length(cipExt_binMid);

%% Loop through each spiral we want to calculate
for ix = loopVctr
	fitSkip = 0;
	lmdaCount = 0;
	
	if mlTopIx(ix) > mlBotIx(ix)
		sprlUp = 1;
	else
		sprlUp = 0;
	end
	
	fprintf('Currently working on Spiral %d\n',ix);

	% Create new extended CIP structs, with NaN arrays for each spiral
	cipConc = cip_concMinR.(sprlNames{ix})(:,cipIncld);
	cipMassTWC = cip_massTWC.(sprlNames{ix})(:,cipIncld);
	cipSmplVol = cip_smplVol.(sprlNames{ix})(:,cipIncld);
	
	cipConc_ext.(sprlNames{ix}) = NaN(size(cipConc,1),numExt_bins); % N(D) using only the exponential fit as a basis
	cipConc_hybrid.(sprlNames{ix}) = NaN(size(cipConc,1),numExt_bins); % Will contain observed CIP data, with fit values beyond D=1.6mm
	
	%{
	cipMassTWC_ext.(sprlNames{ix}) = NaN(size(cipMassTWC,1),numExt_bins);
	cipMassTWC_hybrid.(sprlNames{ix}) = NaN(size(cipMassTWC,1),numExt_bins);
	cipTWC_ext.(sprlNames{ix}) = NaN(size(cip_TWC));
	cipTWC_hybrid.(sprlNames{ix}) = NaN(size(cip_TWC));
	cipSmplVol_ext.(sprlNames{ix}) = NaN(size(cipSmplVol,1),numExt_bins);
	cipSmplVol_hybrid.(sprlNames{ix}) = NaN(size(cipSmplVol,1),numExt_bins);
	%}
	cip_coeff_all.(sprlNames{ix}) = NaN(size(cipConc,1),2);
	cipPSDobsFit.(sprlNames{ix}) = NaN(size(cipConc));
	
	for ii = 1:size(cipConc,1)
% 	for ii = 10
		
		
		% Pull out the distribution at each time step ii
		cipPSD = cipConc(ii,:);
		cipMTWC = cipMassTWC(ii,:);
		
		% Check to see if the entire distribution is NaN to reduce issues with the fit
		% Also, only attempt the fit if we have 10 or more values to fit to
		% Actually perform the fit and get N0/lambda values for each time step
		if (all(isnan(cipConc(ii,:))) || length(find(cipConc(ii,:) > 0)) < 3)
			fprintf('Sizing %d/%d\tInsufficient data\n',ii,size(cipConc,1))
			cip_coeff_all.(sprlNames{ix})(ii,:) = NaN;
			cipPSDobsFit.(sprlNames{ix})(ii,:) = NaN;
			gof.(sprlNames{ix})(ii) = struct('sse',{NaN},'rsquare',{NaN},'dfe',{NaN},'adjrsquare',{NaN},'rmse',{NaN});
			
			fitSkip = fitSkip+1;
			fileTail = '_lt3';
		else		
			[cip_coeff_all.(sprlNames{ix})(ii,:), cipPSDobsFit.(sprlNames{ix})(ii,:), fitresult, gof.(sprlNames{ix})(ii)] = robustLMExpFit(cip_binMid, cipPSD');
				
			if (cip_coeff_all.(sprlNames{ix})(ii,2) < 0)
				fprintf('Sizing %d/%d\tLambda < 0\n',ii,size(cipConc,1))
				lmdaCount = lmdaCount + 1;
				fileTail = '_X';
			else
				fprintf('Sizing %d/%d\n',ii,size(cipConc,1))
				fileTail = '';
			end
		end
		
		% Construct extended CIP concentration distributions for each time step
		% Hybrid: first length(cip_binMid) bins contain observed CIP PSD; bins length(cip_binMid)+1:end use N0/lambda values to compute extended distribution for each assocaited time
		% Ext: Use the N0/lambda values to distribution over all bins
		cipConc_hybrid.(sprlNames{ix})(ii,1:length(cip_binMid)) = cipPSD;
		n0 = cip_coeff_all.(sprlNames{ix})(ii,1);
		lmda = cip_coeff_all.(sprlNames{ix})(ii,2);

		cipConc_hybrid.(sprlNames{ix})(ii,length(cip_binMid)+1:end) = n0.*exp(-lmda.*cipExt_binMid(length(cip_binMid)+1:end));
		cipConc_ext.(sprlNames{ix})(ii,:) = n0.*exp(-lmda.*cipExt_binMid);

		%{
		% Construct extended CIP sample volume for each time step
		% Hybrid: first length(cip_binMid) bins contain observed CIP sample volumes; bins length(cip_binMid)+1:end use average of observed sample volumes for each assocaited time
		% Ext: Use average of observed sample volumes for all bins
		cipSmplVol_hybrid.(sprlNames{ix})(ii,1:length(cip_binMid)) = cipSmplVol(ii,:);
		cipSmplVol_hybrid.(sprlNames{ix})(ii,length(cip_binMid)+1:end) = nanmean(cipSmplVol(ii,:));
		cipSmplVol_ext.(sprlNames{ix})(ii,:) = ones(size(cipSmplVol_ext.(sprlNames{ix})(ii,:))).*nanmean(cipSmplVol(ii,:)); % Assume the average sample volume for entire bin range
		
		% Construct extended CIP mass distribution for each time step
		%	Final result should be g/cm4, and testing suggests that multiplying by sample volume and binwidth is not needed here
		% Hybrid: first length(cip_binMid) bins contain observed CIP mass(TWC); bins length(cip_binMid)+1:end use either mass of liquid sphere below melting layer,
		%	or Brown&Francis95 a-b coefficients above the melting layer **(RELATIONSHIP REQUIRES Diam have units of um)**
		cipMassTWC_hybrid.(sprlNames{ix})(ii,1:length(cip_binMid)) = cipMTWC;
		for iii = length(cip_binMid)+1:length(cipExt_binMid)
			if ((sprlUp && ii < mlBotIx(ix)) || (~sprlUp && ii > mlBotIx(ix))) % below melting layer
				cipMassTWC_hybrid.(sprlNames{ix})(ii,iii) = ((pi/6)*(cipExt_binMid(iii)/10)^3) * (cipConc_ext.(sprlNames{ix})(ii,iii));% * cipSmplVol_hybrid.(sprlNames{ix})(ii,iii) * (cipExt_binwidth(iii)/10));
			else
				cipMassTWC_hybrid.(sprlNames{ix})(ii,iii) = (a*(cipExt_binMid(iii)*10000)^b) * (cipConc_ext.(sprlNames{ix})(ii,iii));% * cipSmplVol_hybrid.(sprlNames{ix})(ii,iii) * (cipExt_binwidth(iii)/10));
			end
		end
		% Ext: Use above melting layer dependencies to calculate mass over whole distribution
		for iv = 1:length(cipExt_binMid)
			if ((sprlUp && ii < mlBotIx(ix)) || (~sprlUp && ii > mlBotIx(ix))) % below melting layer
				cipMassTWC_ext.(sprlNames{ix})(ii,iv) = ((pi/6)*(cipExt_binMid(iii)/10)^3) * (cipConc_ext.(sprlNames{ix})(ii,iv));% * cipSmplVol_hybrid.(sprlNames{ix})(ii,iv) * (cipExt_binwidth(iii)/10));
			else
				cipMassTWC_ext.(sprlNames{ix})(ii,iv) = (a*(cipExt_binMid(iii)*10000)^b) * (cipConc_ext.(sprlNames{ix})(ii,iv));%  * cipSmplVol_hybrid.(sprlNames{ix})(ii,iv) * (cipExt_binwidth(iii)/10));
			end
		end
		
		% Calculate TWC
		% Convert hybrid and extended dists to g/cm3 from g/cm4
		% Then, sum total mass in each bin, and convert to g/m3
		cipMassHybd_t = cipMassTWC_hybrid.(sprlNames{ix})(ii,:).*(cipExt_binwidth'); % Convert to g/cm3
		cipMassExt_t = cipMassTWC_ext.(sprlNames{ix})(ii,:).*(cipExt_binwidth');
		cipTWC_hybrid.(sprlNames{ix})(ii) = nansum(cipMassHybd_t,2)*1e6; % Sum masses in each bin and convert to g/m3
		cipTWC_ext.(sprlNames{ix})(ii) = nansum(cipMassExt_t,2)*1e6;
		%}
		
		
		if plotNDintvl
			cipTime = insec2hhmmss(cip_timeSecs.(sprlNames{ix})(ii));
			if saveFigs && noDisp
				figure('visible','off','Position', [10,10,1200,700]);
			else
				figure('Position', [10,10,1200,700]);
			end
			
			stairs(cip_binMin*10,cipPSD,'r','LineWidth', 1.5);
			hold on
			plot(cipExt_binMin*10,cipConc_ext.(sprlNames{ix})(ii,:),'b','LineWidth', 1.5);
% 			plot(cip_binMin*10,cipPSDobsFit.(sprlNames{ix})(ii,:),'k','LineWidth', 1.5);
			title([flight ' - CIP Exp. Fit - Spiral ' num2str(ix) titlStr ' - ' num2str(cipTime) ' - ' num2str(ii)],'FontSize',24);
			
			lgnd = legend('CIP Obs','CIP Extnd expFit'); %,'CIP expFit');
			lgnd.FontSize = 18;
			xlabel('D [mm]');
			ylabel('N(D) [cm^{-4}]');
			set(gca,'yscale','log');
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
				print([saveDir '/expFit/10sFits' fileStr '/' flight '_expFit-CIP_S' num2str(ix) '-' sprintf('%2.2d',ii) '-' num2str(cipTime) fileStr fileTail],Ftype,'-r0')
			end
		end
	end
	
	fprintf('Fitting skipped due to insufficient data %d times\n',fitSkip);
	fprintf('Lambda parameter was negative %d/%d times\n\n',lmdaCount,size(cipConc,1));
	
	% Plot N(D) averaged over whole spiral(s)
	if plotND
		cipConcExt = cipConc_ext_new.(sprlNames{ix});
		cipObsFit = cipPSDobsFit.(sprlNames{ix});
		
		if saveFigs && noDisp
			figure('visible','off','Position', [10,10,1200,700]);
		else
			figure('Position', [10,10,1200,700]);
		end
		
		stairs(cip_binMin*10,nanmean(cipConc,1),'r','LineWidth', 1.5);
		hold on
		plot(cipExt_binMin*10,nanmean(cipConcExt,1),'b','LineWidth', 1.5);
		plot(cip_binMin*10,nanmean(cipObsFit,1),'k','LineWidth', 1.5);
		title([flight ' - CIP Exp. Fit - Spiral ' num2str(ix) titlStr],'FontSize',24);
		
		lgnd = legend('CIP Obs','CIP Extnd expFit','CIP expFit');
		lgnd.FontSize = 18;
		xlabel('D [mm]');
		ylabel('N(D) [cm^{-4}]');
		set(gca,'yscale','log');
		grid;
		
		ax = ancestor(gca, 'axes');
		xRule = ax.XAxis;
		yRule = ax.YAxis;
		xRule.FontSize = 24;
		yRule.FontSize = 24;
		
		if saveFigs
			set(gcf,'Units','Inches');
			pos = get(gcf,'Position');
			set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
			print([saveDir '/expFit/' flight '_expFit-CIP_S' num2str(ix) fileStr],Ftype,'-r0')
		end
	end
end

if saveMat
	save([saveDir '/expFit/' flight '_expFit-CIP' fileStr])
end

