clearvars; close all;

%% Specify various plotting/calculation parameters
flight = '20150706';

avgTime = 1;

numValid = 29;

cipIncld = 8:34;
idStr = ' - CIP Obs > 200 \mum';
fileIdStr = ['_ObsDgt200um_gte' num2str(numValid) 'bins'];


doFitting	= 1;

doWhole		= 0;
doTempBins	= 1;
doEvryTStp	= 0;

doIGF		= 1;
doExpLMFit	= 0;

saveMat		= 1;
mFileId = '_Fit-CIP_TempBins';



doPlotting	= 1;

plotWhole		= 0;
plotTempBins	= 1;
plotEvryTStp	= 0;
plotAllGood		= 0;
plotAllSkip		= 0;
plotAllNegLmda	= 0;
saveFigs	= 1;
noDisp		= 1;
Ftype		= '-dpng';


savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';
dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';


% Brown and Francis 1995 coefficients
a = 7.38e-11;
b = 1.9;


% Create save directories if they don't exist
saveDir = [savePath flight];
if saveFigs
	if (exist(saveDir, 'dir') ~= 7)
        mkdir(saveDir)
	end
	if (exist([saveDir '/CIP-Fitting'], 'dir') ~= 7)
		mkdir([saveDir '/CIP-Fitting'])
	end
end


%% Load in CIP SD file and then extract only the variables we need from it
if avgTime == 1 || avgTime == 10
	cipDataF = load([dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.CIP.10secAvg.mat']);
elseif avgTime == 5
	cipDataF = load([dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.CIP.5secAvg.mat']);
end

if avgTime == 1
	cip_concMinR = cipDataF.conc_minR_orig;
	cip_massTWC = cipDataF.mass_twc_orig;
	cip_timeSecs = cipDataF.time_secs_orig;
	tempC = cipDataF.tempC_orig;
else
	cip_concMinR = cipDataF.conc_minR_avg;
	cip_massTWC = cipDataF.mass_twc_avg;
	cip_timeSecs = cipDataF.time_secs_avg;
	tempC = cipDataF.tempC_avg;
end
cip_binMin = (cipDataF.bin_min)./10; % convert mm to cm
cip_binMax = (cipDataF.bin_max)./10;
cip_binMid = (cipDataF.bin_mid)./10;
cip_binEdges = [cip_binMin; cip_binMax(end)]; 

pip_binMidEnd = 7300.0; %um

%% Load in various PECAN parameters
mlBotTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlBotTime');
mlTopTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlTopTime');
mlBotTemp = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlBotTemp');
mlTopTemp = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlTopTemp');

%% Create array of extended CIP bins
% Get increment between bin mids for CIP and then create an extended bin mids array
% spanning from the beginning of the CIP bins to the end of the PIP bins
pip_prtlBinEdges = [2200.0	2400.0	2600.0	2800.0	3000.0	3200.0	3400.0	3600.0	3800.0	4000.0	4200.0	4400.0...
						4600.0	4800.0	5000.0	5200.0	5500.0	5800.0	6100.0	6400.0	6700.0	7000.0	7300.0]/10000; %cm

cipExt_binEdges = [cip_binEdges; pip_prtlBinEdges'];
cipExt_binMin = cipExt_binEdges(1:end-1);
cipExt_binMax = cipExt_binEdges(2:end);
cipExt_binMid = (cipExt_binMin+cipExt_binMax)/2;
cipExt_binwidth = diff(cipExt_binEdges); 
numExt_bins = length(cipExt_binMid);

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


% loopVctr = 1:length(sprlNames);
loopVctr = 2;

%% Run the fitting analyses
if doFitting
	for ix = loopVctr
		fprintf('Currently working on Spiral %d\n',ix);
		
		if mlTopIx(ix) > mlBotIx(ix)
			sprlUp = 1;
		else
			sprlUp = 0;
		end
		
		tempC_sprlOrig = tempC.(sprlNames{ix});
		cipConc = cip_concMinR.(sprlNames{ix});
		cipConcSprl = nanmean(cip_concMinR.(sprlNames{ix}),1);
		cipMass = cip_massTWC.(sprlNames{ix});
		cipMassSprl = nanmean(cip_massTWC.(sprlNames{ix}),1);
		
		% Calculate fits over entire spiral
		if doWhole
			if doIGF
				cipConc_ext_igfWhl.(sprlNames{ix}) = NaN(1,numExt_bins);
				cipConc_hybrid_igfWhl.(sprlNames{ix}) = NaN(1,numExt_bins);
				cipMass_ext_igfWhl.(sprlNames{ix}) = NaN(1,numExt_bins);
				cipMass_hybrid_igfWhl.(sprlNames{ix}) = NaN(1,numExt_bins);
				
				[cip_igf_nmlWhl.(sprlNames{ix}), cip_igf_mmtsObsWhl.(sprlNames{ix}), cip_igf_mmtsFitWhl.(sprlNames{ix}),...
					cip_igf_chisquareWhl.(sprlNames{ix}), cipObs_igf_fitWhl.(sprlNames{ix}), cip_igf_fitFlagWhl.(sprlNames{ix}),...
					cip_igf_fitDetailsWhl] = igfFit(cipConcSprl(cipIncld), cip_binEdges(cipIncld(1):cipIncld(end)+1)', [0 2 3]);
				
				n0_tmp = cip_igf_nmlWhl.(sprlNames{ix})(1);
				mu_tmp = cip_igf_nmlWhl.(sprlNames{ix})(2);
				lmda_tmp = cip_igf_nmlWhl.(sprlNames{ix})(3);
				
				cipConc_hybrid_igfWhl.(sprlNames{ix})(1:length(cip_binMid)) = cipConcSprl;
				cipConc_hybrid_igfWhl.(sprlNames{ix})(length(cip_binMid)+1:end) = ...
					10^n0_tmp.*cipExt_binMid(length(cip_binMid)+1:end).^mu_tmp.*exp(-lmda_tmp.*cipExt_binMid(length(cip_binMid)+1:end));
				cipConc_ext_igfWhl.(sprlNames{ix}) = 10^n0_tmp.*cipExt_binMid.^mu_tmp.*exp(-lmda_tmp.*cipExt_binMid);
				
				% Calculate extended distribution of average mass_twc over whole spiral, assuming Brown and Francis (requires D be in um)
				cipMass_hybrid_igfWhl.(sprlNames{ix})(1:length(cip_binMid)) = cipMassSprl;
				cipMass_hybrid_igfWhl.(sprlNames{ix})(length(cip_binMid)+1:end) = ...
					(a.*(cipExt_binMid(length(cip_binMid)+1:end)*10000).^b)' .* (cipConc_hybrid_igfWhl.(sprlNames{ix})(length(cip_binMid)+1:end));
				cipMass_ext_igfWhl.(sprlNames{ix}) = (a.*(cipExt_binMid*10000).^b) .* cipConc_ext_igfWhl.(sprlNames{ix});
			end
			
			if doExpLMFit
				cipConc_ext_expLMWhl.(sprlNames{ix}) = NaN(1,numExt_bins);
				cipConc_hybrid_expLMWhl.(sprlNames{ix}) = NaN(1,numExt_bins);
				cipMass_ext_expLMWhl.(sprlNames{ix}) = NaN(1,numExt_bins);
				cipMass_hybrid_expLMWhl.(sprlNames{ix}) = NaN(1,numExt_bins);
				
				[cip_expLM_nlWhl.(sprlNames{ix}), cipPSDobsFitWhl.(sprlNames{ix}), fitresultWhl,...
					gof_expLMWhl.(sprlNames{ix})] = robustLMExpFit(cip_binMid(cipIncld), cipConcSprl(cipIncld)');
				
				n0_tmp = cip_expLM_nlWhl.(sprlNames{ix})(1);
				lmda_tmp = cip_expLM_nlWhl.(sprlNames{ix})(2);
				
				cipConc_hybrid_expLMWhl.(sprlNames{ix})(1:length(cip_binMid)) = cipConcSprl;
				cipConc_hybrid_expLMWhl.(sprlNames{ix})(length(cip_binMid)+1:end) = n0_tmp.*exp(-lmda_tmp.*cipExt_binMid(length(cip_binMid)+1:end));
				cipConc_ext_expLMWhl.(sprlNames{ix}) = n0_tmp.*exp(-lmda_tmp.*cipExt_binMid);
				
				% Calculate extended distribution of average mass_twc over whole spiral, assuming Brown and Francis (requires D be in um)
				cipMass_hybrid_expLMWhl.(sprlNames{ix})(1:length(cip_binMid)) = cipMassSprl;
				cipMass_hybrid_expLMWhl.(sprlNames{ix})(length(cip_binMid)+1:end) = ...
					(a.*(cipExt_binMid(length(cip_binMid)+1:end)*10000).^b)' .* (cipConc_hybrid_expLMWhl.(sprlNames{ix})(length(cip_binMid)+1:end));
				cipMass_ext_expLMWhl.(sprlNames{ix}) = (a.*(cipExt_binMid*10000).^b) .* cipConc_ext_expLMWhl.(sprlNames{ix});
			end
		end
		
		if doTempBins
			tempCsprl.(sprlNames{ix}) = NaN(length(tempC_sprlOrig),1);
			for ii=1:length(tempC_sprlOrig)
				if tempC_sprlOrig(ii) < 0
					tempCsprl.(sprlNames{ix})(ii) = ceil(tempC_sprlOrig(ii));
				elseif tempC_sprlOrig(ii) > 0
					tempCsprl.(sprlNames{ix})(ii) = floor(tempC_sprlOrig(ii));
				end
			end
			
			tempBinsAll.(sprlNames{ix}) = min(tempCsprl.(sprlNames{ix})):max(tempCsprl.(sprlNames{ix}));
			tempBins = tempBinsAll.(sprlNames{ix});
			
			if doIGF
				cipConc_ext_igfTb.(sprlNames{ix}) = NaN(length(tempBins),numExt_bins); % N(D) using only the IGF as a basis
				cipConc_hybrid_igfTb.(sprlNames{ix}) = NaN(length(tempBins),numExt_bins); % Will contain observed CIP data, with IGF values beyond D=2mm
				cipMass_ext_igfTb.(sprlNames{ix}) = NaN(length(tempBins),numExt_bins);
				cipMass_hybrid_igfTb.(sprlNames{ix}) = NaN(length(tempBins),numExt_bins);
			end
			if doExpLMFit
				cipConc_ext_expLMTb.(sprlNames{ix}) = NaN(length(tempBins),numExt_bins); % N(D) using only the L-M exponential fit as a basis
				cipConc_hybrid_expLMTb.(sprlNames{ix}) = NaN(length(tempBins),numExt_bins);
				cipMass_ext_expLMTb.(sprlNames{ix}) = NaN(length(tempBins),numExt_bins);
				cipMass_hybrid_expLMTb.(sprlNames{ix}) = NaN(length(tempBins),numExt_bins);
			end
			
			fitSkipTmp.(sprlNames{ix}) = [];
			negLmdaTmp.(sprlNames{ix}) = [];
			for tmp = 1:length(tempBins)
				cipObsConc = nanmean(cipConc(tempCsprl.(sprlNames{ix}) == tempBins(tmp),:),1);
				cipObsMass = nanmean(cipMass(tempCsprl.(sprlNames{ix}) == tempBins(tmp),:),1);
				fprintf('\tFitting %d C (MaxT = %d C)\n',tempBins(tmp),tempBins(end))
				
				
				if doIGF
					if (all(isnan(cipObsConc)) || length(find(cipObsConc > 0)) < numValid)
						fprintf('\tInsufficient data\n')
						cip_igf_nmlTb.(sprlNames{ix})(tmp,:) = [NaN,NaN,NaN];
						cip_igf_mmtsObsTb.(sprlNames{ix})(tmp,:) = [NaN,NaN,NaN];
						cip_igf_mmtsFitTb.(sprlNames{ix})(tmp,:) = [NaN,NaN,NaN];
						cip_igf_chisquareTb.(sprlNames{ix})(tmp) = NaN;
						cipObs_igf_fitTb.(sprlNames{ix})(tmp,:) = NaN(1,size(cipObsConc(cipIncld),2));
						cip_igf_fitFlagTb.(sprlNames{ix})(tmp) = NaN;
						
						fitSkipTmp.(sprlNames{ix}) = [fitSkipTmp.(sprlNames{ix}); tmp]; %Temperature bin indices with too few valid points in SD
					else
						[cip_igf_nmlTb.(sprlNames{ix})(tmp,:), cip_igf_mmtsObsTb.(sprlNames{ix})(tmp,:), cip_igf_mmtsFitTb.(sprlNames{ix})(tmp,:),...
							cip_igf_chisquareTb.(sprlNames{ix})(tmp), cipObs_igf_fitTb.(sprlNames{ix})(tmp,:),...
							cip_igf_fitFlagTb.(sprlNames{ix})(tmp), cip_igf_fitDetails] = igfFit(cipObsConc(cipIncld), cip_binEdges(cipIncld(1):cipIncld(end)+1)', [0 2 3]);
					end
					
					n0_tmp = cip_igf_nmlTb.(sprlNames{ix})(tmp,1);
					mu_tmp = cip_igf_nmlTb.(sprlNames{ix})(tmp,2);
					lmda_tmp = cip_igf_nmlTb.(sprlNames{ix})(tmp,3);
					
					% We want to create the extended distributions if lambda is positive or NaN (NaN yields the obs PSD with NaNs
					% in the extended size range
					if lmda_tmp > 0 || isnan(lmda_tmp)
						cipConc_hybrid_igfTb.(sprlNames{ix})(tmp,1:length(cip_binMid)) = cipObsConc;
						cipConc_hybrid_igfTb.(sprlNames{ix})(tmp,length(cip_binMid)+1:end) = ...
							10^n0_tmp.*cipExt_binMid(length(cip_binMid)+1:end).^mu_tmp.*exp(-lmda_tmp.*cipExt_binMid(length(cip_binMid)+1:end));
						
						cipConc_ext_igfTb.(sprlNames{ix})(tmp,:) = 10^n0_tmp.*cipExt_binMid.^mu_tmp.*exp(-lmda_tmp.*cipExt_binMid);
						
						cipMass_hybrid_igfTb.(sprlNames{ix})(tmp,1:length(cip_binMid)) = cipObsMass;
						
						if tempBins(tmp) > mlBotTemp(ix) % below melting layer - assume liquid water spheres
							cipMass_hybrid_igfTb.(sprlNames{ix})(tmp,length(cip_binMid)+1:end) = ...
								((pi/6).*(cipExt_binMid(length(cip_binMid)+1:end)./10).^3)' .* (cipConc_hybrid_igfTb.(sprlNames{ix})(tmp,length(cip_binMid)+1:end));
							cipMass_ext_igfTb.(sprlNames{ix})(tmp,:) = ((pi/6).*(cipExt_binMid./10).^3)' .* cipConc_ext_igfTb.(sprlNames{ix})(tmp,:);
						else % above melting layer - assume ice
							cipMass_hybrid_igfTb.(sprlNames{ix})(tmp,length(cip_binMid)+1:end) = ...
								(a.*(cipExt_binMid(length(cip_binMid)+1:end)*10000).^b)' .* (cipConc_hybrid_igfTb.(sprlNames{ix})(tmp,length(cip_binMid)+1:end));
							cipMass_ext_igfTb.(sprlNames{ix})(tmp,:) = (a.*(cipExt_binMid*10000).^b)' .* cipConc_ext_igfTb.(sprlNames{ix})(tmp,:);
						end
					else
						fprintf('\tLambda < 0 - Skipping iteration\n')
						negLmdaIx.(sprlNames{ix}) = [negLmdaIx.(sprlNames{ix}); tmp]; %Indices of times with negative lambda values
					end
				end
				
				if doExpLMFit
					if (all(isnan(cipObsConc)) || length(find(cipObsConc > 0)) < numValid)
						cip_expLM_nlTb.(sprlNames{ix})(tmp,:) = [NaN,NaN];
						cipPSDobsFitTb.(sprlNames{ix})(tmp,:) = NaN(1,size(cipObsConc(cipIncld),2));
						gof_expLMTb.(sprlNames{ix})(tmp) = struct('sse',{NaN},'rsquare',{NaN},'dfe',{NaN},'adjrsquare',{NaN},'rmse',{NaN});
					else
						[cip_expLM_nlTb.(sprlNames{ix})(tmp,:), cipPSDobsFitTb.(sprlNames{ix})(tmp,:),...
							fitresult, gof_expLMTb.(sprlNames{ix})(tmp)] = robustLMExpFit(cip_binMid(cipIncld), cipObsConc(cipIncld)');
					end
					
					n0_tmp = cip_expLM_nlTb.(sprlNames{ix})(tmp,1);
					lmda_tmp = cip_expLM_nlTb.(sprlNames{ix})(tmp,2);
					
					cipConc_hybrid_expLMTb.(sprlNames{ix})(tmp,1:length(cip_binMid)) = cipObsConc;
					cipConc_hybrid_expLMTb.(sprlNames{ix})(tmp,length(cip_binMid)+1:end) = ...
						n0_tmp.*exp(-lmda_tmp.*cipExt_binMid(length(cip_binMid)+1:end));
					
					cipConc_ext_expLMTb.(sprlNames{ix})(tmp,:) = n0_tmp.*exp(-lmda_tmp.*cipExt_binMid);
					
					cipMass_hybrid_expLMTb.(sprlNames{ix})(tmp,1:length(cip_binMid)) = cipObsMass;
					
					if tempBins(tmp) > mlBotTemp(ix) % below melting layer - assume liquid water spheres
						cipMass_hybrid_expLMTb.(sprlNames{ix})(tmp,length(cip_binMid)+1:end) = ...
							((pi/6).*(cipExt_binMid(length(cip_binMid)+1:end)./10).^3)' .* (cipConc_hybrid_expLMTb.(sprlNames{ix})(tmp,length(cip_binMid)+1:end));
						cipMass_ext_expLMTb.(sprlNames{ix})(tmp,:) = ((pi/6).*(cipExt_binMid./10).^3)' .* cipConc_ext_expLMTb.(sprlNames{ix})(tmp,:);
					else % above melting layer - assume ice
						cipMass_hybrid_expLMTb.(sprlNames{ix})(tmp,length(cip_binMid)+1:end) = ...
							(a.*(cipExt_binMid(length(cip_binMid)+1:end)*10000).^b)' .* (cipConc_hybrid_expLMTb.(sprlNames{ix})(tmp,length(cip_binMid)+1:end));
						cipMass_ext_expLMTb.(sprlNames{ix})(tmp,:) = (a.*(cipExt_binMid*10000).^b)' .* cipConc_ext_expLMTb.(sprlNames{ix})(tmp,:);
					end
				end
				
			end
			fitSkipTmpRatio.(sprlNames{ix}) = length(fitSkipTmp.(sprlNames{ix}))/length(tempBins);
			fprintf('For Spiral %d:\n\tFitting skipped %d/%d times due to insufficient data\n\tLambda < 0 %d/%d times\n',...
				ix,length(fitSkipTmp.(sprlNames{ix})),length(tempBins),length(negLmdaTmp.(sprlNames{ix})),length(tempBins))
		end
		
		% Create new extended CIP structs, with NaN arrays for each spiral
		if doEvryTStp
			if doIGF
				cipConc_ext_igf.(sprlNames{ix}) = NaN(size(cipConc,1),numExt_bins); % N(D) using only the IGF as a basis
				cipConc_hybrid_igf.(sprlNames{ix}) = NaN(size(cipConc,1),numExt_bins); % Will contain observed CIP data, with IGF values beyond D=2mm
				cipMass_ext_igf.(sprlNames{ix}) = NaN(size(cipMass,1),numExt_bins);
				cipMass_hybrid_igf.(sprlNames{ix}) = NaN(size(cipMass,1),numExt_bins);
			end
			if doExpLMFit
				cipConc_ext_expLM.(sprlNames{ix}) = NaN(size(cipConc,1),numExt_bins); % N(D) using only the L-M exponential fit as a basis
				cipConc_hybrid_expLM.(sprlNames{ix}) = NaN(size(cipConc,1),numExt_bins);
				cipMass_ext_expLM.(sprlNames{ix}) = NaN(size(cipMass,1),numExt_bins);
				cipMass_hybrid_expLM.(sprlNames{ix}) = NaN(size(cipMass,1),numExt_bins);
			end
			
			fitSkipIx.(sprlNames{ix}) = [];
			negLmdaIx.(sprlNames{ix}) = [];
			for time = 1:size(cipConc,1)
				fprintf('\tFitting %d/%d\n',time,size(cipConc,1))
				cipObsConc = cipConc(time,:);
				cipObsMass = cipMass(time,:);
				
				if doIGF
					if (all(isnan(cipObsConc)) || length(find(cipObsConc > 0)) < numValid)
						fprintf('\tInsufficient data\n')
						cip_igf_nml.(sprlNames{ix})(time,:) = [NaN,NaN,NaN];
						cip_igf_mmtsObs.(sprlNames{ix})(time,:) = [NaN,NaN,NaN];
						cip_igf_mmtsFit.(sprlNames{ix})(time,:) = [NaN,NaN,NaN];
						cip_igf_chisquare.(sprlNames{ix})(time) = NaN;
						cipObs_igf_fit.(sprlNames{ix})(time,:) = NaN(1,size(cipObsConc(cipIncld),2));
						cip_igf_fitFlag.(sprlNames{ix})(time) = NaN;
						
						fitSkipIx.(sprlNames{ix}) = [fitSkipIx.(sprlNames{ix}); time]; %Indices of times with too few valid points in SD
					else
						[cip_igf_nml.(sprlNames{ix})(time,:), cip_igf_mmtsObs.(sprlNames{ix})(time,:), cip_igf_mmtsFit.(sprlNames{ix})(time,:),...
							cip_igf_chisquare.(sprlNames{ix})(time), cipObs_igf_fit.(sprlNames{ix})(time,:),...
							cip_igf_fitFlag.(sprlNames{ix})(time), cip_igf_fitDetails] = igfFit(cipObsConc(cipIncld), cip_binEdges(cipIncld(1):cipIncld(end)+1)', [0 2 3]);
					end
					
					n0_tmp = cip_igf_nml.(sprlNames{ix})(time,1);
					mu_tmp = cip_igf_nml.(sprlNames{ix})(time,2);
					lmda_tmp = cip_igf_nml.(sprlNames{ix})(time,3);
					
					% We want to create the extended distributions if lambda is positive or NaN (NaN yields the obs. PSD with NaNs
					% in the extended size range
					if lmda_tmp > 0 || isnan(lmda_tmp)
						cipConc_hybrid_igf.(sprlNames{ix})(time,1:length(cip_binMid)) = cipObsConc;
						cipConc_hybrid_igf.(sprlNames{ix})(time,length(cip_binMid)+1:end) = ...
							10^n0_tmp.*cipExt_binMid(length(cip_binMid)+1:end).^mu_tmp.*exp(-lmda_tmp.*cipExt_binMid(length(cip_binMid)+1:end));
						cipConc_ext_igf.(sprlNames{ix})(time,:) = 10^n0_tmp.*cipExt_binMid.^mu_tmp.*exp(-lmda_tmp.*cipExt_binMid);

						cipMass_hybrid_igf.(sprlNames{ix})(time,1:length(cip_binMid)) = cipObsMass;
						if ((sprlUp && time < mlBotIx(ix)) || (~sprlUp && time > mlBotIx(ix))) % below melting layer - assume liquid water spheres
							cipMass_hybrid_igf.(sprlNames{ix})(time,length(cip_binMid)+1:end) = ...
								((pi/6).*(cipExt_binMid(length(cip_binMid)+1:end)./10).^3)' .* (cipConc_hybrid_igf.(sprlNames{ix})(time,length(cip_binMid)+1:end));
							cipMass_ext_igf.(sprlNames{ix})(time,:) = ((pi/6).*(cipExt_binMid./10).^3)' .* cipConc_ext_igf.(sprlNames{ix})(time,:);
						else % above melting layer - assume ice
							cipMass_hybrid_igf.(sprlNames{ix})(time,length(cip_binMid)+1:end) = ...
								(a.*(cipExt_binMid(length(cip_binMid)+1:end)*10000).^b)' .* (cipConc_hybrid_igf.(sprlNames{ix})(time,length(cip_binMid)+1:end));
							cipMass_ext_igf.(sprlNames{ix})(time,:) = (a.*(cipExt_binMid*10000).^b)' .* cipConc_ext_igf.(sprlNames{ix})(time,:);
						end
					else
						fprintf('\tLambda < 0 - Skipping iteration\n')
						negLmdaIx.(sprlNames{ix}) = [negLmdaIx.(sprlNames{ix}); time]; %Indices of times with negative lambda values
					end
				end
				
				if doExpLMFit
					if (all(isnan(cipObsConc)) || length(find(cipObsConc > 0)) < numValid)
						cip_expLM_nl.(sprlNames{ix})(time,:) = [NaN,NaN];
						cipPSDobsFit.(sprlNames{ix})(time,:) = NaN(1,size(cipObsConc(cipIncld),2));
						gof_expLM.(sprlNames{ix})(time) = struct('sse',{NaN},'rsquare',{NaN},'dfe',{NaN},'adjrsquare',{NaN},'rmse',{NaN});
					else
						[cip_expLM_nl.(sprlNames{ix})(time,:), cipPSDobsFit.(sprlNames{ix})(time,:),...
							fitresult, gof_expLM.(sprlNames{ix})(time)] = robustLMExpFit(cip_binMid(cipIncld), cipObsConc(cipIncld)');
					end
					
					n0_tmp = cip_expLM_nl.(sprlNames{ix})(time,1);
					lmda_tmp = cip_expLM_nl.(sprlNames{ix})(time,2);
					
					cipConc_hybrid_expLM.(sprlNames{ix})(time,1:length(cip_binMid)) = cipObsConc;
					cipConc_hybrid_expLM.(sprlNames{ix})(time,length(cip_binMid)+1:end) = ...
						n0_tmp.*exp(-lmda_tmp.*cipExt_binMid(length(cip_binMid)+1:end));
					cipConc_ext_expLM.(sprlNames{ix})(time,:) = n0_tmp.*exp(-lmda_tmp.*cipExt_binMid);
					
					cipMass_hybrid_expLM.(sprlNames{ix})(time,1:length(cip_binMid)) = cipObsMass;
					if ((sprlUp && time < mlBotIx(ix)) || (~sprlUp && time > mlBotIx(ix))) % below melting layer - assume liquid water spheres
						cipMass_hybrid_expLM.(sprlNames{ix})(time,length(cip_binMid)+1:end) = ...
							((pi/6).*(cipExt_binMid(length(cip_binMid)+1:end)./10).^3)' .* (cipConc_hybrid_expLM.(sprlNames{ix})(time,length(cip_binMid)+1:end));
						cipMass_ext_expLM.(sprlNames{ix})(time,:) = ((pi/6).*(cipExt_binMid./10).^3)' .* cipConc_ext_expLM.(sprlNames{ix})(time,:);
					else % above melting layer - assume ice
						cipMass_hybrid_expLM.(sprlNames{ix})(time,length(cip_binMid)+1:end) = ...
							(a.*(cipExt_binMid(length(cip_binMid)+1:end)*10000).^b)' .* (cipConc_hybrid_expLM.(sprlNames{ix})(time,length(cip_binMid)+1:end));
						cipMass_ext_expLM.(sprlNames{ix})(time,:) = (a.*(cipExt_binMid*10000).^b)' .* cipConc_ext_expLM.(sprlNames{ix})(time,:);
					end
				end
				
			end
			fitSkipRatio.(sprlNames{ix}) = length(fitSkipIx.(sprlNames{ix}))/time;
			fprintf('For Spiral %d:\n\tFitting skipped %d/%d times due to insufficient data\n\tLambda < 0 %d/%d times\n',...
				ix,length(fitSkipIx.(sprlNames{ix})),time,length(negLmdaIx.(sprlNames{ix})),time)
		end
	end
	
	if saveMat
		save([saveDir '/CIP-Fitting/' flight mFileId fileIdStr '.mat']);
	end
end

if ~doFitting && doPlotting
	load([saveDir '/CIP-Fitting/' flight mFileId fileIdStr '.mat'],'-regexp',...
		'^(?!doFitting$|doPlotting$|plotWhole$|plotTempBins$|plotEvryTStp$|plotAllGood$|plotAllSkip$|',...
		'plotAllNegLmda$|doExpLMFit$|doIGF$|doWhole$|doTempBins$|doEvryTStp$|saveMat$|mFileId$|saveFigs$|',...
		'noDisp$|Ftype$|savePath$)\w');
end


%% Plot fit results
if doPlotting
	if plotWhole
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
			if doIGF
				plot(cipExt_binMid(cipIncld(1):end)*10,cipConc_ext_igfWhl.(sprlNames{ix})(cipIncld(1):end),'b','LineWidth',2,'DisplayName','IGF');
			end
			if doExpLMFit
				plot(cipExt_binMid(cipIncld(1):end)*10,cipConc_ext_expLMWhl.(sprlNames{ix})(cipIncld(1):end),'m','LineWidth',2,'DisplayName','L-M Exp');
			end
			title(sprintf('%s - CIP N(D) Fits - Spiral %d - %d sec avg%s',flight,ix,avgTime,idStr),'FontSize',24);
			
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
				print([saveDir '/CIP-Fitting/' flight '_FitComparison-CIP_S' num2str(ix) '_Whole' fileIdStr],Ftype,'-r0')
			end
			
		end
	end
	
	if plotTempBins
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
				if doIGF
					plot(cipExt_binMid(cipIncld(1):end)*10,cipConc_ext_igfTb.(sprlNames{ix})(tmp,cipIncld(1):end),'b','LineWidth',2,'DisplayName','IGF');
				end
				if doExpLMFit
					plot(cipExt_binMid(cipIncld(1):end)*10,cipConc_ext_expLMTb.(sprlNames{ix})(tmp,cipIncld(1):end),'m','LineWidth',2,'DisplayName','L-M Exp');
				end
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
					print([saveDir '/CIP-Fitting/' flight '_FitComparison-CIP_S' num2str(ix) '_' num2str(tempBins(tmp)) 'degC' fileIdStr],Ftype,'-r0')
				end
			end
			
		end
	end
	
	if plotEvryTStp
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
				if doIGF
					plot(cipExt_binMid(cipIncld(1):end)*10,cipConc_ext_igf.(sprlNames{ix})(time,cipIncld(1):end),'b','LineWidth',2,'DisplayName','IGF');
				end
				if doExpLMFit
					plot(cipExt_binMid(cipIncld(1):end)*10,cipConc_ext_expLM.(sprlNames{ix})(time,cipIncld(1):end),'m','LineWidth',2,'DisplayName','L-M Exp');
				end
% 				title(sprintf('%s - CIP N(D) Fits - Spiral %d - #%d - %s    %.2f%cC avg T%s',flight,ix,time,pltT.str,crntTempC,char(176),idStr),'FontSize',24);
				title({sprintf('%s - CIP N(D) Fits - Spiral %d - #%d - %s - %s',flight,ix,time,pltT.str,idStr), sprintf('%.2f%cC avg Temp',crntTempC,char(176))},'FontSize',24);


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
					print([saveDir '/CIP-Fitting/' flight '_FitComparison-CIP_S' num2str(ix) '_' sprintf('%2.2d',time) fileIdStr],Ftype,'-r0')
				end
			end
			
		end
	end
	
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
			
			title(sprintf('%s - CIP N(D) Good - Spiral %d - %d sec avg%s',flight,ix,avgTime,idStr),'FontSize',24);
			
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
				print([saveDir '/CIP-Fitting/' flight '_ND-GoodAll-CIP_S' num2str(ix) fileIdStr],Ftype,'-r0')
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
			
			title(sprintf('%s - CIP M(D) Good - Spiral %d - %d sec avg%s',flight,ix,avgTime,idStr),'FontSize',24);
			
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
				print([saveDir '/CIP-Fitting/' flight '_MD-GoodAll-CIP_S' num2str(ix) fileIdStr],Ftype,'-r0')
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
			
			title({sprintf('%s - CIP N(D) Fit Skips - Spiral %d - %d sec avg%s',flight,ix,avgTime,idStr),...
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
				print([saveDir '/CIP-Fitting/' flight '_ND-FitSkipsAll-CIP_S' num2str(ix) fileIdStr],Ftype,'-r0')
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
			
			title({sprintf('%s - CIP M(D) Fit Skips - Spiral %d - %d sec avg%s',flight,ix,avgTime,idStr),...
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
				print([saveDir '/CIP-Fitting/' flight '_MD-FitSkipsAll-CIP_S' num2str(ix) fileIdStr],Ftype,'-r0')
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
			title({sprintf('%s - CIP N(D) w/%s<0 - Spiral %d - %d sec avg%s',flight,lmdaStr,ix,avgTime,idStr),...
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
				print([saveDir '/CIP-Fitting/' flight '_ND-NegLmdaAll-CIP_S' num2str(ix) fileIdStr],Ftype,'-r0')
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
			
			title({sprintf('%s - CIP M(D) w/%s<0 - Spiral %d - %d sec avg%s',flight,lmdaStr,ix,avgTime,idStr),...
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
				print([saveDir '/CIP-Fitting/' flight '_MD-NegLmdaAll-CIP_S' num2str(ix) fileIdStr],Ftype,'-r0')
			end
			
		end
	end
end