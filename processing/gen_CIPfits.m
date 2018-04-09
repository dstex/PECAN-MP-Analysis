clearvars

%% Specify various plotting/calculation parameters
flight = '20150709';

avgTime = 10;

extnd12cm = 1; %1.2cm
extnd55cm = 0; %5.5cm

if extnd12cm
	numValid = 17; % Required number of bins with valid data needed for given SD to be fit
	fileIdStr = ['_Fit-CIP_' num2str(avgTime) 'secAvg_1.2cm'];
end
if extnd55cm
	numValid = 24;
	fileIdStr = ['_Fit-CIP_' num2str(avgTime) 'secAvg_5cm'];
end

% Check ratio of TWC from the extended portion of the PSD to that from the observed portion
% If this ratio exceeds the given threshold, the extended portion of the PSD will be set to NaN
chkTWCratio = 1;
twcRatioThresh = 0.5;

cipIncld = 8:34; % Bins to use in fitting

doSprl		= 1;
doTempBins	= 0;
doEvryTStp	= 1;

saveMat	= 1;

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';


% Brown and Francis 1995 coefficients
a = 7.38e-11;
b = 1.9;

diary([dataPath 'mp-data/' flight '/sDist/' flight fileIdStr '.log'])

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


%% Load in various PECAN parameters
mlBotTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlBotTime');
mlTopTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlTopTime');
mlBotTemp = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlBotTemp');
mlTopTemp = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlTopTemp');

%% Create array of extended CIP bins
% Get increment between bin mids for CIP and then create an extended bin mids array
% spanning from the beginning of the CIP bins to the end of the PIP bins

if extnd12cm
	pip_prtlBinEdges = [2200.0	2400.0	2600.0	2800.0	3000.0	3200.0	3400.0	3600.0	3800.0	4000.0	4200.0	4400.0...
		4600.0	4800.0	5000.0	5200.0	5500.0	5800.0	6100.0	6400.0	6700.0	7000.0	7300.0 ...
		7800.0	8300.0	8800.0	9500.0	12000.0]/10000; %cm
end

if extnd55cm
	pip_prtlBinEdges = [2200.0	2400.0	2600.0	2800.0	3000.0	3200.0	3400.0	3600.0	3800.0	4000.0	4200.0	4400.0...
		4600.0	4800.0	5000.0	5200.0	5500.0	5800.0	6100.0	6400.0	6700.0	7000.0	7300.0 ...
		8000.0:1000.0:20000.0 22000.0:2000.0:36000.0 40000.0 45000.0 50000.0 55000.0]/10000; %cm
end

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


loopVctr = 1:length(sprlNames);
% loopVctr = 2;

%% Run the fitting analyses
for ix = loopVctr
	
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
	
	%% Calculate fits of spiral-averaged SDs
	if doSprl
		fprintf('\nFitting Spiral-Averaged SD - Spiral %d\n',ix);
		cipConc_ext_igfWhl.(sprlNames{ix}) = NaN(1,numExt_bins);
		cipConc_hybrid_igfWhl.(sprlNames{ix}) = NaN(1,numExt_bins);
		cipMass_ext_igfWhl.(sprlNames{ix}) = NaN(1,numExt_bins);
		cipMass_hybrid_igfWhl.(sprlNames{ix}) = NaN(1,numExt_bins);
		
		[cip_igf_nmlWhl.(sprlNames{ix}), cip_igf_mmtsObsWhl.(sprlNames{ix}), cip_igf_mmtsFitWhl.(sprlNames{ix}),...
			cip_igf_chisquareWhl.(sprlNames{ix}), cipObs_igf_fitWhl.(sprlNames{ix}), cip_igf_fitFlagWhl.(sprlNames{ix}),...
			cip_igf_fitDetailsWhl] = igfFit(cipConcSprl(cipIncld), cip_binEdges(cipIncld(1):cipIncld(end)+1)', [0 2 3], 0);
		
		n0_tmp = cip_igf_nmlWhl.(sprlNames{ix})(1);
		mu_tmp = cip_igf_nmlWhl.(sprlNames{ix})(2);
		lmda_tmp = cip_igf_nmlWhl.(sprlNames{ix})(3);
		
		cipConc_hybrid_igfWhl.(sprlNames{ix})(1:length(cip_binMid)) = cipConcSprl;
		cipConc_hybrid_igfWhl.(sprlNames{ix})(length(cip_binMid)+1:end) = ...
			10^n0_tmp.*cipExt_binMid(length(cip_binMid)+1:end).^mu_tmp.*exp(-lmda_tmp.*cipExt_binMid(length(cip_binMid)+1:end));
		cipConc_ext_igfWhl.(sprlNames{ix}) = 10^n0_tmp.*cipExt_binMid.^mu_tmp.*exp(-lmda_tmp.*cipExt_binMid);
		
		% Calculate extended distribution of average mass_twc over whole spiral, assuming Brown and Francis (requires D be in um)
		cipMass_hybrid_igfWhl.(sprlNames{ix})(1:length(cip_binMid)) = cipMassSprl; % g cm-4
		cipMass_hybrid_igfWhl.(sprlNames{ix})(length(cip_binMid)+1:end) = ...
			(a.*(cipExt_binMid(length(cip_binMid)+1:end)*10000).^b)' .* (cipConc_hybrid_igfWhl.(sprlNames{ix})(length(cip_binMid)+1:end));
		cipMass_ext_igfWhl.(sprlNames{ix}) = (a.*(cipExt_binMid*10000).^b) .* cipConc_ext_igfWhl.(sprlNames{ix});
		
		% Calculate TWC for extended distribution
		cipTWC_hybrid_igfWhl.(sprlNames{ix}) = nansum((cipMass_hybrid_igfWhl.(sprlNames{ix}).*cipExt_binwidth').*1e6,2); % g m-3
		cipTWC_ext_igfWhl.(sprlNames{ix}) = nansum((cipMass_ext_igfWhl.(sprlNames{ix}).*cipExt_binwidth').*1e6,2); % g m-3
		
		% Calculate total number concentration for extended distribution
		cipNt_hybrid_igfWhl.(sprlNames{ix}) = nansum(cipConc_hybrid_igfWhl.(sprlNames{ix}).*cipExt_binwidth',2); % cm-3
		cipNt_ext_igfWhl.(sprlNames{ix}) = nansum(cipConc_ext_igfWhl.(sprlNames{ix}).*cipExt_binwidth',2); % cm-3
		
		% Calculate median mass diameter
		cipDmm_hybrid_igfWhl.(sprlNames{ix}) = calc_mmd(cipExt_binMid,cipMass_hybrid_igfWhl.(sprlNames{ix}).*cipExt_binwidth',cipTWC_hybrid_igfWhl.(sprlNames{ix})./1e6); % cm
		cipDmm_ext_igfWhl.(sprlNames{ix}) = calc_mmd(cipExt_binMid,cipMass_ext_igfWhl.(sprlNames{ix}).*cipExt_binwidth',cipTWC_ext_igfWhl.(sprlNames{ix})./1e6);
	end
	
	%% Calculate fits of temperature-binned SDs
	if doTempBins
		fprintf('\nFitting Temperature-Binned SDs - Spiral %d\n',ix);
		tempCrnd.(sprlNames{ix}) = NaN(length(tempC_sprlOrig),1);
		for ii=1:length(tempC_sprlOrig)
			if tempC_sprlOrig(ii) < 0
				tempCrnd.(sprlNames{ix})(ii) = ceil(tempC_sprlOrig(ii));
			elseif tempC_sprlOrig(ii) > 0
				tempCrnd.(sprlNames{ix})(ii) = floor(tempC_sprlOrig(ii));
			end
		end
		
		tempBinsAll.(sprlNames{ix}) = min(tempCrnd.(sprlNames{ix})):max(tempCrnd.(sprlNames{ix}));
		tempBins = tempBinsAll.(sprlNames{ix});
		
		if doIGF
			cipConc_ext_igfTb.(sprlNames{ix}) = NaN(length(tempBins),numExt_bins); % N(D) using only the IGF as a basis
			cipConc_hybrid_igfTb.(sprlNames{ix}) = NaN(length(tempBins),numExt_bins); % Will contain observed CIP data, with IGF values beyond D=2mm
			cipMass_ext_igfTb.(sprlNames{ix}) = NaN(length(tempBins),numExt_bins);
			cipMass_hybrid_igfTb.(sprlNames{ix}) = NaN(length(tempBins),numExt_bins);
		end
		
		fitSkipTmp.(sprlNames{ix}) = [];
		negLmdaTmp.(sprlNames{ix}) = [];
		for tmp = 1:length(tempBins)
			cipObsConc = nanmean(cipConc(tempCrnd.(sprlNames{ix}) == tempBins(tmp),:),1);
			cipObsMass = nanmean(cipMass(tempCrnd.(sprlNames{ix}) == tempBins(tmp),:),1);
			fprintf('\tFitting %d C (MaxT = %d C)\n',tempBins(tmp),tempBins(end))
			
			if (all(isnan(cipObsConc)) || length(find(cipObsConc > 0)) < numValid)
				fprintf('\t\tInsufficient data\n')
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
					cip_igf_fitFlagTb.(sprlNames{ix})(tmp), cip_igf_fitDetails] = igfFit(cipObsConc(cipIncld), cip_binEdges(cipIncld(1):cipIncld(end)+1)', [0 2 3], 0);
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
				fprintf('\t\tLambda < 0 - Skipping iteration\n')
				negLmdaTmp.(sprlNames{ix}) = [negLmdaTmp.(sprlNames{ix}); tmp]; %Indices of temps with negative lambda values
			end
			
			
		end

		% Calculate TWC for extended distribution
		cipTWC_hybrid_igfTb.(sprlNames{ix}) = nansum((cipMass_hybrid_igfTb.(sprlNames{ix}).*cipExt_binwidth').*1e6,2); % g m-3
		cipTWC_ext_igfTb.(sprlNames{ix}) = nansum((cipMass_ext_igfTb.(sprlNames{ix}).*cipExt_binwidth').*1e6,2); % g m-3
		
		% Calculate total number concentration for extended distribution
		cipNt_hybrid_igfTb.(sprlNames{ix}) = nansum(cipConc_hybrid_igfTb.(sprlNames{ix}).*cipExt_binwidth',2); % cm-3
		cipNt_ext_igfTb.(sprlNames{ix}) = nansum(cipConc_ext_igfTb.(sprlNames{ix}).*cipExt_binwidth',2); % cm-3
		
		% Calculate median mass diameter
		cipDmm_hybrid_igfTb.(sprlNames{ix}) = calc_mmd(cipExt_binMid,cipMass_hybrid_igfTb.(sprlNames{ix}).*cipExt_binwidth',cipTWC_hybrid_igfTb.(sprlNames{ix})./1e6); % cm
		cipDmm_ext_igfTb.(sprlNames{ix}) = calc_mmd(cipExt_binMid,cipMass_ext_igfTb.(sprlNames{ix}).*cipExt_binwidth',cipTWC_ext_igfTb.(sprlNames{ix})./1e6);
		
		
		fitSkipTmpRatio.(sprlNames{ix}) = length(fitSkipTmp.(sprlNames{ix}))/length(tempBins);
		fprintf('\t***Fitting skipped %d/%d times due to insufficient data\n\t***Lambda < 0 %d/%d times\n',...
			length(fitSkipTmp.(sprlNames{ix})),length(tempBins),length(negLmdaTmp.(sprlNames{ix})),length(tempBins))
	end
	
	%% Calculate fits for every time step (dependent on averaging time selected)
	if doEvryTStp
		fprintf('\nFitting Individual %d-sec SDs - Spiral %d\n',avgTime,ix);

		cipConc_ext_igf.(sprlNames{ix}) = NaN(size(cipConc,1),numExt_bins); % N(D) using only the IGF as a basis
		cipConc_hybrid_igf.(sprlNames{ix}) = NaN(size(cipConc,1),numExt_bins); % Will contain observed CIP data, with IGF values beyond D=2mm
		cipMass_ext_igf.(sprlNames{ix}) = NaN(size(cipMass,1),numExt_bins);
		cipMass_hybrid_igf.(sprlNames{ix}) = NaN(size(cipMass,1),numExt_bins);
		twcRatio.(sprlNames{ix}) = NaN(size(cipMass,1),1);
		
		fitSkipIx.(sprlNames{ix}) = [];
		negLmdaIx.(sprlNames{ix}) = [];
		for time = 1:size(cipConc,1)
			fprintf('\tFitting %d/%d\n',time,size(cipConc,1))
			cipObsConc = cipConc(time,:);
			cipObsMass = cipMass(time,:);
			
			
			if (all(isnan(cipObsConc)) || length(find(cipObsConc > 0)) < numValid)
				fprintf('\t\tInsufficient data\n')
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
					cip_igf_fitFlag.(sprlNames{ix})(time), cip_igf_fitDetails] = igfFit(cipObsConc(cipIncld), cip_binEdges(cipIncld(1):cipIncld(end)+1)', [0 2 3], 0);
			end
			
			n0_tmp = cip_igf_nml.(sprlNames{ix})(time,1);
			mu_tmp = cip_igf_nml.(sprlNames{ix})(time,2);
			lmda_tmp = cip_igf_nml.(sprlNames{ix})(time,3);
			
			% We want to create the extended distributions if lambda is positive or NaN (NaN yields the obs. PSD with NaNs
			% in the extended size range
			% First if-statement below accounts for single oddball case from PECAN
			%   where a postive derived lambda value was less than 1.0, and caused
			%   the resulting fit N(D) to decrease too slowly in the extended range,
			%   which in turn allowed for an M(D) which increased in the extended range...
			if strcmp(flight,'20150709') && ix == 8
				lmdaThresh = 1.0;
			else
				lmdaThresh = 0.0;
			end
			if lmda_tmp > lmdaThresh || isnan(lmda_tmp)
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
				
				
				% Check to see if total mass in extended portion of PSD exceeds 50% of the mass in the observed portion
				% If so, we'll set the extended portion to NaNs as the fit was likely not realistic
				if chkTWCratio
					cipTWC_obsOnly = nansum(cipMass_hybrid_igf.(sprlNames{ix})(time,1:length(cip_binMid)),2);
					cipTWC_extOnly = nansum(cipMass_hybrid_igf.(sprlNames{ix})(time,length(cip_binMid)+1:end),2);
					
					if all(isnan(cipMass_hybrid_igf.(sprlNames{ix})(time,1:length(cip_binMid))))
						cipTWC_obsOnly = NaN;
					end
					if all(isnan(cipMass_hybrid_igf.(sprlNames{ix})(time,length(cip_binMid)+1:end)))
						cipTWC_extOnly = NaN;
					end
					
					twcRatio.(sprlNames{ix})(time) = cipTWC_extOnly/cipTWC_obsOnly;
					
					if twcRatio.(sprlNames{ix})(time) > twcRatioThresh
						cipConc_hybrid_igf.(sprlNames{ix})(time,length(cip_binMid)+1:end) = NaN;
						cipConc_ext_igf.(sprlNames{ix})(time,:) = NaN;
						cipMass_hybrid_igf.(sprlNames{ix})(time,length(cip_binMid)+1:end) = NaN;
						cipMass_ext_igf.(sprlNames{ix})(time,:) = NaN;
						fprintf('\t\tIGF TWC ratio > %.2f (%.2f). Extended distribution set to NaN\n',twcRatioThresh,twcRatio.(sprlNames{ix})(time));
					end
				end
				
				
				
			else
				fprintf('\t\tLambda < 0 - Skipping iteration\n')
				negLmdaIx.(sprlNames{ix}) = [negLmdaIx.(sprlNames{ix}); time]; %Indices of times with negative lambda values
			end
			
			
		end
		

		% Calculate TWC for extended distribution
		nanMassIx_hybrid_igf = find(all(isnan(cipMass_hybrid_igf.(sprlNames{ix})),2));
		nanMassIx_ext_igf = find(all(isnan(cipMass_ext_igf.(sprlNames{ix})),2));
		cipTWC_hybrid_igf.(sprlNames{ix}) = nansum((cipMass_hybrid_igf.(sprlNames{ix}).*cipExt_binwidth').*1e6,2); % g m-3
		cipTWC_ext_igf.(sprlNames{ix}) = nansum((cipMass_ext_igf.(sprlNames{ix}).*cipExt_binwidth').*1e6,2); % g m-3
		cipTWC_hybrid_igf.(sprlNames{ix})(nanMassIx_hybrid_igf) = NaN; % Set times with NaN in all contributing mass bins to NaN 
		cipTWC_ext_igf.(sprlNames{ix})(nanMassIx_ext_igf) = NaN;
		
		% Calculate total number concentration for extended distribution
		nanConcIx_hybrid_igf = find(all(isnan(cipConc_hybrid_igf.(sprlNames{ix})),2));
		nanConcIx_ext_igf = find(all(isnan(cipConc_ext_igf.(sprlNames{ix})),2));
		cipNt_hybrid_igf.(sprlNames{ix}) = nansum(cipConc_hybrid_igf.(sprlNames{ix}).*cipExt_binwidth',2); % cm-3
		cipNt_ext_igf.(sprlNames{ix}) = nansum(cipConc_ext_igf.(sprlNames{ix}).*cipExt_binwidth',2); % cm-3
		cipNt_hybrid_igf.(sprlNames{ix})(nanConcIx_hybrid_igf) = NaN; % Set times with NaN in all contributing concentration bins to NaN
		cipNt_ext_igf.(sprlNames{ix})(nanConcIx_ext_igf) = NaN;
		
		% Calculate median mass diameter
		cipDmm_hybrid_igf.(sprlNames{ix}) = calc_mmd(cipExt_binMid,cipMass_hybrid_igf.(sprlNames{ix}).*cipExt_binwidth',cipTWC_hybrid_igf.(sprlNames{ix})./1e6); % cm
		cipDmm_ext_igf.(sprlNames{ix}) = calc_mmd(cipExt_binMid,cipMass_ext_igf.(sprlNames{ix}).*cipExt_binwidth',cipTWC_ext_igf.(sprlNames{ix})./1e6);
		cipDmm_hybrid_igf.(sprlNames{ix})(nanMassIx_hybrid_igf) = NaN; % Set times with NaN in all contributing mass bins to NaN
		cipDmm_ext_igf.(sprlNames{ix})(nanMassIx_ext_igf) = NaN;
		
		
		fitSkipRatio.(sprlNames{ix}) = length(fitSkipIx.(sprlNames{ix}))/time;
		fprintf('\t***Fitting skipped %d/%d times due to insufficient data\n\t***Lambda < 0 %d/%d times\n',...
			length(fitSkipIx.(sprlNames{ix})),time,length(negLmdaIx.(sprlNames{ix})),time)
	end
end

if saveMat
	clearvars('chkTWCratio', 'cipConc', 'cipConcSprl', 'cipDataF', 'cipMass', 'cipMassSprl', 'cipObsConc', 'cipObsMass',... 
		'cipTWC_extOnly', 'cipTWC_obsOnly', 'doEvryTStp', 'doSprl', 'doTempBins', 'ii', 'ix', 'iz',...
		'lmda_tmp', 'mu_tmp', 'n0_tmp','saveMat', 'sprlUp', 'tempBins', 'tempC_sprlOrig', 'time', 'tmp',...
		'nanConcIx_ext_igf','nanConcIx_hybrid_igf','nanMassIx_ext_igf','nanMassIx_hybrid_igf');
	save([dataPath 'mp-data/' flight '/sDist/' flight fileIdStr '.mat'],'-regexp','^(?!dataPath$|fileIdStr$)\w');
end

diary off