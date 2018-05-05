clearvars

%% Specify various plotting/calculation parameters
flights = {'20150617','20150620','20150701','20150702','20150706','20150709'};
% flights = {'20150709'};

avgTime = 10;

extnd12mm = 1; %12mm
extnd55mm = 0; %55mm

numValid = 17; % Required number of bins with valid data needed for given SD to be fit

if extnd12mm
	fileIdStr = ['_Fit-CIP_' num2str(avgTime) 'secAvg_12mm'];
	statStr = '12 mm';
	pip_prtlBinEdges = [2200.0	2400.0	2600.0	2800.0	3000.0	3200.0	3400.0	3600.0	3800.0	4000.0	4200.0	4400.0...
		4600.0	4800.0	5000.0	5200.0	5500.0	5800.0	6100.0	6400.0	6700.0	7000.0	7300.0 ...
		7800.0	8300.0	8800.0	9500.0	12000.0]/1000; % [mm]
end
if extnd55mm
	fileIdStr = ['_Fit-CIP_' num2str(avgTime) 'secAvg_55mm'];
	statStr = '55 mm';
	pip_prtlBinEdges = [2200.0	2400.0	2600.0	2800.0	3000.0	3200.0	3400.0	3600.0	3800.0	4000.0	4200.0	4400.0...
		4600.0	4800.0	5000.0	5200.0	5500.0	5800.0	6100.0	6400.0	6700.0	7000.0	7300.0 ...
		8000.0:1000.0:20000.0 22000.0:2000.0:36000.0 40000.0 45000.0 50000.0 55000.0]/1000; % [mm]
end

% Check ratio of TWC from the extended portion of the PSD to that from the observed portion
% If this ratio exceeds the given threshold, the PSD index is logged for later examination
chkTWCratio = 1;
twcRatioThresh = 7.0;

cipIncld = 8:34; % Bins to use in fitting

saveMat	= 1;

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';


% Brown and Francis 1995 coefficients
% Using adjustment factor of 1.53 for the published value of 'a' (0.00294),
% as proposed by Hogan et al. (2012) to account for the fact that
% BF95 used Dmean and not Dmax as we did for PECAN
a = 0.00192; % g cm^-b
b = 1.9;

initialVars = who;
initialVars{end+1} = 'initialVars';
initialVars{end+1} = 'iFlt';

for iFlt = 1:length(flights)
	flight = flights{iFlt};
	fprintf('\nExtending %s PSDs to %s...\n',flight,statStr);
	diaryF = [dataPath 'mp-data/' flight '/sDist/' flight fileIdStr '.log'];
	if (exist(diaryF, 'file') == 2)
		delete(diaryF);
	end
	diary(diaryF)
	
	%% Load in CIP SD file and then extract only the variables we need from it
	if avgTime == 1 || avgTime == 10
		cipDataF = load([dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.CIP.10secAvg.mat']);
	elseif avgTime == 5
		cipDataF = load([dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.CIP.5secAvg.mat']);
	end
	
	if avgTime == 1
		cip_concMinR_cm4 = cipDataF.conc_minR_cm4_orig; % [cm-4]
		cip_massTWC_gcm4 = cipDataF.mass_twc_gcm4_orig; % [g cm-4]
		cip_massIWC_gcm4 = cipDataF.mass_ice_gcm4_orig; % [g cm-4]
		cip_massLWC_gcm4 = cipDataF.mass_lw_gcm4_orig; % [g cm-4]
		cip_timeSecs = cipDataF.time_secs_orig;
		iceFlag = cipDataF.ice_flag_orig;
		sprlMeanAR_lw = cipDataF.sprlMeanAspR_lw_orig;
	else
		cip_concMinR_cm4 = cipDataF.conc_minR_cm4_avg; % [cm-4]
		cip_massTWC_gcm4 = cipDataF.mass_twc_gcm4_avg; % [g cm-4]
		cip_massIWC_gcm4 = cipDataF.mass_ice_gcm4_avg; % [g cm-4]
		cip_massLWC_gcm4 = cipDataF.mass_lw_gcm4_avg; % [g cm-4]
		cip_timeSecs = cipDataF.time_secs_avg;
		iceFlag = cipDataF.ice_flag_avg;
		sprlMeanAR_lw = cipDataF.sprlMeanAspR_lw_avg;
	end
	cip_binMin_mm = cipDataF.bin_min_mm; % mm
	cip_binMax_mm = cipDataF.bin_max_mm;
	cip_binMid_mm = cipDataF.bin_mid_mm;
	cip_binEdges_mm = [cip_binMin_mm; cip_binMax_mm(end)];
	
	cip_binMin_cm = cip_binMin_mm./10;
	cip_binMax_cm = cip_binMax_mm./10;
	cip_binMid_cm = cip_binMid_mm./10;
	cip_binEdges_cm = cip_binEdges_mm./10;
	
	numObs_bins = length(cip_binMid_cm);
	
	
	%% Load in various PECAN parameters
	mlBotTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlBotTime');
	mlTopTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlTopTime');
	mlBotTemp = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlBotTemp');
	mlTopTemp = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlTopTemp');
	
	%% Create array of extended CIP bins
	% Get increment between bin mids for CIP and then create an extended bin mids array
	% spanning from the beginning of the CIP bins to the end of the PIP bins
	
	cipExt_binEdges_mm = [cip_binEdges_mm; pip_prtlBinEdges'];
	cipExt_binMin_mm = cipExt_binEdges_mm(1:end-1);
	cipExt_binMax_mm = cipExt_binEdges_mm(2:end);
	cipExt_binMid_mm = (cipExt_binMin_mm+cipExt_binMax_mm)/2;
	cipExt_binwidth_mm = diff(cipExt_binEdges_mm);
	
	cipExt_binEdges_cm = cipExt_binEdges_mm./10;
	cipExt_binMin_cm = cipExt_binMin_mm./10;
	cipExt_binMax_cm = cipExt_binMax_mm./10;
	cipExt_binMid_cm = cipExt_binMid_mm./10;
	cipExt_binwidth_cm = cipExt_binwidth_mm./10;
	
	numExt_bins = length(cipExt_binMid_cm);
	
	sprlNames = fieldnames(cip_concMinR_cm4); % Variable used unimportant - just needs to be one of the structs
	
	
	loopVctr = 1:length(sprlNames);
	% loopVctr = 2;
	
	%% Run the fitting analyses
	for ix = loopVctr
		cipConc_cm4 = cip_concMinR_cm4.(sprlNames{ix}); % [cm-4]
		cipMassTWC_gcm4 = cip_massTWC_gcm4.(sprlNames{ix}); % [g cm-4]
		cipMassIWC_gcm4 = cip_massIWC_gcm4.(sprlNames{ix}); % [g cm-4]
		cipMassLWC_gcm4 = cip_massLWC_gcm4.(sprlNames{ix}); % [g cm-4]
		iceFlg = iceFlag.(sprlNames{ix});
		
		sprlLen = size(cipConc_cm4,1);
		
		fprintf('\nFitting Individual %d-sec SDs - Spiral %d\n',avgTime,ix);
		
		cipConc_cm4_ext_igf.(sprlNames{ix}) = NaN(sprlLen,numExt_bins); % N(D) using only the IGF as a basis
		cipConc_cm4_hybrid_igf.(sprlNames{ix}) = NaN(sprlLen,numExt_bins); % Will contain observed CIP data, with IGF values beyond D=2mm
		cipMass_gcm4_ext_igf.(sprlNames{ix}) = NaN(sprlLen,numExt_bins);
		cipMass_gcm4_hybrid_igf.(sprlNames{ix}) = NaN(sprlLen,numExt_bins);
		cipMassIWC_gcm4_hybrid_igf.(sprlNames{ix}) = NaN(sprlLen,numExt_bins);
		cipMassLWC_gcm4_hybrid_igf.(sprlNames{ix}) = NaN(sprlLen,numExt_bins);
		twcRatio.(sprlNames{ix}) = NaN(sprlLen,1);
		
		fitSkipIx.(sprlNames{ix}) = [];
		lowLmdaIx.(sprlNames{ix}) = [];
		twcRatioExcdIx.(sprlNames{ix}) = [];
		
		for time = 1:size(cipConc_cm4,1)
			fprintf('\tFitting %d/%d\n',time,sprlLen)
			cipObsConc_cm4 = cipConc_cm4(time,:);
			cipObsMassTWC_gcm4 = cipMassTWC_gcm4(time,:);
			cipObsMassIWC_gcm4 = cipMassIWC_gcm4(time,:);
			cipObsMassLWC_gcm4 = cipMassLWC_gcm4(time,:);
			
			
			if (all(isnan(cipObsConc_cm4)) || length(find(cipObsConc_cm4 > 0)) < numValid)
				fprintf('\t\tInsufficient data\n')
				cip_igf_nml.(sprlNames{ix})(time,:) = [NaN,NaN,NaN];
				cip_igf_mmtsObs.(sprlNames{ix})(time,:) = [NaN,NaN,NaN];
				cip_igf_mmtsFit.(sprlNames{ix})(time,:) = [NaN,NaN,NaN];
				cip_igf_chisquare.(sprlNames{ix})(time) = NaN;
				cipObs_igf_fit.(sprlNames{ix})(time,:) = NaN(1,size(cipObsConc_cm4(cipIncld),2));
				cip_igf_fitFlag.(sprlNames{ix})(time) = NaN;
				
				fitSkipIx.(sprlNames{ix}) = [fitSkipIx.(sprlNames{ix}); time]; %Indices of times with too few valid points in SD
			else
				[cip_igf_nml.(sprlNames{ix})(time,:), cip_igf_mmtsObs.(sprlNames{ix})(time,:), cip_igf_mmtsFit.(sprlNames{ix})(time,:),...
					cip_igf_chisquare.(sprlNames{ix})(time), cipObs_igf_fit.(sprlNames{ix})(time,:),...
					cip_igf_fitFlag.(sprlNames{ix})(time), cip_igf_fitDetails] = igfFit(cipObsConc_cm4(cipIncld), cip_binEdges_cm(cipIncld(1):cipIncld(end)+1)', [0 2 3], 0);
			end
			
			n0_tmp = cip_igf_nml.(sprlNames{ix})(time,1);
			mu_tmp = cip_igf_nml.(sprlNames{ix})(time,2);
			lmda_tmp = cip_igf_nml.(sprlNames{ix})(time,3);
			
			% We want to create the extended distributions if lambda is positive or NaN (NaN yields the observed PSD with NaNs
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
				cipConc_cm4_hybrid_igf.(sprlNames{ix})(time,1:numObs_bins) = cipObsConc_cm4;
				cipConc_cm4_hybrid_igf.(sprlNames{ix})(time,numObs_bins+1:end) = ...
					10^n0_tmp.*cipExt_binMid_cm(numObs_bins+1:end).^mu_tmp.*exp(-lmda_tmp.*cipExt_binMid_cm(numObs_bins+1:end));
				cipConc_cm4_ext_igf.(sprlNames{ix})(time,:) = 10^n0_tmp.*cipExt_binMid_cm.^mu_tmp.*exp(-lmda_tmp.*cipExt_binMid_cm);
				
				cipMass_gcm4_hybrid_igf.(sprlNames{ix})(time,1:numObs_bins) = cipObsMassTWC_gcm4;
				cipMassIWC_gcm4_hybrid_igf.(sprlNames{ix})(time,1:numObs_bins) = cipObsMassIWC_gcm4;
				cipMassLWC_gcm4_hybrid_igf.(sprlNames{ix})(time,1:numObs_bins) = cipObsMassLWC_gcm4;
				
				if iceFlg(time) % Above/in melting layer - ice
					cipMass_gcm4_hybrid_igf.(sprlNames{ix})(time,numObs_bins+1:end) = ...
						(a.*cipExt_binMid_cm(numObs_bins+1:end).^b)' .* (cipConc_cm4_hybrid_igf.(sprlNames{ix})(time,numObs_bins+1:end));
					cipMass_gcm4_ext_igf.(sprlNames{ix})(time,:) = (a.*cipExt_binMid_cm.^b)' .* cipConc_cm4_ext_igf.(sprlNames{ix})(time,:);
				else % Below melting layer - assume liquid water spheres
					cipMass_gcm4_hybrid_igf.(sprlNames{ix})(time,numObs_bins+1:end) = ...
						sprlMeanAR_lw.(sprlNames{ix}).*( ((pi/6).*cipExt_binMid_cm(numObs_bins+1:end).^3)' .* (cipConc_cm4_hybrid_igf.(sprlNames{ix})(time,numObs_bins+1:end)) );
					cipMass_gcm4_ext_igf.(sprlNames{ix})(time,:) = sprlMeanAR_lw.(sprlNames{ix}).*( ((pi/6).*cipExt_binMid_cm.^3)' .* cipConc_cm4_ext_igf.(sprlNames{ix})(time,:) );
				end
				
				cipMassIWC_gcm4_hybrid_igf.(sprlNames{ix})(time,numObs_bins+1:end) = ...
					(a.*cipExt_binMid_cm(numObs_bins+1:end).^b)' .* (cipConc_cm4_hybrid_igf.(sprlNames{ix})(time,numObs_bins+1:end));
				
				if iceFlg(time) % Use overall average for adjusting LW values in regions of ice (ultimately not used)
					cipMassLWC_gcm4_hybrid_igf.(sprlNames{ix})(time,numObs_bins+1:end) = ...
						0.7.*( ((pi/6).*cipExt_binMid_cm(numObs_bins+1:end).^3)' .* (cipConc_cm4_hybrid_igf.(sprlNames{ix})(time,numObs_bins+1:end)) );
				else
					cipMassLWC_gcm4_hybrid_igf.(sprlNames{ix})(time,numObs_bins+1:end) = ...
						sprlMeanAR_lw.(sprlNames{ix}).*( ((pi/6).*cipExt_binMid_cm(numObs_bins+1:end).^3)' .* (cipConc_cm4_hybrid_igf.(sprlNames{ix})(time,numObs_bins+1:end)) );
				end
				
				% Check to see if total mass in extended portion of PSD exceeds some percentage of the mass in the observed portion
				% If so, we'll make note of where this occurs and determine if the resultant fits are potentially problematic
				if chkTWCratio
					cipTWC_obsOnly = nansum((cipMass_gcm4_hybrid_igf.(sprlNames{ix})(time,1:numObs_bins).*cipExt_binwidth_cm(1:numObs_bins)').*1e6,2);
					cipTWC_extOnly = nansum((cipMass_gcm4_hybrid_igf.(sprlNames{ix})(time,numObs_bins+1:end).*cipExt_binwidth_cm(numObs_bins+1:end)').*1e6,2);
					
					if all(isnan(cipMass_gcm4_hybrid_igf.(sprlNames{ix})(time,1:numObs_bins)))
						cipTWC_obsOnly = NaN;
					end
					if all(isnan(cipMass_gcm4_hybrid_igf.(sprlNames{ix})(time,numObs_bins+1:end)))
						cipTWC_extOnly = NaN;
					end
					
					twcRatio.(sprlNames{ix})(time) = cipTWC_extOnly/cipTWC_obsOnly;
					
					if twcRatio.(sprlNames{ix})(time) > twcRatioThresh
						cipConc_cm4_hybrid_igf.(sprlNames{ix})(time,:) = NaN;
						cipConc_cm4_ext_igf.(sprlNames{ix})(time,:) = NaN;
						cipMass_gcm4_hybrid_igf.(sprlNames{ix})(time,:) = NaN;
						cipMass_gcm4_ext_igf.(sprlNames{ix})(time,:) = NaN;
						
						twcRatioExcdIx.(sprlNames{ix}) = [twcRatioExcdIx.(sprlNames{ix}); time]; %Indices of times with potentially poor repr. of extended mass
						
						fprintf('\t\tIGF TWC ratio > %.2f (%.2f)\n',twcRatioThresh,twcRatio.(sprlNames{ix})(time));
					end
				end
				
				
				
			else
				fprintf('\t\tLambda < %.1f - Skipping iteration\n',lmdaThresh)
				lowLmdaIx.(sprlNames{ix}) = [lowLmdaIx.(sprlNames{ix}); time]; %Indices of times with negative/low lambda values
			end
			
			
		end
		
		% Calculate TWC for extended distribution
		nanMassIx_hybrid_igf = find(all(isnan(cipMass_gcm4_hybrid_igf.(sprlNames{ix})),2));
		nanMassIx_ext_igf = find(all(isnan(cipMass_gcm4_ext_igf.(sprlNames{ix})),2));
		cipTWC_gm3_hybrid_igf.(sprlNames{ix}) = nansum((cipMass_gcm4_hybrid_igf.(sprlNames{ix}).*cipExt_binwidth_cm').*1e6,2); % g m-3
		cipTWC_gm3_ext_igf.(sprlNames{ix}) = nansum((cipMass_gcm4_ext_igf.(sprlNames{ix}).*cipExt_binwidth_cm').*1e6,2); % g m-3
		cipTWC_gm3_hybrid_igf.(sprlNames{ix})(nanMassIx_hybrid_igf) = NaN; % Set times with NaN in all contributing mass bins to NaN
		cipTWC_gm3_ext_igf.(sprlNames{ix})(nanMassIx_ext_igf) = NaN;
		
		nanMassIWCIx_hybrid_igf = find(all(isnan(cipMassIWC_gcm4_hybrid_igf.(sprlNames{ix})),2));
		cipIWC_gm3_hybrid_igf.(sprlNames{ix}) = nansum((cipMassIWC_gcm4_hybrid_igf.(sprlNames{ix}).*cipExt_binwidth_cm').*1e6,2); % g m-3
		cipIWC_gm3_hybrid_igf.(sprlNames{ix})(nanMassIWCIx_hybrid_igf) = NaN; % Set times with NaN in all contributing mass bins to NaN
		
		nanMassLWCIx_hybrid_igf = find(all(isnan(cipMassLWC_gcm4_hybrid_igf.(sprlNames{ix})),2));
		cipLWC_gm3_hybrid_igf.(sprlNames{ix}) = nansum((cipMassLWC_gcm4_hybrid_igf.(sprlNames{ix}).*cipExt_binwidth_cm').*1e6,2); % g m-3
		cipLWC_gm3_hybrid_igf.(sprlNames{ix})(nanMassLWCIx_hybrid_igf) = NaN; % Set times with NaN in all contributing mass bins to NaN
		
		% Calculate total number concentration for extended distribution
		nanConcIx_hybrid_igf = find(all(isnan(cipConc_cm4_hybrid_igf.(sprlNames{ix})),2));
		nanConcIx_ext_igf = find(all(isnan(cipConc_cm4_ext_igf.(sprlNames{ix})),2));
		cipNt_cm3_hybrid_igf.(sprlNames{ix}) = nansum(cipConc_cm4_hybrid_igf.(sprlNames{ix}).*cipExt_binwidth_cm',2); % cm-3
		cipNt_cm3_ext_igf.(sprlNames{ix}) = nansum(cipConc_cm4_ext_igf.(sprlNames{ix}).*cipExt_binwidth_cm',2); % cm-3
		cipNt_cm3_hybrid_igf.(sprlNames{ix})(nanConcIx_hybrid_igf) = NaN; % Set times with NaN in all contributing concentration bins to NaN
		cipNt_cm3_ext_igf.(sprlNames{ix})(nanConcIx_ext_igf) = NaN;
		
		% Calculate median mass diameter
		cipDmm_mm_hybrid_igf.(sprlNames{ix}) = (calc_mmd(cipExt_binMid_cm,cipMass_gcm4_hybrid_igf.(sprlNames{ix}).*cipExt_binwidth_cm',cipTWC_gm3_hybrid_igf.(sprlNames{ix})./1e6)).*10; % mm
		cipDmm_mm_ext_igf.(sprlNames{ix}) = (calc_mmd(cipExt_binMid_cm,cipMass_gcm4_ext_igf.(sprlNames{ix}).*cipExt_binwidth_cm',cipTWC_gm3_ext_igf.(sprlNames{ix})./1e6)).*10;
		cipDmm_mm_hybrid_igf.(sprlNames{ix})(nanMassIx_hybrid_igf) = NaN; % Set times with NaN in all contributing mass bins to NaN
		cipDmm_mm_ext_igf.(sprlNames{ix})(nanMassIx_ext_igf) = NaN;
		
		
		totalSkip = length(fitSkipIx.(sprlNames{ix}))+length(lowLmdaIx.(sprlNames{ix}));
		fitSkipRatio.(sprlNames{ix}) = totalSkip/time;
		fprintf('\t***Fitting skipped %d/%d times (%.3f%%) due to insufficient data\n',length(fitSkipIx.(sprlNames{ix})),...
			time,(length(fitSkipIx.(sprlNames{ix}))/time)*100)
		fprintf('\t***Lambda < %.1f %d/%d times (%.3f%%)\n',lmdaThresh,length(lowLmdaIx.(sprlNames{ix})),time,...
			(length(lowLmdaIx.(sprlNames{ix}))/time)*100)
		fprintf('\t***TWC ratio exceeded %d/%d times (%.3f%%)\n',...
			length(twcRatioExcdIx.(sprlNames{ix})),time,(length(twcRatioExcdIx.(sprlNames{ix}))/time)*100)
		fprintf('\t***Total number discarded: %d (%.3f%%)\n',...
			totalSkip,fitSkipRatio.(sprlNames{ix})*100)
	end
	
	
	if saveMat
		clearvars cipConc_cm4 cipMassIWC_gcm4 cipMassLWC_gcm4 cipMassTWC_gcm4...
			cipObsConc_cm4 cipObsMassIWC_gcm4 cipObsMassLWC_gcm4 cipObsMassTWC_gcm4 cipTWC_extOnly...
			cipTWC_obsOnly diaryF iceFlg ix lmda_tmp loopVctr mu_tmp...
			n0_tmp nanMassLWCIx_hybrid_igf nanMassIx_hybrid_igf nanMassIx_ext_igf nanMassIWCIx_hybrid_igf...
			nanConcIx_hybrid_igf nanConcIx_ext_igf sprlLen sprlMeanAspctR_lw time totalSkip
		
		save([dataPath 'mp-data/' flight '/sDist/' flight fileIdStr '.mat'],'-regexp',...
			'^(?!dataPath$|fileIdStr$|iFlt$|initialVars$|chkTWCratio$|saveMat$|flights$|extnd12mm$|extnd55mm$|statStr$)\w');
	end
	
	diary off
	clearvars('-except',initialVars{:});
end