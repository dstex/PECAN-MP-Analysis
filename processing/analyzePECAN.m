clearvars;

%% Specify various plotting/calculation parameters
flights = {'20150617','20150620','20150701','20150702','20150706','20150709'};
% flights = {'20150617'};

avgTime = 10;

fileIdStr = ['_Fit-CIP_' num2str(avgTime) 'secAvg_1.2cm'];
% fileIdStr = ['_Fit-CIP_' num2str(avgTime) 'secAvg_5cm'];

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';


% Initialize variables to be analyzed over all flights/spirals
cipAR_sprlAvgs = NaN(42,34); % 42 spirals, with 34 diameter bins each

iSprl = 1; % Counter for total number of spirals (used in whole-project variable concat)




for iFlt = 1:length(flights)
	flight = flights{iFlt};
	%% Load in struct of all original SD data and the CIP Fitted Data
	sDistFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.CIP.10secAvg.mat'];
	sDistF = load(sDistFile);

	load([dataPath 'mp-data/' flight '/sDist/' flight fileIdStr '.mat']);

	
	loopVctr = 1:length(sprlNames);
	
	%% Concatenate multi-flight variables
	for ix = loopVctr
		if ~isnan(mlTopTime(ix))
			[~, topIx] = min(abs(sDistF.time_secs_avg.(sprlNames{ix}) - mlTopTime(ix))); % ix of nearest time matching ML top
		else
			[~, topIx] = min(abs(sDistF.tempC_avg.(sprlNames{ix}))); % ix of nearest temp nearest 0 deg C
		end
		
		% Calculate diameter-averaged area ratio for temps < 0
		if sDistF.tempC_avg.(sprlNames{ix})(1) < sDistF.tempC_avg.(sprlNames{ix})(end) % Spiral down
			cipAR_sprlAvgs(iSprl,:) = nanmean(sDistF.mean_areaRatio_avg.(sprlNames{ix})(1:topIx-1,:),1);
		else % Spiral up
			cipAR_sprlAvgs(iSprl,:) = nanmean(sDistF.mean_areaRatio_avg.(sprlNames{ix})(topIx+1:end,:),1);
		end
		
		iSprl = iSprl+1;
	end
end