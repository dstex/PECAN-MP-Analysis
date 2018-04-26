% Script for determining number of fits where TWC from extended
% distribution exceeded the given TWC threshold (usually > 50% of
% the TWC from the observed portion of the distribution).

clearvars;

flights = {'20150617','20150620','20150701','20150702','20150706','20150709'};

fileIdStr = '_Fit-CIP_10secAvg_5cm';
dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

numTWCdiscardT = 0;
numPSDsT = 0;
numAllSkipsT = 0;
numInsufDatT = 0;
numLowLmdaT = 0;
excdTWCratios = [];
allSkipRatios = [];
for iFlt = 1:length(flights)
	flight = flights{iFlt};
	fprintf('\nCounting %s...\n',flight);
	dataFlt = load([dataPath 'mp-data/' flight '/sDist/' flight fileIdStr '.mat']);
	
	numTWCdiscardFlt = 0;
	numInsufDatFlt = 0;
	numLowLmdaFlt = 0;
	numPSDsFlt = 0;
	for ix=1:length(dataFlt.sprlNames)
		temp = dataFlt.twcRatioExcdIx.(dataFlt.sprlNames{ix});
		excdTWCratios = [excdTWCratios; dataFlt.twcRatio.(dataFlt.sprlNames{ix})(temp)];
		allSkipRatios = [allSkipRatios; dataFlt.fitSkipRatio.(dataFlt.sprlNames{ix})];
		numTWCdiscardFlt=numTWCdiscardFlt+length(temp);
		numInsufDatFlt = numInsufDatFlt+length(dataFlt.fitSkipIx.(dataFlt.sprlNames{ix}));
		numLowLmdaFlt = numLowLmdaFlt+length(dataFlt.lowLmdaIx.(dataFlt.sprlNames{ix}));
		numPSDsFlt = numPSDsFlt+length(dataFlt.twcRatio.(dataFlt.sprlNames{ix}));
% 		if ~isempty(temp)
% 			fprintf('Spiral %d:\n',ix);
% 			disp(temp);
% 		end
	end
	fprintf('\nDiscarded %d/%d (%.5f) PSDs for TWC ratio exceedance\n',numTWCdiscardFlt,numPSDsFlt,(numTWCdiscardFlt/numPSDsFlt));
	numAllSkipsT = numAllSkipsT+numInsufDatFlt+numLowLmdaFlt;
	numTWCdiscardT = numTWCdiscardT+numTWCdiscardFlt;
	numInsufDatT = numInsufDatT+numInsufDatFlt;
	numLowLmdaT = numLowLmdaT+numLowLmdaFlt;
	numPSDsT = numPSDsT+numPSDsFlt;
	clearvars dataFlt
end
fprintf('\nSkipped fitting %d/%d (%.5f) PSDs TOTAL due to insufficient data',numInsufDatT,numPSDsT,(numInsufDatT/numPSDsT));
fprintf('\nSkipped fitting %d/%d (%.5f) PSDs TOTAL due to low lambda values',numLowLmdaT,numPSDsT,(numLowLmdaT/numPSDsT));
fprintf('\nDiscarded %d/%d (%.5f) PSDs TOTAL for TWC ratio exceedance',numTWCdiscardT,numPSDsT,(numTWCdiscardT/numPSDsT));
fprintf('\nSkipped fitting %d/%d (%.5f) PSDs TOTAL due to all factors\n',numAllSkipsT,numPSDsT,(numAllSkipsT/numPSDsT));