clearvars;

flights = {'20150709'};%,'20150620','20150701','20150702','20150706','20150709'};

fileIdStr = '.CIP.10secAvg.mat';
dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

meanAspcts = [];
for iFlt = 1:length(flights)
	flight = flights{iFlt};
	load([dataPath 'mp-data/' flight '/sDist/sdistCI.' flight fileIdStr],'mean_aspectRatio_elps_avg');
	load([dataPath 'mp-data/' flight '/sDist/sdistCI.' flight fileIdStr],'ice_flag_avg');
	sprlNames = fieldnames(mean_aspectRatio_elps_avg);
	
	for ix=1:length(sprlNames)
		meanAspcts = [meanAspcts; nanmean(nanmean(mean_aspectRatio_elps_avg.(sprlNames{ix})(~ice_flag_avg.(sprlNames{ix}),:),1))];
	end
end
nanmean(meanAspcts)