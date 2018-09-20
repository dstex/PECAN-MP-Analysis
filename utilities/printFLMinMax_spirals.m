% Used particularly for quickly getting max/min values for spirals
% to fill in P-3 Spirals spreadsheet

flight = '20150620';

probe = 'CIP';

avgTime = 10; % Not relevant for this script as we're looking at 1-sec FL data

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

startT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'startT');
endT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'endT');

sDistFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.' probe '.' num2str(avgTime) 'secAvg.mat'];

load(sDistFile);


sprlNames = fieldnames(time_secs_orig);

for ix=1:length(startT)
	spiralNum = ix;
	time_secs = time_secsFL_orig.(sprlNames{spiralNum});
	time_hhmmss = insec2hhmmss(time_secs);

	tempC = tempC_orig.(sprlNames{spiralNum});
	RH = RH_orig.(sprlNames{spiralNum});
	alt = alt_orig.(sprlNames{spiralNum});


	fprintf('\nData from FL on %s for spiral %d\n\n',flight,spiralNum);

	fprintf('startT\t = %d (hhmmss)\t%.2f (sec)\n',time_hhmmss(1),time_secs(1));
	fprintf('endT\t = %d (hhmmss)\t%.2f (sec)\n\n',time_hhmmss(end),time_secs(end));


	fprintf('Altitude Start\t = %.4f m MSL\n',alt(1));
	fprintf('Altitude End\t = %.4f m MSL\n\n',alt(end));

	fprintf('Min Temp = %.2f C\n',min(tempC));
	fprintf('Max Temp = %.2f C\n\n',max(tempC));

	fprintf('Min hybrid RH = %.2f%%\n',min(RH));
	fprintf('Max hybrid RH = %.2f%%\n\n',max(RH));
end