% This script can be used to output a number of parameters for a given
% flight, probe, averaging time, and time (hhmmss UTC)
% Written by Dan Stechman
% Can be run with varying degrees of interactiveness, namely
% the selection of flight, probe and averaging times


% flight = input('Flight: ','s'); % Uncomment to prompt for flight
flight = '20150706';

% probe = input('Probe: ','s'); % Uncomment to prompt for probe
probe = 'CIP';

% avgTime = input('Averaging time: '); % Uncomment to prompt for averaging time
avgTime = 5;

InTimehhmmss = input('Time (hhmmss): ');

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

startT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'startT');
endT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'endT');

sDistFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.' probe '.' num2str(avgTime) 'secAvg.mat'];

FLfile = [dataPath 'FlightLevelData/Processed/' flight '_FltLvl_Processed.mat'];

load(sDistFile);
load(FLfile, '-regexp', '^(?!flight$|dataPath$)\w');


inTime_secs = hhmmss2insec(InTimehhmmss);

% Determine which spiral the given time occurs in
for ix=1:length(startT)
	if (startT(ix) <= inTime_secs && endT(ix) >= inTime_secs)
		spiralNum = ix;
	end
end

sprlNames = fieldnames(time_secs_orig);

% Find the index of the time closest to the input
[~, inqIx] = min(abs(time_secs_avg.(sprlNames{spiralNum}) - inTime_secs));
[~, inqIxFL] = min(abs(time_secsFL_orig.(sprlNames{spiralNum}) - inTime_secs));
[~, inqIxFL_all] = min(abs(time_secsFL_all - inTime_secs));

Nt = n_avg.(sprlNames{spiralNum})(inqIx);
TWC = twc_avg.(sprlNames{spiralNum})(inqIx);
Dmm = Dmm_twc_avg.(sprlNames{spiralNum})(inqIx);

temp = tempC_orig.(sprlNames{spiralNum})(inqIxFL);
rh = RH_orig.(sprlNames{spiralNum})(inqIxFL);
alt = alt_orig.(sprlNames{spiralNum})(inqIxFL);

latitude = lat(inqIxFL_all);
longitude = lon(inqIxFL_all);


fprintf('\nData from FL and %s probe on %s at/near %.2f (hhmmss):\n',probe,flight,InTimehhmmss);
fprintf('Averaging time: %d sec\n',avgTime);
fprintf('Time occurs during spiral %d\n',spiralNum);
fprintf('Lat = %.3f\tLon = %.3f\n',latitude,longitude);
fprintf('Altitude = %.3f m MSL\n\n',alt);

fprintf('Closest probe time to input:\t%.2f (hhmmss)\t%.2f (sec)\tProbeIX = %d\n',insec2hhmmss(time_secs_avg.(sprlNames{spiralNum})(inqIx)),...
	time_secs_avg.(sprlNames{spiralNum})(inqIx),inqIx);
fprintf('Closest FL time to input:\t%.2f (hhmmss)\t%.2f (sec)\tFLix = %d\n\n',insec2hhmmss(time_secsFL_orig.(sprlNames{spiralNum})(inqIxFL)),...
	time_secsFL_orig.(sprlNames{spiralNum})(inqIxFL),inqIxFL);

fprintf('Nt = %f cm-3\n',Nt);
fprintf('TWC = %f g m-3\n',TWC);
fprintf('Dmm = %.3f mm\n',Dmm);

fprintf('Temp = %f C\n',temp);
fprintf('RH = %f %%\n',rh);


