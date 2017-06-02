% This is used to create a cell array of times (in HHMMSS)
% pertaining to a set temperature bins (edges given by tempBinEdges) array. 
% Primarily used for determining which particles to include in which temperature
% bins when making representative particle image diagrams with the use of the ImgView
% program.

clearvars

flight = '20150702';


	
switch flight
	case '20150617'
		tempBinEdges = -19.5:1.0:16.5;
	
	case '20150620'
		tempBinEdges = -19.5:1.0:20.5;
	
	case '20150701'
		tempBinEdges = -19.5:1.0:20.5;
	
	case '20150702'
		tempBinEdges = -19.5:1.0:23.5;
	
	case '20150706'
		tempBinEdges = -19.5:1.0:19.5;

    case '20150709'
        tempBinEdges = -19.5:1.0:19.5;

end

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

startT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'startT');
endT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'endT');
startFLix = nc_varget([dataPath '/' flight '_PECANparams.nc'],'startFLix');
endFLix = nc_varget([dataPath '/' flight '_PECANparams.nc'],'endFLix');
startFLix = startFLix+1; % Indices are 0-based in parameter file
endFLix = endFLix+1;

fltLvlFile = ['/Users/danstechman/GoogleDrive/PECAN-Data/FlightLevelData/Processed/' flight '_FltLvl_Processed.mat'];
savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/';

load(fltLvlFile)

for spiralNum=1:length(startT)
	%% Identify which times pertain to which temperatures
	timehhmmss_spiral = insec2hhmmss(time_secs_FL(startFLix(spiralNum):endFLix(spiralNum)));

	numBins = length(tempBinEdges)-1;
	binIxs = cell(numBins,1);
	binTemps = cell(numBins,1);

	for ix=1:numBins
		if ix < numBins
			binIxs{ix} = find((TA(startFLix(spiralNum):endFLix(spiralNum)) >= tempBinEdges(ix)) & (TA(startFLix(spiralNum):endFLix(spiralNum)) < tempBinEdges(ix+1)));
		else
			binIxs{ix} = find((TA(startFLix(spiralNum):endFLix(spiralNum)) >= tempBinEdges(ix)) & (TA(startFLix(spiralNum):endFLix(spiralNum)) <= tempBinEdges(ix+1)));
		end

		binTemps{ix} = [mean(tempBinEdges(ix:ix+1));timehhmmss_spiral(binIxs{ix})];
	end

	table = cell2table(binTemps);
	writetable(table,[savePath flight '_Spiral' num2str(spiralNum) '_tempBinnedTimes.csv']);
end