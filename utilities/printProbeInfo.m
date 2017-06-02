% Used for getting diode voltages and air temperature at a specific time
% Specifically used in diagnosing PIP issues
clearvars

flight = '20150706';

probe = 'PIP';

timeInHHMMSS = 60100;
timeInSec = hhmmss2insec(timeInHHMMSS);

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

fltLvlFile = ['/Users/danstechman/GoogleDrive/PECAN-Data/FlightLevelData/Processed/' flight '_FltLvl_Processed.mat'];
flData = load(fltLvlFile);
TA = flData.TA;
timeSecsFL = flData.time_secs_FL;

switch flight
	case '20150617'
		csvCIPStr = '01CIP20150617005743.csv';
		csvPIPStr = '00PIP20150617005743.csv';
	case '20150620'
		csvCIPStr = '01CIP20150620005009.csv';
		csvPIPStr = '00PIP20150620005009.csv';
	case '20150701'
		csvCIPStr = '01CIP20150701034344.csv';
		csvPIPStr = '00PIP20150701034344.csv';
	case '20150702'
		csvCIPStr = '01CIP20150702015858.csv';
		csvPIPStr = '00PIP20150702015858.csv';
	case '20150706'
		csvCIPStr = '01CIP20150706003806.csv';
		csvPIPStr = '00PIP20150706003806.csv';
	case '20150709'
		csvCIPStr = '01CIP20150709002356.csv';
		csvPIPStr = '00PIP20150709002356.csv';
end

if strcmp(probe,'CIP')
	tasfilename = ['/Users/danstechman/GoogleDrive/PECAN-Data/mp-data/' flight '/' csvCIPStr];
else
	tasfilename = ['/Users/danstechman/GoogleDrive/PECAN-Data/mp-data/' flight '/' csvPIPStr];
end

loadTASinfo

[~, flIx] = min(abs(timeSecsFL - timeInSec));
[~, tasIx] = min(abs(Time - timeInSec));

fprintf('Searching for time: %.3f (sec)\t%s (hhmmss)\n',timeInSec,datestr(timeInSec/3600/24,'HH:MM:SS.FFF'));

fprintf('D1 voltage\t = %.2f (volts)\n',Diode_1_Volts(tasIx));
fprintf('D32 voltage\t = %.2f (volts)\n',Diode_32_Volts(tasIx));
fprintf('D64 voltage\t = %.2f (volts)\n',Diode_64_Volts(tasIx));
fprintf('Laser current\t = %.2f (mA)\n',Laser_Current(tasIx));
fprintf('Temperature\t = %.2f (%cC)\n',TA(flIx),char(176));