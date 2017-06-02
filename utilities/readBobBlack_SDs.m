% Reads in data from Bob Black's PECAN output files
% Also reads in variables from our own processed SD datasets and does any necessary conversions
% (such as unnormalizing by sample volume or binwidth) to line up with Black's data

close all;clearvars;

flight = '20150706';

avgTime = 6;

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';


%% Specify file names and read in raw Black data and UIUC structures
blackCIPFile = [dataPath 'mp-data/' flight '/Black_SDs/' flight '_BobBlack.CIP'];
uiucCIPFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.CIP.' num2str(avgTime) 'secAvg.mat'];
blackPIPFile = [dataPath 'mp-data/' flight '/Black_SDs/' flight '_BobBlack.PIP'];
uiucPIPFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.PIP.' num2str(avgTime) 'secAvg.mat'];

startT = nc_varget([dataPath flight '_PECANparams.nc'],'startT');
endT = nc_varget([dataPath flight '_PECANparams.nc'],'endT');
mlBotTime = nc_varget([dataPath flight '_PECANparams.nc'],'mlBotTime');

uiucCIPSDall = load(uiucCIPFile);
uiucPIPSDall = load(uiucPIPFile);

fidCIP = fopen(blackCIPFile,'r','ieee-be');
rawBlackCIPData = fread(fidCIP,'single');
blackCIPSDall = reshape(rawBlackCIPData(1:end-64),[64 7 floor(length(rawBlackCIPData)/64/7)]);
numRecsCIP = size(blackCIPSDall,3);

fidPIP = fopen(blackPIPFile,'r','ieee-be');
rawBlackPIPData = fread(fidPIP,'single');
blackPIPSDall = reshape(rawBlackPIPData(1:end-64),[64 7 floor(length(rawBlackPIPData)/64/7)]);
numRecsPIP = size(blackPIPSDall,3);


	
%% Read in all of the Black PIP variables
bb_pip_date_all = NaN(numRecsPIP,1); % YYMMDD
bb_pip_timeHHMMSS_all = NaN(numRecsPIP,1); %HHMMSS.m
bb_pip_tas_all = NaN(numRecsPIP,1);
bb_pip_avgTime_all = NaN(numRecsPIP,1);
bb_pip_lwc_all = NaN(numRecsPIP,1); % g/m3
bb_pip_iwc_all = NaN(numRecsPIP,1); % g/m3
bb_pip_smplVol_L_all = NaN(numRecsPIP,1); % Liters
bb_pip_numRjct_all = NaN(numRecsPIP,1);
bb_pip_sumPixelA_all = NaN(numRecsPIP,1); % mm^2
bb_pip_avgPixelA_all = NaN(numRecsPIP,1); % mm^2

bb_pip_conc_total_all = NaN(numRecsPIP,64);
bb_pip_conc_frac_all = NaN(numRecsPIP,64);

bb_pip_binEdges = 0.05:0.1:6.45;
bb_pip_binwidth = diff(bb_pip_binEdges); % mm
bb_pip_binMin = bb_pip_binEdges(1:end-1);
bb_pip_binMax = bb_pip_binEdges(2:end);
bb_pip_binMid = (bb_pip_binMin + bb_pip_binMax)/2;

for ix = 1:numRecsPIP
	bb_pip_date_all(ix) = blackPIPSDall(1,1,ix);
	bb_pip_timeHHMMSS_all(ix) = blackPIPSDall(2,1,ix);
	bb_pip_tas_all(ix) = blackPIPSDall(9,1,ix);
	bb_pip_avgTime_all(ix) = blackPIPSDall(8,1,ix);
	bb_pip_lwc_all(ix) = blackPIPSDall(10,1,ix);
	bb_pip_iwc_all(ix) = blackPIPSDall(11,1,ix);
	bb_pip_smplVol_L_all(ix) = blackPIPSDall(15,1,ix);
	bb_pip_numRjct_all(ix) = blackPIPSDall(18,1,ix);
	bb_pip_sumPixelA_all(ix) = blackPIPSDall(27,1,ix);
	bb_pip_avgPixelA_all(ix) = blackPIPSDall(28,1,ix);
	
	bb_pip_conc_total_all(ix,:) = (blackPIPSDall(:,2,ix))'; % (L-1) (0.1 mm-1)
	bb_pip_conc_frac_all(ix,:) = (blackPIPSDall(:,3,ix))';
% 	bb_pip_conc_centerIn_all(ix,:) = bb_pip_conc_total_all(ix,:) - bb_pip_conc_frac_all(ix,:);
% 	bb_pip_conc_centerIn_all(ix,:) = bb_pip_conc_frac_all(ix,:) - bb_pip_conc_total_all(ix,:);
end

% Convert PIP concentrations from L^-1 0.1mm^-1 to cm^-4
bb_pip_conc_total_all = bb_pip_conc_total_all.*0.1; % cm^-4
bb_pip_conc_frac_all = bb_pip_conc_frac_all.*0.1;
% bb_pip_conc_centerIn_all = bb_pip_conc_total_all - bb_pip_conc_frac_all;
bb_pip_conc_centerIn_all = bb_pip_conc_frac_all - bb_pip_conc_total_all;


% Convert timeHHMMSS to time in seconds
bb_pip_timeSecs_all = hhmmss2insec(bb_pip_timeHHMMSS_all);

for iz=1:length(startT)
	sprlLocs = find(bb_pip_timeSecs_all >= startT(iz) & bb_pip_timeSecs_all <= endT(iz));
	
	spiralStr = ['sprl' num2str(iz)];
	
	% Initialize spiral structs then sort Black variables into them
	bb_pip_timeSecsAvg.(spiralStr) = bb_pip_timeSecs_all(sprlLocs);
	bb_pip_tasAvg.(spiralStr) = bb_pip_tas_all(sprlLocs);
	bb_pip_lwcAvg.(spiralStr) = bb_pip_lwc_all(sprlLocs);
	bb_pip_iwcAvg.(spiralStr) = bb_pip_iwc_all(sprlLocs);
	bb_pip_smplVolAvg.(spiralStr) = bb_pip_smplVol_L_all(sprlLocs);
	bb_pip_numRjctAvg.(spiralStr) = bb_pip_numRjct_all(sprlLocs);
	bb_pip_sumPixelAreaAvg.(spiralStr) = bb_pip_sumPixelA_all(sprlLocs);
	bb_pip_avgPixelAreaAvg.(spiralStr) = bb_pip_avgPixelA_all(sprlLocs);
	
	bb_pip_concTAvg.(spiralStr) = bb_pip_conc_total_all(sprlLocs,:);
	bb_pip_concFAvg.(spiralStr) = bb_pip_conc_frac_all(sprlLocs,:);
	bb_pip_concAvg.(spiralStr) = bb_pip_conc_centerIn_all(sprlLocs,:);
	
end

%% Read in all UIUC PIP variables
uiuc_pip_timeSecsAvg = uiucPIPSDall.time_secs_avg;
uiuc_pip_lwcAvg = uiucPIPSDall.lwc_avg; %g cm-4
uiuc_pip_iwcAvg = uiucPIPSDall.iwc_avg; % g cm-4
uiuc_pip_smplVolAvg = uiucPIPSDall.sampleVol_avg; % cm-3
uiuc_pip_rjctRatioAvg = uiucPIPSDall.reject_ratio_avg;
uiuc_pip_areaAvg = uiucPIPSDall.area_avg; 
uiuc_pip_calcdAreaAvg = uiucPIPSDall.area_calcd_avg;

uiuc_pip_concAvg = uiucPIPSDall.conc_minR_avg; % cm-4

uiuc_pip_binwidth = uiucPIPSDall.bin_size; % mm
uiuc_pip_binMin = uiucPIPSDall.bin_min;
uiuc_pip_binMax = uiucPIPSDall.bin_max;
uiuc_pip_binMid = uiucPIPSDall.bin_mid;
	
sprlNames = fieldnames(uiuc_pip_concAvg);


%% Read in all of the Black CIP variables
bb_cip_date_all = NaN(numRecsCIP,1);
bb_cip_timeHHMMSS_all = NaN(numRecsCIP,1);
bb_cip_tas_all = NaN(numRecsCIP,1);
bb_cip_avgTime_all = NaN(numRecsCIP,1);
bb_cip_lwc_all = NaN(numRecsCIP,1);
bb_cip_iwc_all = NaN(numRecsCIP,1);
bb_cip_smplVol_L_all = NaN(numRecsCIP,1);
bb_cip_numRjct_all = NaN(numRecsCIP,1);
bb_cip_sumPixelA_all = NaN(numRecsCIP,1);
bb_cip_avgPixelA_all = NaN(numRecsCIP,1);

bb_cip_conc_total_all = NaN(numRecsCIP,64);
bb_cip_conc_frac_all = NaN(numRecsCIP,64);
bb_cip_conc_centerIn_all = NaN(numRecsCIP,64);

bb_cip_binEdges = 0.0125:0.025:1.6125;
bb_cip_binwidth = diff(bb_cip_binEdges); % mm
bb_cip_binMin = bb_cip_binEdges(1:end-1);
bb_cip_binMax = bb_cip_binEdges(2:end);
bb_cip_binMid = (bb_cip_binMin + bb_cip_binMax)/2;



for ix = 1:numRecsCIP
	bb_cip_date_all(ix) = blackCIPSDall(1,1,ix);
	bb_cip_timeHHMMSS_all(ix) = blackCIPSDall(2,1,ix);
	bb_cip_tas_all(ix) = blackCIPSDall(9,1,ix);
	bb_cip_avgTime_all(ix) = blackCIPSDall(8,1,ix);
	bb_cip_lwc_all(ix) = blackCIPSDall(12,1,ix); %g m-3
	bb_cip_iwc_all(ix) = blackCIPSDall(13,1,ix); %g m-3
	bb_cip_smplVol_L_all(ix) = blackCIPSDall(14,1,ix); %L
	bb_cip_numRjct_all(ix) = blackCIPSDall(18,1,ix);
	bb_cip_sumPixelA_all(ix) = blackCIPSDall(25,1,ix); %mm2
	bb_cip_avgPixelA_all(ix) = blackCIPSDall(26,1,ix); %mm2
	
	bb_cip_conc_total_all(ix,:) = (blackCIPSDall(:,4,ix))'; % (L-1) (0.025 mm-1)
	bb_cip_conc_frac_all(ix,:) = (blackCIPSDall(:,5,ix))';
end

% Convert CIP concentrations from L^-1 0.025mm^-1 to cm^-4
bb_cip_conc_total_all = bb_cip_conc_total_all.*0.4; % cm^-4
bb_cip_conc_frac_all = bb_cip_conc_frac_all.*0.4;
% bb_cip_conc_centerIn_all = bb_cip_conc_total_all - bb_cip_conc_frac_all;
bb_cip_conc_centerIn_all = bb_cip_conc_frac_all - bb_cip_conc_total_all;

% Convert timeHHMMSS to time in seconds
bb_cip_timeSecs_all = hhmmss2insec(bb_cip_timeHHMMSS_all);

for iz=1:length(startT)
	sprlLocs = find(bb_cip_timeSecs_all >= startT(iz) & bb_cip_timeSecs_all <= endT(iz));
	
	spiralStr = ['sprl' num2str(iz)];
	
	% Initialize spiral structs then sort Black variables into them
	bb_cip_timeSecsAvg.(spiralStr) = bb_cip_timeSecs_all(sprlLocs);
	bb_cip_tasAvg.(spiralStr) = bb_cip_tas_all(sprlLocs);
	bb_cip_lwcAvg.(spiralStr) = bb_cip_lwc_all(sprlLocs);
	bb_cip_iwcAvg.(spiralStr) = bb_cip_iwc_all(sprlLocs);
	bb_cip_smplVolAvg.(spiralStr) = bb_cip_smplVol_L_all(sprlLocs);
	bb_cip_numRjctAvg.(spiralStr) = bb_cip_numRjct_all(sprlLocs);
	bb_cip_sumPixelAreaAvg.(spiralStr) = bb_cip_sumPixelA_all(sprlLocs);
	bb_cip_avgPixelAreaAvg.(spiralStr) = bb_cip_avgPixelA_all(sprlLocs);
	
	bb_cip_concTAvg.(spiralStr) = bb_cip_conc_total_all(sprlLocs,:);
	bb_cip_concFAvg.(spiralStr) = bb_cip_conc_frac_all(sprlLocs,:);
	bb_cip_concAvg.(spiralStr) = bb_cip_conc_centerIn_all(sprlLocs,:);
	
end

%% Read in all UIUC CIP variables
uiuc_cip_timeSecsAvg = uiucCIPSDall.time_secs_avg;
uiuc_cip_lwcAvg = uiucCIPSDall.lwc_avg;
uiuc_cip_iwcAvg = uiucCIPSDall.iwc_avg;
uiuc_cip_smplVolAvg = uiucCIPSDall.sampleVol_avg;
uiuc_cip_rjctRatioAvg = uiucCIPSDall.reject_ratio_avg;
uiuc_cip_areaAvg = uiucCIPSDall.area_avg;
uiuc_cip_calcdAreaAvg = uiucCIPSDall.area_calcd_avg;

uiuc_cip_concAvg = uiucCIPSDall.conc_minR_avg;

uiuc_cip_binwidth = uiucCIPSDall.bin_size;
uiuc_cip_binMin = uiucCIPSDall.bin_min;
uiuc_cip_binMax = uiucCIPSDall.bin_max;
uiuc_cip_binMid = uiucCIPSDall.bin_mid;