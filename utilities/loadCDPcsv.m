%function loadCDPcsv(cdpTASf, keepAll)
%% Import data from text file
% Almost exactly the same functionality as loadTASinfo.m, though with probe specific variable
% names and the option to only import key variables (i.e., TAS)

%% Initialize variables
keepAll = 0;
delimiter = ',';
startRow = 124;

%% Open the text file
[fileID, errmsg] = fopen(cdpTASf,'r');

%% Format Spec
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Read columns of data according to format string
dataArray = textscan(fileID, formatSpec, 'HeaderLines', startRow-1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
if keepAll
	CDP_Time = dataArray{:, 1};
	CDP_Laser_Current = dataArray{:, 2};
	CDP_Laser_Monitor = dataArray{:, 3};
	CDP_Wing_Board_Temp = dataArray{:, 4};
	CDP_Laser_Temp = dataArray{:, 5};
	CDP_HK4 = dataArray{:, 6};
	CDP_HK5 = dataArray{:, 7};
	CDP_HK6 = dataArray{:, 8};
	CDP_HK7 = dataArray{:, 9};
	CDP_DOF_Reject_Cnt = dataArray{:, 10};
	CDP_Unused = dataArray{:, 11};
	CDP_Unused1 = dataArray{:, 12};
	CDP_Unused2 = dataArray{:, 13};
	CDP_Unused3 = dataArray{:, 14};
	CDP_Unused4 = dataArray{:, 15};
	CDP_Bin_1 = dataArray{:, 16};
	CDP_Bin_2 = dataArray{:, 17};
	CDP_Bin_3 = dataArray{:, 18};
	CDP_Bin_4 = dataArray{:, 19};
	CDP_Bin_5 = dataArray{:, 20};
	CDP_Bin_6 = dataArray{:, 21};
	CDP_Bin_7 = dataArray{:, 22};
	CDP_Bin_8 = dataArray{:, 23};
	CDP_Bin_9 = dataArray{:, 24};
	CDP_Bin_10 = dataArray{:, 25};
	CDP_Bin_11 = dataArray{:, 26};
	CDP_Bin_12 = dataArray{:, 27};
	CDP_Bin_13 = dataArray{:, 28};
	CDP_Bin_14 = dataArray{:, 29};
	CDP_Bin_15 = dataArray{:, 30};
	CDP_Bin_16 = dataArray{:, 31};
	CDP_Bin_17 = dataArray{:, 32};
	CDP_Bin_18 = dataArray{:, 33};
	CDP_Bin_19 = dataArray{:, 34};
	CDP_Bin_20 = dataArray{:, 35};
	CDP_Bin_21 = dataArray{:, 36};
	CDP_Bin_22 = dataArray{:, 37};
	CDP_Bin_23 = dataArray{:, 38};
	CDP_Bin_24 = dataArray{:, 39};
	CDP_Bin_25 = dataArray{:, 40};
	CDP_Bin_26 = dataArray{:, 41};
	CDP_Bin_27 = dataArray{:, 42};
	CDP_Bin_28 = dataArray{:, 43};
	CDP_Bin_29 = dataArray{:, 44};
	CDP_Bin_30 = dataArray{:, 45};
	CDP_Number_Conc = dataArray{:, 46};
	CDP_LWC = dataArray{:, 47};
	CDP_MVD = dataArray{:, 48};
	CDP_ED = dataArray{:, 49};
	CDP_Status = dataArray{:, 50};
else
	CDP_Time = dataArray{:, 1};
	CDP_Laser_Current = dataArray{:, 2};
	CDP_Number_Conc = dataArray{:, 46};
	CDP_LWC = dataArray{:, 47};
end


%% Clear temporary variables
clearvars cdpTASf delimiter startRow formatSpec fileID dataArray ans keepAll;
%end