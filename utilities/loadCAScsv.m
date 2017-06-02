%function loadCAScsv(cdpTASf, keepAll)
%% Import data from text file
% Almost exactly the same functionality as loadTASinfo.m, though with probe specific variable
% names and the option to only import key variables (i.e., TAS)

%% Initialize variables
keepAll = 0;
delimiter = ',';
startRow = 38;

%% Open the text file
[fileID, errmsg] = fopen(casTASf,'r');

%% Format Spec
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Read columns of data according to format string
dataArray = textscan(fileID, formatSpec, 'HeaderLines', startRow-1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
if keepAll
	CAS_Time = dataArray{:, 1};
	CAS_Sum_of_Transit = dataArray{:, 2};
	CAS_Sum_of_Particles = dataArray{:, 3};
	CAS_Fifo_Full = dataArray{:, 4};
	CAS_Reset_Flag = dataArray{:, 5};
	CAS_Forward_Overflow = dataArray{:, 6};
	CAS_Backward_Overflow = dataArray{:, 7};
	CAS_IAC_1 = dataArray{:, 8};
	CAS_IAC_2 = dataArray{:, 9};
	CAS_IAC_3 = dataArray{:, 10};
	CAS_IAC_4 = dataArray{:, 11};
	CAS_IAC_5 = dataArray{:, 12};
	CAS_IAC_6 = dataArray{:, 13};
	CAS_IAC_7 = dataArray{:, 14};
	CAS_IAC_8 = dataArray{:, 15};
	CAS_IAC_9 = dataArray{:, 16};
	CAS_IAC_10 = dataArray{:, 17};
	CAS_IAC_11 = dataArray{:, 18};
	CAS_IAC_12 = dataArray{:, 19};
	CAS_IAC_13 = dataArray{:, 20};
	CAS_IAC_14 = dataArray{:, 21};
	CAS_IAC_15 = dataArray{:, 22};
	CAS_IAC_16 = dataArray{:, 23};
	CAS_IAC_17 = dataArray{:, 24};
	CAS_IAC_18 = dataArray{:, 25};
	CAS_IAC_19 = dataArray{:, 26};
	CAS_IAC_20 = dataArray{:, 27};
	CAS_IAC_21 = dataArray{:, 28};
	CAS_IAC_22 = dataArray{:, 29};
	CAS_IAC_23 = dataArray{:, 30};
	CAS_IAC_24 = dataArray{:, 31};
	CAS_IAC_25 = dataArray{:, 32};
	CAS_IAC_26 = dataArray{:, 33};
	CAS_IAC_27 = dataArray{:, 34};
	CAS_IAC_28 = dataArray{:, 35};
	CAS_IAC_29 = dataArray{:, 36};
	CAS_IAC_30 = dataArray{:, 37};
	CAS_IAC_31 = dataArray{:, 38};
	CAS_IAC_32 = dataArray{:, 39};
	CAS_IAC_33 = dataArray{:, 40};
	CAS_IAC_34 = dataArray{:, 41};
	CAS_IAC_35 = dataArray{:, 42};
	CAS_IAC_36 = dataArray{:, 43};
	CAS_IAC_37 = dataArray{:, 44};
	CAS_IAC_38 = dataArray{:, 45};
	CAS_IAC_39 = dataArray{:, 46};
	CAS_IAC_40 = dataArray{:, 47};
	CAS_IAC_41 = dataArray{:, 48};
	CAS_IAC_42 = dataArray{:, 49};
	CAS_IAC_43 = dataArray{:, 50};
	CAS_IAC_44 = dataArray{:, 51};
	CAS_IAC_45 = dataArray{:, 52};
	CAS_IAC_46 = dataArray{:, 53};
	CAS_IAC_47 = dataArray{:, 54};
	CAS_IAC_48 = dataArray{:, 55};
	CAS_IAC_49 = dataArray{:, 56};
	CAS_IAC_50 = dataArray{:, 57};
	CAS_IAC_51 = dataArray{:, 58};
	CAS_IAC_52 = dataArray{:, 59};
	CAS_IAC_53 = dataArray{:, 60};
	CAS_IAC_54 = dataArray{:, 61};
	CAS_IAC_55 = dataArray{:, 62};
	CAS_IAC_56 = dataArray{:, 63};
	CAS_IAC_57 = dataArray{:, 64};
	CAS_IAC_58 = dataArray{:, 65};
	CAS_IAC_59 = dataArray{:, 66};
	CAS_IAC_60 = dataArray{:, 67};
	CAS_IAC_61 = dataArray{:, 68};
	CAS_IAC_62 = dataArray{:, 69};
	CAS_IAC_63 = dataArray{:, 70};
	CAS_IAC_64 = dataArray{:, 71};
	CAS_Dynamic_Pressure = dataArray{:, 72};
	CAS_Static_Pressure = dataArray{:, 73};
	CAS_Ambient_Temp = dataArray{:, 74};
	CAS_Forward_Heat_Sink_T = dataArray{:, 75};
	CAS_Back_Heat_Sink_T = dataArray{:, 76};
	CAS_Forward_Block_T = dataArray{:, 77};
	CAS_Backward_Block_T = dataArray{:, 78};
	CAS_Photodiode_1 = dataArray{:, 79};
	CAS_Photodiode_2 = dataArray{:, 80};
	CAS_Photodiode_3 = dataArray{:, 81};
	CAS_Photodiode_4 = dataArray{:, 82};
	CAS_Qualifier_TEC_Temp = dataArray{:, 83};
	CAS_Forward_TEC_Temp = dataArray{:, 84};
	CAS_Backward_TEC_T = dataArray{:, 85};
	CAS_Qual_Heat_Sink_T = dataArray{:, 86};
	CAS_Qual_Hi_Gain_Volt = dataArray{:, 87};
	CAS_Qual_Mid_Gain_Volt = dataArray{:, 88};
	CAS_Qual_Lo_Gain_Volt = dataArray{:, 89};
	CAS_Fwd_Hi_Gain_Volt = dataArray{:, 90};
	CAS_Fwd_Mid_Gain_Volt = dataArray{:, 91};
	CAS_Fwd_Lo_Gain_Volt = dataArray{:, 92};
	CAS_Back_Hi_Gain_Volt = dataArray{:, 93};
	CAS_Back_Mid_Gain_Volt = dataArray{:, 94};
	CAS_Back_Lo_Gain_Volt = dataArray{:, 95};
	CAS_Internal_Temp = dataArray{:, 96};
	CAS_RH_ = dataArray{:, 97};
	CAS_Spare_Analog = dataArray{:, 98};
	CAS_LWC_Hotwire = dataArray{:, 99};
	CAS_LWC_Slave_Monitor = dataArray{:, 100};
	CAS_Laser_Curr_Mon_mA = dataArray{:, 101};
	CAS_Laser_Monitor = dataArray{:, 102};
	CAS_Forw_ch0 = dataArray{:, 103};
	CAS_Forw_ch1 = dataArray{:, 104};
	CAS_Forw_ch2 = dataArray{:, 105};
	CAS_Forw_ch3 = dataArray{:, 106};
	CAS_Forw_ch4 = dataArray{:, 107};
	CAS_Forw_ch5 = dataArray{:, 108};
	CAS_Forw_ch6 = dataArray{:, 109};
	CAS_Forw_ch7 = dataArray{:, 110};
	CAS_Forw_ch8 = dataArray{:, 111};
	CAS_Forw_ch9 = dataArray{:, 112};
	CAS_Forw_ch10 = dataArray{:, 113};
	CAS_Forw_ch11 = dataArray{:, 114};
	CAS_Forw_ch12 = dataArray{:, 115};
	CAS_Forw_ch13 = dataArray{:, 116};
	CAS_Forw_ch14 = dataArray{:, 117};
	CAS_Forw_ch15 = dataArray{:, 118};
	CAS_Forw_ch16 = dataArray{:, 119};
	CAS_Forw_ch17 = dataArray{:, 120};
	CAS_Forw_ch18 = dataArray{:, 121};
	CAS_Forw_ch19 = dataArray{:, 122};
	CAS_Forw_ch20 = dataArray{:, 123};
	CAS_Forw_ch21 = dataArray{:, 124};
	CAS_Forw_ch22 = dataArray{:, 125};
	CAS_Forw_ch23 = dataArray{:, 126};
	CAS_Forw_ch24 = dataArray{:, 127};
	CAS_Forw_ch25 = dataArray{:, 128};
	CAS_Forw_ch26 = dataArray{:, 129};
	CAS_Forw_ch27 = dataArray{:, 130};
	CAS_Forw_ch28 = dataArray{:, 131};
	CAS_Forw_ch29 = dataArray{:, 132};
	CAS_Back_ch0 = dataArray{:, 133};
	CAS_Back_ch1 = dataArray{:, 134};
	CAS_Back_ch2 = dataArray{:, 135};
	CAS_Back_ch3 = dataArray{:, 136};
	CAS_Back_ch4 = dataArray{:, 137};
	CAS_Back_ch5 = dataArray{:, 138};
	CAS_Back_ch6 = dataArray{:, 139};
	CAS_Back_ch7 = dataArray{:, 140};
	CAS_Back_ch8 = dataArray{:, 141};
	CAS_Back_ch9 = dataArray{:, 142};
	CAS_Back_ch10 = dataArray{:, 143};
	CAS_Back_ch11 = dataArray{:, 144};
	CAS_Back_ch12 = dataArray{:, 145};
	CAS_Back_ch13 = dataArray{:, 146};
	CAS_Back_ch14 = dataArray{:, 147};
	CAS_Back_ch15 = dataArray{:, 148};
	CAS_Back_ch16 = dataArray{:, 149};
	CAS_Back_ch17 = dataArray{:, 150};
	CAS_Back_ch18 = dataArray{:, 151};
	CAS_Back_ch19 = dataArray{:, 152};
	CAS_Back_ch20 = dataArray{:, 153};
	CAS_Back_ch21 = dataArray{:, 154};
	CAS_Back_ch22 = dataArray{:, 155};
	CAS_Back_ch23 = dataArray{:, 156};
	CAS_Back_ch24 = dataArray{:, 157};
	CAS_Back_ch25 = dataArray{:, 158};
	CAS_Back_ch26 = dataArray{:, 159};
	CAS_Back_ch27 = dataArray{:, 160};
	CAS_Back_ch28 = dataArray{:, 161};
	CAS_Back_ch29 = dataArray{:, 162};
	CAS_True_Air_Speed = dataArray{:, 163};
	CAS_Number_Conc = dataArray{:, 164};
	CAS_LWC = dataArray{:, 165};
	CAS_MVD = dataArray{:, 166};
	CAS_ED = dataArray{:, 167};
	CAS_Status = dataArray{:, 168};
else
	CAS_Time = dataArray{:, 1};
	CAS_Sum_of_Particles = dataArray{:, 3};
	CAS_True_Air_Speed = dataArray{:, 163};
	CAS_Number_Conc = dataArray{:, 164};
	CAS_LWC = dataArray{:, 165};
end


%% Clear temporary variables
clearvars casTASf delimiter startRow formatSpec fileID dataArray ans keepAll;
%end