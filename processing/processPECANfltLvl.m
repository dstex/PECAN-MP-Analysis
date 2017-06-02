% This script will apply any necessary corrections to the PECAN flight-level
% data and output a new file to be used in future plotting and calculations.
%
% Dewpoint and altitude variables considered most reliable in 
% NOAA Flight QC Summaries are declared on a per-flight basis.

clear all; close all;

flight = '20150617';

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

latSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_latSens');
lonSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_lonSens');
tempSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_tempSens');
rhSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_rhSens');
dewPtSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_dwptSens');
altSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_altSens');
FLstr = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_rawFile');

if strcmp(flight,'20150706')	
	badRH = 21898:22141;
	badTd = 21898:22141;
end
	
fltLvlFile = [dataPath 'FlightLevelData/' FLstr];

saveDir = [dataPath 'FlightLevelData/Processed/'];
if (exist(saveDir, 'dir') ~= 7)
	mkdir(saveDir)
end

Torig = nc_varget(fltLvlFile,tempSens); %Ambient Temperature (C)
TDorig = nc_varget(fltLvlFile,dewPtSens);
Hum_Rel_orig = nc_varget(fltLvlFile,rhSens);
Alt = nc_varget(fltLvlFile,altSens);
HH = nc_varget(fltLvlFile,'HH');
MM = nc_varget(fltLvlFile,'MM');
SS = nc_varget(fltLvlFile,'SS');
lat = nc_varget(fltLvlFile,latSens);
lon = nc_varget(fltLvlFile,lonSens);

time_secs_FL = (HH*3600 + MM*60 + SS);

TA = Torig;
TD = TDorig;

%% Remove any bad data
if exist('badRH','var')
	Hum_Rel_orig(badRH) = NaN;
end

if exist('badTd','var')
	TD(badTd) = NaN;
end

lat(lat<30 | lat>50) = NaN;
lon(lon<-110 | lon>-85) = NaN;
Alt(Alt<0 | Alt>10000) = NaN;

%% Temperture sensor wetting correction
% Following Zipser et al. (1981)

dewExcd = find(TD > TA);

if ~isempty(dewExcd)
	newSatValue = ((TD(dewExcd) - TA(dewExcd))./2) + TA(dewExcd); % Value at midpoint between dewpoint and temperature values
	TA(dewExcd) = newSatValue;
	TD(dewExcd) = newSatValue;
end


%% Relative Humidity Calculation
% From Lowe and Ficke (1974) and Lowe (1981 pers. comm.)

% For water:
a0w = 6.107799961;
a1w = 4.436518521e-1;
a2w = 1.428945805e-2;
a3w = 2.650648471e-4;
a4w = 3.031240396e-6;
a5w = 2.034080948e-8;
a6w = 6.136820929e-11;

% For ice:
a0i = 6.10690449;
a1i = 5.02660639e-1;
a2i = 1.87743264e-2;
a3i = 4.13476180e-4;
a4i = 5.72333773e-6;
a5i = 4.71651246e-8;
a6i = 1.78086695e-10;

vp_lw = a0w+TD.*(a1w+TD.*(a2w+TD.*(a3w+TD.*(a4w+TD.*(a5w+(TD.*a6w))))));
svp_ice = a0i+TA.*(a1i+TA.*(a2i+TA.*(a3i+TA.*(a4i+TA.*(a5i+(TA.*a6i))))));
svp_lw = a0w+TA.*(a1w+TA.*(a2w+TA.*(a3w+TA.*(a4w+TA.*(a5w+(TA.*a6w))))));

RH_ice = (vp_lw./svp_ice).*100;
RH_lw = (vp_lw./svp_lw).*100;

RH_hybrid = zeros(size(Hum_Rel_orig));
for ix=1:length(RH_hybrid)
	if TA(ix) <= 0
		RH_hybrid(ix) = RH_ice(ix);
	else
		RH_hybrid(ix) = RH_lw(ix);
	end
end

% Create netCDF file and initialize dimensions
setpref('SNCTOOLS','PRESERVE_FVD',true);
ncFile = ([saveDir flight '_FltLvl_Processed.nc']);
nc_create_empty(ncFile,'netcdf4-classic');
nc_add_dimension(ncFile,'time_secs_FL',0);

% Create variable names, types, and fill values
timeStrct.Name = 'time_secs_FL';
timeStrct.Dimension = {'time_secs_FL'};
timeStrct.Datatype = 'float';
timeStrct.Attribute.Name = '_FillValue';
timeStrct.Attribute.Value = NaN;

latStrct.Name = 'lat';
latStrct.Dimension = {'time_secs_FL'};
latStrct.Datatype = 'float';
latStrct.Attribute.Name = '_FillValue';
latStrct.Attribute.Value = NaN;

lonStrct.Name = 'lon';
lonStrct.Dimension = {'time_secs_FL'};
lonStrct.Datatype = 'float';
lonStrct.Attribute.Name = '_FillValue';
lonStrct.Attribute.Value = NaN;

altStrct.Name = 'Alt';
altStrct.Dimension = {'time_secs_FL'};
altStrct.Datatype = 'float';
altStrct.Attribute.Name = '_FillValue';
altStrct.Attribute.Value = NaN;

TAStrct.Name = 'TA';
TAStrct.Dimension = {'time_secs_FL'};
TAStrct.Datatype = 'float';
TAStrct.Attribute.Name = '_FillValue';
TAStrct.Attribute.Value = NaN;

TDStrct.Name = 'TD';
TDStrct.Dimension = {'time_secs_FL'};
TDStrct.Datatype = 'float';
TDStrct.Attribute.Name = '_FillValue';
TDStrct.Attribute.Value = NaN;

RHlStrct.Name = 'RH_lw';
RHlStrct.Dimension = {'time_secs_FL'};
RHlStrct.Datatype = 'float';
RHlStrct.Attribute.Name = '_FillValue';
RHlStrct.Attribute.Value = NaN;

RHiStrct.Name = 'RH_ice';
RHiStrct.Dimension = {'time_secs_FL'};
RHiStrct.Datatype = 'float';
RHiStrct.Attribute.Name = '_FillValue';
RHiStrct.Attribute.Value = NaN;

RHhStrct.Name = 'RH_hybrid';
RHhStrct.Dimension = {'time_secs_FL'};
RHhStrct.Datatype = 'float';
RHhStrct.Attribute.Name = '_FillValue';
RHhStrct.Attribute.Value = NaN;

RHrStrct.Name = 'RH_raw';
RHrStrct.Dimension = {'time_secs_FL'};
RHrStrct.Datatype = 'float';
RHrStrct.Attribute.Name = '_FillValue';
RHrStrct.Attribute.Value = NaN;

TArStrct.Name = 'TA_raw';
TArStrct.Dimension = {'time_secs_FL'};
TArStrct.Datatype = 'float';
TArStrct.Attribute.Name = '_FillValue';
TArStrct.Attribute.Value = NaN;

TDrStrct.Name = 'TD_raw';
TDrStrct.Dimension = {'time_secs_FL'};
TDrStrct.Datatype = 'float';
TDrStrct.Attribute.Name = '_FillValue';
TDrStrct.Attribute.Value = NaN;

% Add each variable to the file
nc_addvar(ncFile,timeStrct);
nc_addvar(ncFile,latStrct);
nc_addvar(ncFile,lonStrct);
nc_addvar(ncFile,altStrct);
nc_addvar(ncFile,TAStrct);
nc_addvar(ncFile,TDStrct);
nc_addvar(ncFile,RHlStrct);
nc_addvar(ncFile,RHiStrct);
nc_addvar(ncFile,RHhStrct);
nc_addvar(ncFile,RHrStrct);
nc_addvar(ncFile,TArStrct);
nc_addvar(ncFile,TDrStrct);

% Assign additional attributes to variables
nc_attput(ncFile,nc_global,'Flight',flight);
nc_attput(ncFile,nc_global,'ProcessDate',[datestr(datetime('now')) ' Central']);
nc_attput(ncFile,nc_global,'Original HRD file',FLstr);
nc_attput(ncFile,nc_global,'Original HRD file',FLstr);
nc_attput(ncFile,'time_secs_FL','units','sec (UTC)');
nc_attput(ncFile,'time_secs_FL','Description','Time in UTC seconds');
nc_attput(ncFile,'time_secs_FL','RawVarName','(''HH''*3600 + ''MM''*60 + ''SS'')');
nc_attput(ncFile,'lat','units','degrees');
nc_attput(ncFile,'lat','Description','Latitude');
nc_attput(ncFile,'lat','RawVarName',latSens);
nc_attput(ncFile,'lon','units','degrees');
nc_attput(ncFile,'lon','Description','Longitude');
nc_attput(ncFile,'lon','RawVarName',lonSens);
nc_attput(ncFile,'Alt','units','meters');
nc_attput(ncFile,'Alt','Description','Meters above mean sea-level');
nc_attput(ncFile,'Alt','RawVarName',altSens);
nc_attput(ncFile,'TA','units','degrees C');
nc_attput(ncFile,'TA','Description','Ambient temperature with sensor wetting correction');
nc_attput(ncFile,'TA','RawVarName',tempSens);
nc_attput(ncFile,'TD','units','degrees C');
nc_attput(ncFile,'TD','Description','Dewpoint temperature with sensor wetting correction');
nc_attput(ncFile,'TD','RawVarName',dewPtSens);
nc_attput(ncFile,'RH_lw','units','%');
nc_attput(ncFile,'RH_lw','Description','Relative humidity w.r.t. liquid water');
nc_attput(ncFile,'RH_ice','units','%');
nc_attput(ncFile,'RH_ice','Description','Relative humidity w.r.t. ice');
nc_attput(ncFile,'RH_hybrid','units','%');
nc_attput(ncFile,'RH_hybrid','Description','Relative humidity (RH_lw for T>0; RH_ice for T<=0)');
nc_attput(ncFile,'RH_raw','units','%');
nc_attput(ncFile,'RH_raw','Description','HRD calculated relative humidity');
nc_attput(ncFile,'RH_raw','RawVarName',rhSens);
nc_attput(ncFile,'TA_raw','units','degrees C');
nc_attput(ncFile,'TA_raw','Description','HRD calculated ambient temperature');
nc_attput(ncFile,'TA_raw','RawVarName',tempSens);
nc_attput(ncFile,'TD_raw','units','degrees C');
nc_attput(ncFile,'TD_raw','Description','HRD calculated dewpoint temperature');
nc_attput(ncFile,'TD_raw','RawVarName',dewPtSens);

% Write actual data to each variable
nc_varput(ncFile,'time_secs_FL',time_secs_FL);
nc_varput(ncFile,'lat',lat);
nc_varput(ncFile,'lon',lon);
nc_varput(ncFile,'Alt',Alt);
nc_varput(ncFile,'TA',TA);
nc_varput(ncFile,'TD',TD);
nc_varput(ncFile,'RH_lw',RH_lw);
nc_varput(ncFile,'RH_ice',RH_ice);
nc_varput(ncFile,'RH_hybrid',RH_hybrid);
nc_varput(ncFile,'RH_raw',Hum_Rel_orig);
nc_varput(ncFile,'TA_raw',Torig);
nc_varput(ncFile,'TD_raw',TDorig);

tmp1 = saveDir;
tmp2 = flight;

clearvars altStrct latStrct lonStrct timeStrct TAStrct TDStrct RHlStrct...
	RHiStrct RHhStrct RHrStrct TArStrct TDrStrct dewPtSens altSens HH MM SS...
	dewExcd newSatValue	a0w a1w a2w a3w a4w a5w a6w a0i a1i a2i a3i a4i a5i a6i...
	fltLvlFile ix latSens lonSens ncFile rhSens dataPath FLstr saveDir flight...
	tempSens svp_ice svp_lw vp_lw

save([tmp1 tmp2 '_FltLvl_Processed.mat']);