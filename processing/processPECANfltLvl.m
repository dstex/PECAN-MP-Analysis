% This script will apply any necessary corrections to the PECAN flight-level
% data and output a new file to be used in future plotting and calculations.
%
% Dewpoint and altitude variables considered most reliable in 
% NOAA Flight QC Summaries are declared on a per-flight basis.

clear all; close all;

runAll = 0;

if runAll
	flights = {'20150617','20150620','20150701','20150702','20150706','20150709'};
else
	flights = {'20150620'};
end


for iiz=1:length(flights)
	flight = flights{iiz};

	dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

	latSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_latSens');
	lonSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_lonSens');
	tempSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_tempSens');
	rhSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_rhSens');
	dewPtSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_dwptSens');
	altSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_altSens');
	spSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_spSens');
	dpSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_dpSens');
	wSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_wSens');
	wind_xSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_wind_xSens');
	wind_ySens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_wind_ySens');
	relWind_xSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_relWind_xSens');
	relWind_ySens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_relWind_ySens');
	wdSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_wdSens');
	wsSens = nc_attget([dataPath '/' flight '_PECANparams.nc'],-1,'FL_wsSens');
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
	staticP = nc_varget(fltLvlFile,spSens);
	dynamicP = nc_varget(fltLvlFile,dpSens);
	w_dpj = nc_varget(fltLvlFile,wSens);
	u = nc_varget(fltLvlFile,wind_xSens);
	v = nc_varget(fltLvlFile,wind_ySens);
	uRel = nc_varget(fltLvlFile,relWind_xSens);
	vRel = nc_varget(fltLvlFile,relWind_ySens);
	windD = nc_varget(fltLvlFile,wdSens);
	windS = nc_varget(fltLvlFile,wsSens);
	TAS = nc_varget(fltLvlFile,'TAS.d');

	Alt = nc_varget(fltLvlFile,altSens);
	HH = nc_varget(fltLvlFile,'HH');
	MM = nc_varget(fltLvlFile,'MM');
	SS = nc_varget(fltLvlFile,'SS');
	lat = nc_varget(fltLvlFile,latSens);
	lon = nc_varget(fltLvlFile,lonSens);

	time_secs_FL = (HH*3600 + MM*60 + SS);


	mDate = datetime(flight,'InputFormat','yyyyMMdd');
	mDatePre = mDate - days(1);
	mDatePreStr = datestr(mDatePre,'yyyymmdd');
	fl_dt = cell(size(HH));
	for idt=1:length(HH)
		if HH(idt) > 18
			fl_dt{idt} = sprintf('%s%02d%02d%02d',mDatePreStr,HH(idt),MM(idt),SS(idt));
		else
			fl_dt{idt} = sprintf('%s%02d%02d%02d',flight,HH(idt),MM(idt),SS(idt));
		end
	end
	fl_dt = str2double(fl_dt);


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
	
	
	% Alternate RH calculation using equations for SVP over water and ice
	% as given in Huang (2018 JAMC)
	
	vp_lw_H18 = exp(34.494 - (4924.99./(TD+237.1)))./((TD+105).^1.57);
	svp_ice_H18 = exp(43.494 - (6545.8./(TA+278)))./((TA+868).^2);
	svp_lw_H18 = exp(34.494 - (4924.99./(TA+237.1)))./((TA+105).^1.57);

	RH_ice_H18 = (vp_lw_H18./svp_ice_H18).*100;
	RH_lw_H18 = (vp_lw_H18./svp_lw_H18).*100;
	
	
	RH_hybrid = zeros(size(Hum_Rel_orig));
	RH_hybrid_H18 = zeros(size(Hum_Rel_orig));
	for ix=1:length(RH_hybrid)
		if TA(ix) <= 0
			RH_hybrid_H18(ix) = RH_ice_H18(ix);
			RH_hybrid(ix) = RH_ice(ix);
		else
			RH_hybrid_H18(ix) = RH_lw_H18(ix);
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

	dtStrct.Name = 'datetime_FL';
	dtStrct.Dimension = {'time_secs_FL'};

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

	statPStrct.Name = 'staticPres';
	statPStrct.Dimension = {'time_secs_FL'};
	statPStrct.Datatype = 'float';
	statPStrct.Attribute.Name = '_FillValue';
	statPStrct.Attribute.Value = NaN;

	dynmPStrct.Name = 'dynamicPres';
	dynmPStrct.Dimension = {'time_secs_FL'};
	dynmPStrct.Datatype = 'float';
	dynmPStrct.Attribute.Name = '_FillValue';
	dynmPStrct.Attribute.Value = NaN;

	wdpjStrct.Name = 'w_dpj';
	wdpjStrct.Dimension = {'time_secs_FL'};
	wdpjStrct.Datatype = 'float';
	wdpjStrct.Attribute.Name = '_FillValue';
	wdpjStrct.Attribute.Value = NaN;

	uStrct.Name = 'u';
	uStrct.Dimension = {'time_secs_FL'};
	uStrct.Datatype = 'float';
	uStrct.Attribute.Name = '_FillValue';
	uStrct.Attribute.Value = NaN;

	vStrct.Name = 'v';
	vStrct.Dimension = {'time_secs_FL'};
	vStrct.Datatype = 'float';
	vStrct.Attribute.Name = '_FillValue';
	vStrct.Attribute.Value = NaN;

	uRStrct.Name = 'u_relative';
	uRStrct.Dimension = {'time_secs_FL'};
	uRStrct.Datatype = 'float';
	uRStrct.Attribute.Name = '_FillValue';
	uRStrct.Attribute.Value = NaN;

	vRStrct.Name = 'v_relative';
	vRStrct.Dimension = {'time_secs_FL'};
	vRStrct.Datatype = 'float';
	vRStrct.Attribute.Name = '_FillValue';
	vRStrct.Attribute.Value = NaN;

	wdStrct.Name = 'windDir';
	wdStrct.Dimension = {'time_secs_FL'};
	wdStrct.Datatype = 'float';
	wdStrct.Attribute.Name = '_FillValue';
	wdStrct.Attribute.Value = NaN;

	wsStrct.Name = 'windSpd';
	wsStrct.Dimension = {'time_secs_FL'};
	wsStrct.Datatype = 'float';
	wsStrct.Attribute.Name = '_FillValue';
	wsStrct.Attribute.Value = NaN;

	tasStrct.Name = 'TAS';
	tasStrct.Dimension = {'time_secs_FL'};
	tasStrct.Datatype = 'float';
	tasStrct.Attribute.Name = '_FillValue';
	tasStrct.Attribute.Value = NaN;

	% Add each variable to the file
	nc_addvar(ncFile,timeStrct);
	nc_addvar(ncFile,dtStrct);
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
	nc_addvar(ncFile,statPStrct);
	nc_addvar(ncFile,dynmPStrct);
	nc_addvar(ncFile,wdpjStrct);
	nc_addvar(ncFile,uStrct);
	nc_addvar(ncFile,vStrct);
	nc_addvar(ncFile,uRStrct);
	nc_addvar(ncFile,vRStrct);
	nc_addvar(ncFile,wdStrct);
	nc_addvar(ncFile,wsStrct);
	nc_addvar(ncFile,tasStrct);

	% Assign additional attributes to variables
	nc_attput(ncFile,nc_global,'Flight',flight);
	nc_attput(ncFile,nc_global,'ProcessDate',[datestr(datetime('now')) ' Central']);
	nc_attput(ncFile,nc_global,'Original HRD file',FLstr);
	nc_attput(ncFile,nc_global,'Original HRD file',FLstr);

	nc_attput(ncFile,'time_secs_FL','units','sec (UTC)');
	nc_attput(ncFile,'time_secs_FL','Description','Time in UTC seconds');
	nc_attput(ncFile,'time_secs_FL','RawVarName','(''HH''*3600 + ''MM''*60 + ''SS'')');

	nc_attput(ncFile,'datetime_FL','Description','Date/time string (YYYYMMDDHHMMSS)');

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

	nc_attput(ncFile,'staticPres','units','mb');
	nc_attput(ncFile,'staticPres','Description','Static Pressure');
	nc_attput(ncFile,'staticPres','RawVarName',spSens);

	nc_attput(ncFile,'dynamicPres','units','mb');
	nc_attput(ncFile,'dynamicPres','Description','Dynamic Pressure');
	nc_attput(ncFile,'dynamicPres','RawVarName',dpSens);

	nc_attput(ncFile,'w_dpj','units','m/s');
	nc_attput(ncFile,'w_dpj','Description','DPJ calculated vertical wind velocity');
	nc_attput(ncFile,'w_dpj','RawVarName',wSens);

	nc_attput(ncFile,'u','units','m/s');
	nc_attput(ncFile,'u','Description','x-component of total wind');
	nc_attput(ncFile,'u','RawVarName',wind_xSens);

	nc_attput(ncFile,'v','units','m/s');
	nc_attput(ncFile,'v','Description','y-component of total wind');
	nc_attput(ncFile,'v','RawVarName',wind_ySens);

	nc_attput(ncFile,'u_relative','units','m/s');
	nc_attput(ncFile,'u_relative','Description','x-component of relative wind');
	nc_attput(ncFile,'u_relative','RawVarName',relWind_xSens);

	nc_attput(ncFile,'v_relative','units','m/s');
	nc_attput(ncFile,'v_relative','Description','y-component of relative wind');
	nc_attput(ncFile,'v_relative','RawVarName',relWind_ySens);

	nc_attput(ncFile,'windDir','units','degrees');
	nc_attput(ncFile,'windDir','Description','Wind direction (from north)');
	nc_attput(ncFile,'windDir','RawVarName',wdSens);

	nc_attput(ncFile,'windSpd','units','m/s');
	nc_attput(ncFile,'windSpd','Description','Wind speed');
	nc_attput(ncFile,'windSpd','RawVarName',wsSens);

	nc_attput(ncFile,'TAS','units','m/s');
	nc_attput(ncFile,'TAS','Description','True Air Speed');
	nc_attput(ncFile,'TAS','RawVarName','TAS.d');

	% Write actual data to each variable
	nc_varput(ncFile,'time_secs_FL',time_secs_FL);
	nc_varput(ncFile,'datetime_FL',fl_dt);
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
	nc_varput(ncFile,'staticPres',staticP);
	nc_varput(ncFile,'dynamicPres',dynamicP);
	nc_varput(ncFile,'w_dpj',w_dpj);
	nc_varput(ncFile,'u',u);
	nc_varput(ncFile,'v',v);
	nc_varput(ncFile,'u_relative',uRel);
	nc_varput(ncFile,'v_relative',vRel);
	nc_varput(ncFile,'windDir',windD);
	nc_varput(ncFile,'windSpd',windS);
	nc_varput(ncFile,'TAS',TAS);

	tmp1 = saveDir;
	tmp2 = flight;

	clearvars altStrct latStrct lonStrct timeStrct TAStrct TDStrct RHlStrct...
		RHiStrct RHhStrct RHrStrct TArStrct TDrStrct statPStrct dynmPStrct wdpjStrct...
		uStrct vStrct uRStrct vRStrct wdStrct wsStrct dtStrct tasStrct dewPtSens altSens rhSens spSens dpSens...
		tempSens latSens lonSens wSens wind_xSens wind_ySens relWind_xSens relWind_ySens wdSens...
		wsSens HH MM SS mDate mDatePre mDatePreStr dewExcd newSatValue	a0w a1w a2w a3w a4w a5w a6w...
		a0i a1i a2i a3i a4i a5i a6i fltLvlFile ix idt ncFile dataPath FLstr saveDir flight svp_ice...
		svp_lw vp_lw

	save([tmp1 tmp2 '_FltLvl_Processed.mat']);
	
	clearvars -except flights runAll iiz
end