% This script is intended to read in all relevant spiral-sorted microphysics/flight-level
% data from PECAN and save it out to a single netCDF file for each flight.
% This file can then be used in Python or any other language.
clearvars;


flight = '20150617';

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

outfile = [dataPath 'mp-data/' flight '/' flight '_CIPobs-spirals-10s1sAvg.nc'];

if (exist(outfile, 'file') == 2)
	delete(outfile);
end

load([dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.CIP.10secAvg.mat'],'-regexp','_orig|_avg|bin_');


startT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'startT');
endT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'endT');
mlBotTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlBotTime');
mlTopTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlTopTime');
mlBotTemp = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlBotTemp');
mlTopTemp = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlTopTemp');
sprlZone = nc_varget([dataPath '/' flight '_PECANparams.nc'],'sprlZone');
mcsType = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mcsType');



ncRoot.Name = '/';
ncRoot.Format = 'netcdf4';
ncRoot.Dimensions(1).Name = 'obsBins';
ncRoot.Dimensions(1).Length = length(bin_min);
ncRoot.Dimensions(2).Name = 'spirals';
ncRoot.Dimensions(2).Length = length(startT);
ncRoot.Attributes(1).Name = 'flight';
ncRoot.Attributes(1).Value = flight;

zoneA(1).Name = 'Units';
zoneA(1).Value = 'T = Transition Zone; S = Stratiform Region; A = Rear Anvil; U = Unclassified';
zoneA(2).Name = 'Description';
zoneA(2).Value = 'Location of spiral relative to MCS structure';
ncRoot.Variables(1).Name = 'sprlZone';
ncRoot.Variables(1).Dimensions(1) = ncRoot.Dimensions(2);
ncRoot.Variables(1).Attributes = zoneA;
ncRoot.Variables(1).Datatype = 'char';

mTypeA(1).Name = 'Units';
mTypeA(1).Value = 'F=Trailing Stratiform (Formative),M=Trailing Stratiform (Mature),L=Leading Stratiform,P=Parallel Stratiform,C=Cluster MCS';
mTypeA(2).Name = 'Description';
mTypeA(2).Value = 'General MCS system type at time of each spiral';
ncRoot.Variables(2).Name = 'mcsType';
ncRoot.Variables(2).Dimensions(1) = ncRoot.Dimensions(2);
ncRoot.Variables(2).Attributes = mTypeA;
ncRoot.Variables(2).Datatype = 'char';

ncwriteschema(outfile,ncRoot);

ncwrite(outfile,'/sprlZone',sprlZone);
ncwrite(outfile,'/mcsType',mcsType);

sprlNames = fieldnames(time_secs_avg);



for ix = 1:length(sprlNames)
	eval(['nc' sprlNames{ix} '.Name = ''spiral_' num2str(ix) ''';']);
	eval(['nc' sprlNames{ix} '.Format = ''netcdf4'';']);
	
	
	eval(['nc' sprlNames{ix} '.Dimensions(1).Name = ''cipTsec_1s'';']);
	eval(['nc' sprlNames{ix} '.Dimensions(1).Length = ' num2str(length(time_secs_orig.(sprlNames{ix}))) ';']);
	eval(['nc' sprlNames{ix} '.Dimensions(2).Name = ''cipTsec_10s'';']);
	eval(['nc' sprlNames{ix} '.Dimensions(2).Length = ' num2str(length(time_secs_avg.(sprlNames{ix}))) ';']);
	eval(['nc' sprlNames{ix} '.Dimensions(3).Name = ''flTsec_1s'';']);
	eval(['nc' sprlNames{ix} '.Dimensions(3).Length = ' num2str(length(time_secsFL_orig.(sprlNames{ix}))) ';']);
	eval(['nc' sprlNames{ix} '.Dimensions(4).Name = ''flTsec_10s'';']);
	eval(['nc' sprlNames{ix} '.Dimensions(4).Length = ' num2str(length(time_secsFL_avg.(sprlNames{ix}))) ';']);
	
	
	timeA1(1).Name = 'Units';
	timeA1(1).Value = 'Seconds since midnight UTC';
	timeA1(2).Name = 'Description';
	timeA1(2).Value = 'Time in seconds since UTC midnight - 1 sec average';
	eval(['nc' sprlNames{ix} '.Variables(1).Name = ''cipTsec_1s'';']);
	eval(['nc' sprlNames{ix} '.Variables(1).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(1).Attributes = timeA1;']);
	eval(['nc' sprlNames{ix} '.Variables(1).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(1).FillValue = NaN;']);
	
	timeA10(1).Name = 'Units';
	timeA10(1).Value = 'Seconds since midnight UTC';
	timeA10(2).Name = 'Description';
	timeA10(2).Value = 'Time in seconds since UTC midnight - 10 sec average';
	eval(['nc' sprlNames{ix} '.Variables(2).Name = ''cipTsec_10s'';']);
	eval(['nc' sprlNames{ix} '.Variables(2).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(2).Attributes = timeA10;']);
	eval(['nc' sprlNames{ix} '.Variables(2).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(2).FillValue = NaN;']);
	
	timeFlA1(1).Name = 'Units';
	timeFlA1(1).Value = 'Seconds since midnight UTC';
	timeFlA1(2).Name = 'Description';
	timeFlA1(2).Value = 'Flight-level time in seconds since UTC midnight - 1 sec average';
	eval(['nc' sprlNames{ix} '.Variables(3).Name = ''flTsec_1s'';']);
	eval(['nc' sprlNames{ix} '.Variables(3).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(3);']);
	eval(['nc' sprlNames{ix} '.Variables(3).Attributes = timeFlA1;']);
	eval(['nc' sprlNames{ix} '.Variables(3).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(3).FillValue = NaN;']);
	
	timeFlA10(1).Name = 'Units';
	timeFlA10(1).Value = 'Seconds since midnight UTC';
	timeFlA10(2).Name = 'Description';
	timeFlA10(2).Value = 'Flight-level time in seconds since UTC midnight - 10 sec average';
	eval(['nc' sprlNames{ix} '.Variables(4).Name = ''flTsec_10s'';']);
	eval(['nc' sprlNames{ix} '.Variables(4).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(4);']);
	eval(['nc' sprlNames{ix} '.Variables(4).Attributes = timeFlA1;']);
	eval(['nc' sprlNames{ix} '.Variables(4).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(4).FillValue = NaN;']);
	
	ntA10(1).Name = 'Units';
	ntA10(1).Value = 'cm-3';
	ntA10(2).Name = 'Description';
	ntA10(2).Value = 'Total number concentration';
	eval(['nc' sprlNames{ix} '.Variables(5).Name = ''cipNt'';']);
	eval(['nc' sprlNames{ix} '.Variables(5).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(5).Attributes = ntA10;']);
	eval(['nc' sprlNames{ix} '.Variables(5).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(5).FillValue = NaN;']);
	
	twcA10(1).Name = 'Units';
	twcA10(1).Value = 'g m-3';
	twcA10(2).Name = 'Description';
	twcA10(2).Value = 'Total water content';
	eval(['nc' sprlNames{ix} '.Variables(6).Name = ''cipTWC'';']);
	eval(['nc' sprlNames{ix} '.Variables(6).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(6).Attributes = twcA10;']);
	eval(['nc' sprlNames{ix} '.Variables(6).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(6).FillValue = NaN;']);
	
	dmmA10(1).Name = 'Units';
	dmmA10(1).Value = 'mm';
	dmmA10(2).Name = 'Description';
	dmmA10(2).Value = 'Median mass diameter';
	eval(['nc' sprlNames{ix} '.Variables(7).Name = ''cipDmm'';']);
	eval(['nc' sprlNames{ix} '.Variables(7).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(7).Attributes = dmmA10;']);
	eval(['nc' sprlNames{ix} '.Variables(7).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(7).FillValue = NaN;']);
	
	tempA1(1).Name = 'Units';
	tempA1(1).Value = 'deg C';
	tempA1(2).Name = 'Description';
	tempA1(2).Value = 'Flight-level temperature (w/Zipser wetting corr) - 1 sec average';
	eval(['nc' sprlNames{ix} '.Variables(8).Name = ''tempC_1s'';']);
	eval(['nc' sprlNames{ix} '.Variables(8).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(3);']);
	eval(['nc' sprlNames{ix} '.Variables(8).Attributes = tempA1;']);
	eval(['nc' sprlNames{ix} '.Variables(8).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(8).FillValue = NaN;']);
	
	tempA10(1).Name = 'Units';
	tempA10(1).Value = 'deg C';
	tempA10(2).Name = 'Description';
	tempA10(2).Value = 'Flight-level temperature (w/Zipser wetting corr) - 10 sec average';
	eval(['nc' sprlNames{ix} '.Variables(9).Name = ''tempC_10s'';']);
	eval(['nc' sprlNames{ix} '.Variables(9).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(4);']);
	eval(['nc' sprlNames{ix} '.Variables(9).Attributes = tempA10;']);
	eval(['nc' sprlNames{ix} '.Variables(9).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(9).FillValue = NaN;']);
	
	rhA1(1).Name = 'Units';
	rhA1(1).Value = '%';
	rhA1(2).Name = 'Description';
	rhA1(2).Value = 'Flight-level relative humidity (w.r.t. water above 0C; w.r.t. ice at/below 0C) - 1 sec average';
	eval(['nc' sprlNames{ix} '.Variables(10).Name = ''rh_1s'';']);
	eval(['nc' sprlNames{ix} '.Variables(10).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(3);']);
	eval(['nc' sprlNames{ix} '.Variables(10).Attributes = rhA1;']);
	eval(['nc' sprlNames{ix} '.Variables(10).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(10).FillValue = NaN;']);
	
	rhA10(1).Name = 'Units';
	rhA10(1).Value = '%';
	rhA10(2).Name = 'Description';
	rhA10(2).Value = 'Flight-level relative humidity (w.r.t. water above 0C; w.r.t. ice at/below 0C) - 10 sec average';
	eval(['nc' sprlNames{ix} '.Variables(11).Name = ''rh_10s'';']);
	eval(['nc' sprlNames{ix} '.Variables(11).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(4);']);
	eval(['nc' sprlNames{ix} '.Variables(11).Attributes = rhA10;']);
	eval(['nc' sprlNames{ix} '.Variables(11).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(11).FillValue = NaN;']);

	
	arA1(1).Name = 'Units';
	arA1(1).Value = '%';
	arA1(2).Name = 'Description';
	arA1(2).Value = 'Ratio of projected area to area of circle with diameter Dmax - 1 sec average';
	eval(['nc' sprlNames{ix} '.Variables(12).Name = ''areaRatio_1s'';']);
	eval(['nc' sprlNames{ix} '.Variables(12).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(12).Attributes = arA1;']);
	eval(['nc' sprlNames{ix} '.Variables(12).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(12).FillValue = NaN;']);
	
	arA10(1).Name = 'Units';
	arA10(1).Value = '%';
	arA10(2).Name = 'Description';
	arA10(2).Value = 'Ratio of projected area to area of circle with diameter Dmax - 10 sec average';
	eval(['nc' sprlNames{ix} '.Variables(13).Name = ''areaRatio_10s'';']);
	eval(['nc' sprlNames{ix} '.Variables(13).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(13).Attributes = arA10;']);
	eval(['nc' sprlNames{ix} '.Variables(13).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(13).FillValue = NaN;']);
	
	erA1(1).Name = 'Units';
	erA1(1).Value = 'um';
	erA1(2).Name = 'Description';
	erA1(2).Value = 'Effective radius - 1 sec average';
	eval(['nc' sprlNames{ix} '.Variables(14).Name = ''efctvRadius_1s'';']);
	eval(['nc' sprlNames{ix} '.Variables(14).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(14).Attributes = erA1;']);
	eval(['nc' sprlNames{ix} '.Variables(14).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(14).FillValue = NaN;']);
	
	erA10(1).Name = 'Units';
	erA10(1).Value = 'um';
	erA10(2).Name = 'Description';
	erA10(2).Value = 'Effective radius - 10 sec average';
	eval(['nc' sprlNames{ix} '.Variables(15).Name = ''efctvRadius_10s'';']);
	eval(['nc' sprlNames{ix} '.Variables(15).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(15).Attributes = erA10;']);
	eval(['nc' sprlNames{ix} '.Variables(15).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(15).FillValue = NaN;']);
	
	arMA1(1).Name = 'Units';
	arMA1(1).Value = '%';
	arMA1(2).Name = 'Description';
	arMA1(2).Value = 'Binned mean ratio of projected area to area of circle with diameter Dmax - 1 sec average';
	eval(['nc' sprlNames{ix} '.Variables(16).Name = ''meanAreaRatio_1s'';']);
	eval(['nc' sprlNames{ix} '.Variables(16).Dimensions(1) = ncRoot.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(16).Dimensions(2) = nc' sprlNames{ix} '.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(16).Attributes = arMA1;']);
	eval(['nc' sprlNames{ix} '.Variables(16).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(16).FillValue = NaN;']);
	
	arMA10(1).Name = 'Units';
	arMA10(1).Value = '%';
	arMA10(2).Name = 'Description';
	arMA10(2).Value = 'Binned mean ratio of projected area to area of circle with diameter Dmax - 10 sec average';
	eval(['nc' sprlNames{ix} '.Variables(17).Name = ''meanAreaRatio_10s'';']);
	eval(['nc' sprlNames{ix} '.Variables(17).Dimensions(1) = ncRoot.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(17).Dimensions(2) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(17).Attributes = arMA10;']);
	eval(['nc' sprlNames{ix} '.Variables(17).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(17).FillValue = NaN;']);
	
	
	eval(['ncwriteschema(outfile, nc' sprlNames{ix} ');']);
	
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/cipTsec_1s'',time_secs_orig.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/cipTsec_10s'',time_secs_avg.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/flTsec_1s'',time_secsFL_orig.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/flTsec_10s'',time_secsFL_avg.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/cipNt'',n_avg.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/cipTWC'',twc_avg.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/cipDmm'',Dmm_twc_avg.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/tempC_1s'',tempC_orig.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/tempC_10s'',tempC_avg.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/rh_1s'',RH_orig.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/rh_10s'',RH_avg.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/areaRatio_1s'',area_ratio_orig.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/areaRatio_10s'',area_ratio_avg.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/efctvRadius_1s'',efct_rad_orig.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/efctvRadius_10s'',efct_rad_avg.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/meanAreaRatio_1s'',permute(mean_areaRatio_orig.(sprlNames{ix}),[2 1]));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/meanAreaRatio_10s'',permute(mean_areaRatio_avg.(sprlNames{ix}),[2 1]));']);
	
end