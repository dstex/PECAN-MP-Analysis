% This script is intended to read in all relevant spiral-sorted microphysics/flight-level
% data from PECAN and save it out to a single netCDF file for each flight.
% This file can then be used in Python or any other language.
clearvars;


flight = '20150617';

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';


bmxComp = 0;

if bmxComp
    load([dataPath 'mp-data/' flight '/sDist/' flight '_Fit-CIP_10secAvg_12mm_wBAMEXmass.mat'],'-regexp','hybrid_igf|cipExt_|cip_igf_|negLmdaIx|fitSkipIx');
    outfile = [dataPath 'mp-data/' flight '/' flight '_CIPfit-spirals-10s1sAvg_wBAMEXmass.nc'];
else
    load([dataPath 'mp-data/' flight '/sDist/' flight '_Fit-CIP_10secAvg_12mm.mat'],'-regexp','hybrid_igf|cipExt_|cip_igf_|negLmdaIx|fitSkipIx');
    outfile = [dataPath 'mp-data/' flight '/' flight '_CIPfit-spirals-10s1sAvg.nc'];
end

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

%% Remove ML from select variables
sprlNames = fieldnames(cipTWC_gm3_hybrid_igf);
for ix = 1:length(sprlNames)
    cipLmda_10s.(sprlNames{ix}) = cip_igf_nml.(sprlNames{ix})(:,3);
	cipTWC_gm3_hybrid_igf_mlr.(sprlNames{ix}) = cipTWC_gm3_hybrid_igf.(sprlNames{ix});
    if bmxComp
        cipTWCBMX29jun_gm3_hybrid_igf_mlr.(sprlNames{ix}) = cipTWCBMX29jun_gm3_hybrid_igf.(sprlNames{ix});
        cipTWCBMX3jul_gm3_hybrid_igf_mlr.(sprlNames{ix}) = cipTWCBMX3jul_gm3_hybrid_igf.(sprlNames{ix});
        cipTWCBMX6jul_gm3_hybrid_igf_mlr.(sprlNames{ix}) = cipTWCBMX6jul_gm3_hybrid_igf.(sprlNames{ix});
    end
	cipDmm_mm_hybrid_igf_mlr.(sprlNames{ix}) = cipDmm_mm_hybrid_igf.(sprlNames{ix});
	efct_rad_um_avg_mlr.(sprlNames{ix}) = efct_rad_um_avg.(sprlNames{ix});
    
    
    tempC_60s.(sprlNames{ix}) = arrayfun(@(i) nanmean(tempC_avg.(sprlNames{ix})(i:i+6-1)),1:6:length(tempC_avg.(sprlNames{ix}))-6+1)';
    cipLmda_60s.(sprlNames{ix}) = arrayfun(@(i) nanmean(cipLmda_10s.(sprlNames{ix})(i:i+6-1)),1:6:length(cipLmda_10s.(sprlNames{ix}))-6+1)';
    
	if ~isnan(mlTopTemp(ix)) && ~isnan(mlBotTemp(ix))
		mlIXs = find(tempC_avg.(sprlNames{ix}) <= mlBotTemp(ix) & tempC_avg.(sprlNames{ix}) >= mlTopTemp(ix));
		cipTWC_gm3_hybrid_igf_mlr.(sprlNames{ix})(mlIXs) = NaN;
        if bmxComp
            cipTWCBMX29jun_gm3_hybrid_igf_mlr.(sprlNames{ix})(mlIXs) = NaN;
            cipTWCBMX3jul_gm3_hybrid_igf_mlr.(sprlNames{ix})(mlIXs) = NaN;
            cipTWCBMX6jul_gm3_hybrid_igf_mlr.(sprlNames{ix})(mlIXs) = NaN;
        end
		cipDmm_mm_hybrid_igf_mlr.(sprlNames{ix})(mlIXs) = NaN;
		efct_rad_um_avg_mlr.(sprlNames{ix})(mlIXs) = NaN;
	end
end


%% netCDF metadata and structure definition
ncRoot.Name = '/';
ncRoot.Format = 'netcdf4';
ncRoot.Dimensions(1).Name = 'obsBins';
ncRoot.Dimensions(1).Length = length(bin_min_mm);
ncRoot.Dimensions(2).Name = 'extndBins';
ncRoot.Dimensions(2).Length = length(cipExt_binMin_mm);
ncRoot.Dimensions(3).Name = 'habits';
ncRoot.Dimensions(3).Length = 10;
ncRoot.Dimensions(4).Name = 'spirals';
ncRoot.Dimensions(4).Length = length(startT);
ncRoot.Attributes(1).Name = 'flight';
ncRoot.Attributes(1).Value = flight;

zoneA(1).Name = 'Units';
zoneA(1).Value = 'T = Transition Zone; S = Enhanced Stratiform Region; A = Anvil Region; X = N/A';
zoneA(2).Name = 'Description';
zoneA(2).Value = 'Location of spiral relative to MCS structure';
ncRoot.Variables(1).Name = 'sprlZone';
ncRoot.Variables(1).Dimensions(1) = ncRoot.Dimensions(4);
ncRoot.Variables(1).Attributes = zoneA;
ncRoot.Variables(1).Datatype = 'char';

mTypeA(1).Name = 'Units';
mTypeA(1).Value = 'T=Trailing Stratiform,,L=Leading Stratiform,P=Parallel Stratiform,C=Cluster MCS,F=Post-Frontal';
mTypeA(2).Name = 'Description';
mTypeA(2).Value = 'General MCS system type at time of each spiral';
ncRoot.Variables(2).Name = 'mcsType';
ncRoot.Variables(2).Dimensions(1) = ncRoot.Dimensions(4);
ncRoot.Variables(2).Attributes = mTypeA;
ncRoot.Variables(2).Datatype = 'char';

ncwriteschema(outfile,ncRoot);

ncwrite(outfile,'/sprlZone',sprlZone);
ncwrite(outfile,'/mcsType',mcsType);


for ix = 1:length(sprlNames)
	eval(['nc' sprlNames{ix} '.Name = ''spiral_' num2str(ix) ''';']);
	eval(['nc' sprlNames{ix} '.Format = ''netcdf4'';']);
	
	
	eval(['nc' sprlNames{ix} '.Dimensions(1).Name = ''cipTsec_1s'';']);
	eval(['nc' sprlNames{ix} '.Dimensions(1).Length = ' num2str(length(time_secs_orig.(sprlNames{ix}))) ';']);
	eval(['nc' sprlNames{ix} '.Dimensions(2).Name = ''cipTsec_10s'';']);
	eval(['nc' sprlNames{ix} '.Dimensions(2).Length = ' num2str(length(time_secs_avg.(sprlNames{ix}))) ';']);
	eval(['nc' sprlNames{ix} '.Dimensions(3).Name = ''flTsec_1s'';']);
	eval(['nc' sprlNames{ix} '.Dimensions(3).Length = ' num2str(length(time_secsFL_orig.(sprlNames{ix}))) ';']);
	eval(['nc' sprlNames{ix} '.Dimensions(4).Name = ''flTemp_10s'';']);
	eval(['nc' sprlNames{ix} '.Dimensions(4).Length = ' num2str(length(time_secsFL_avg.(sprlNames{ix}))) ';']);
    eval(['nc' sprlNames{ix} '.Dimensions(5).Name = ''tempC_60s'';']);
	eval(['nc' sprlNames{ix} '.Dimensions(5).Length = ' num2str(length(tempC_60s.(sprlNames{ix}))) ';']);
    eval(['nc' sprlNames{ix} '.Dimensions(6).Name = ''cipLmda_60s'';']);
	eval(['nc' sprlNames{ix} '.Dimensions(6).Length = ' num2str(length(cipLmda_60s.(sprlNames{ix}))) ';']);
	
	%% CIP Time
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
	%% FL Time
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
	%% Nt
	ntA10(1).Name = 'Units';
	ntA10(1).Value = 'cm-3';
	ntA10(2).Name = 'Description';
	ntA10(2).Value = 'Total number concentration - CIP obs + IGF fit';
	eval(['nc' sprlNames{ix} '.Variables(5).Name = ''cipNt_hybrid_igf'';']);
	eval(['nc' sprlNames{ix} '.Variables(5).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(5).Attributes = ntA10;']);
	eval(['nc' sprlNames{ix} '.Variables(5).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(5).FillValue = NaN;']);
	%% TWC
	twcA10(1).Name = 'Units';
	twcA10(1).Value = 'g m-3';
	twcA10(2).Name = 'Description';
	twcA10(2).Value = 'Total water content - CIP obs + IGF fit';
	eval(['nc' sprlNames{ix} '.Variables(6).Name = ''cipTWC_hybrid_igf'';']);
	eval(['nc' sprlNames{ix} '.Variables(6).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(6).Attributes = twcA10;']);
	eval(['nc' sprlNames{ix} '.Variables(6).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(6).FillValue = NaN;']);
	%% Dmm
	dmmA10(1).Name = 'Units';
	dmmA10(1).Value = 'mm';
	dmmA10(2).Name = 'Description';
	dmmA10(2).Value = 'Median mass diameter - CIP obs + IGF fit';
	eval(['nc' sprlNames{ix} '.Variables(7).Name = ''cipDmm_hybrid_igf'';']);
	eval(['nc' sprlNames{ix} '.Variables(7).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(7).Attributes = dmmA10;']);
	eval(['nc' sprlNames{ix} '.Variables(7).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(7).FillValue = NaN;']);
	%% Temperature
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
	%% RH
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
	%% Area
	areaA1(1).Name = 'Units';
	areaA1(1).Value = 'mm2/cm4';
	areaA1(2).Name = 'Description';
	areaA1(2).Value = 'Projected area (extinction) - 1 sec average';
	eval(['nc' sprlNames{ix} '.Variables(12).Name = ''area_1s'';']);
	eval(['nc' sprlNames{ix} '.Variables(12).Dimensions(1) = ncRoot.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(12).Dimensions(2) = nc' sprlNames{ix} '.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(12).Attributes = areaA1;']);
	eval(['nc' sprlNames{ix} '.Variables(12).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(12).FillValue = NaN;']);
	
	areaA10(1).Name = 'Units';
	areaA10(1).Value = 'mm2/cm4';
	areaA10(2).Name = 'Description';
	areaA10(2).Value = 'Projected area (extinction) - 10 sec average';
	eval(['nc' sprlNames{ix} '.Variables(13).Name = ''area_10s'';']);
	eval(['nc' sprlNames{ix} '.Variables(13).Dimensions(1) = ncRoot.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(13).Dimensions(2) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(13).Attributes = areaA10;']);
	eval(['nc' sprlNames{ix} '.Variables(13).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(13).FillValue = NaN;']);
	
	areaCA1(1).Name = 'Units';
	areaCA1(1).Value = 'mm2/cm4';
	areaCA1(2).Name = 'Description';
	areaCA1(2).Value = 'Particle area calculated using A-D relations - 1 sec average';
	eval(['nc' sprlNames{ix} '.Variables(14).Name = ''area_calcd_1s'';']);
	eval(['nc' sprlNames{ix} '.Variables(14).Dimensions(1) = ncRoot.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(14).Dimensions(2) = nc' sprlNames{ix} '.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(14).Attributes = areaCA1;']);
	eval(['nc' sprlNames{ix} '.Variables(14).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(14).FillValue = NaN;']);
	
	areaCA10(1).Name = 'Units';
	areaCA10(1).Value = 'mm2/cm4';
	areaCA10(2).Name = 'Description';
	areaCA10(2).Value = 'Particle area calculated using A-D relations - 10 sec average';
	eval(['nc' sprlNames{ix} '.Variables(15).Name = ''area_calcd_10s'';']);
	eval(['nc' sprlNames{ix} '.Variables(15).Dimensions(1) = ncRoot.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(15).Dimensions(2) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(15).Attributes = areaCA10;']);
	eval(['nc' sprlNames{ix} '.Variables(15).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(15).FillValue = NaN;']);
	%% Area Ratio
	arA1(1).Name = 'Units';
	arA1(1).Value = '%';
	arA1(2).Name = 'Description';
	arA1(2).Value = 'Ratio of projected area to area of circle with diameter Dmax - 1 sec average';
	eval(['nc' sprlNames{ix} '.Variables(16).Name = ''areaRatio_1s'';']);
	eval(['nc' sprlNames{ix} '.Variables(16).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(16).Attributes = arA1;']);
	eval(['nc' sprlNames{ix} '.Variables(16).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(16).FillValue = NaN;']);
	
	arA10(1).Name = 'Units';
	arA10(1).Value = '%';
	arA10(2).Name = 'Description';
	arA10(2).Value = 'Ratio of projected area to area of circle with diameter Dmax - 10 sec average';
	eval(['nc' sprlNames{ix} '.Variables(17).Name = ''areaRatio_10s'';']);
	eval(['nc' sprlNames{ix} '.Variables(17).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(17).Attributes = arA10;']);
	eval(['nc' sprlNames{ix} '.Variables(17).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(17).FillValue = NaN;']);
	
	arMA1(1).Name = 'Units';
	arMA1(1).Value = '%';
	arMA1(2).Name = 'Description';
	arMA1(2).Value = 'Binned mean ratio of projected area to area of circle with diameter Dmax - 1 sec average';
	eval(['nc' sprlNames{ix} '.Variables(20).Name = ''meanAreaRatio_1s'';']);
	eval(['nc' sprlNames{ix} '.Variables(20).Dimensions(1) = ncRoot.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(20).Dimensions(2) = nc' sprlNames{ix} '.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(20).Attributes = arMA1;']);
	eval(['nc' sprlNames{ix} '.Variables(20).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(20).FillValue = NaN;']);
	
	arMA10(1).Name = 'Units';
	arMA10(1).Value = '%';
	arMA10(2).Name = 'Description';
	arMA10(2).Value = 'Binned mean ratio of projected area to area of circle with diameter Dmax - 10 sec average';
	eval(['nc' sprlNames{ix} '.Variables(21).Name = ''meanAreaRatio_10s'';']);
	eval(['nc' sprlNames{ix} '.Variables(21).Dimensions(1) = ncRoot.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(21).Dimensions(2) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(21).Attributes = arMA10;']);
	eval(['nc' sprlNames{ix} '.Variables(21).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(21).FillValue = NaN;']);
	%% Effective Radius
	erA1(1).Name = 'Units';
	erA1(1).Value = 'um';
	erA1(2).Name = 'Description';
	erA1(2).Value = 'Effective radius - 1 sec average';
	eval(['nc' sprlNames{ix} '.Variables(18).Name = ''efctvRadius_1s'';']);
	eval(['nc' sprlNames{ix} '.Variables(18).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(18).Attributes = erA1;']);
	eval(['nc' sprlNames{ix} '.Variables(18).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(18).FillValue = NaN;']);
	
	erA10(1).Name = 'Units';
	erA10(1).Value = 'um';
	erA10(2).Name = 'Description';
	erA10(2).Value = 'Effective radius - 10 sec average';
	eval(['nc' sprlNames{ix} '.Variables(19).Name = ''efctvRadius_10s'';']);
	eval(['nc' sprlNames{ix} '.Variables(19).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(19).Attributes = erA10;']);
	eval(['nc' sprlNames{ix} '.Variables(19).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(19).FillValue = NaN;']);
	%% Perimeter
	prmMA1(1).Name = 'Units';
	prmMA1(1).Value = 'um';
	prmMA1(2).Name = 'Description';
	prmMA1(2).Value = 'Binned mean perimeter - 1 sec average';
	eval(['nc' sprlNames{ix} '.Variables(22).Name = ''meanPerim_1s'';']);
	eval(['nc' sprlNames{ix} '.Variables(22).Dimensions(1) = ncRoot.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(22).Dimensions(2) = nc' sprlNames{ix} '.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(22).Attributes = prmMA1;']);
	eval(['nc' sprlNames{ix} '.Variables(22).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(22).FillValue = NaN;']);
	
	prmMA10(1).Name = 'Units';
	prmMA10(1).Value = 'um';
	prmMA10(2).Name = 'Description';
	prmMA10(2).Value = 'Binned mean perimeter - 10 sec average';
	eval(['nc' sprlNames{ix} '.Variables(23).Name = ''meanPerim_10s'';']);
	eval(['nc' sprlNames{ix} '.Variables(23).Dimensions(1) = ncRoot.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(23).Dimensions(2) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(23).Attributes = prmMA10;']);
	eval(['nc' sprlNames{ix} '.Variables(23).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(23).FillValue = NaN;']);
	%% Reject Ratio
	rrA1(1).Name = 'Units';
	rrA1(1).Value = '%';
	rrA1(2).Name = 'Description';
	rrA1(2).Value = 'Reject ratio - 1 sec average';
	eval(['nc' sprlNames{ix} '.Variables(24).Name = ''rjctRatio_1s'';']);
	eval(['nc' sprlNames{ix} '.Variables(24).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(1);']);
	eval(['nc' sprlNames{ix} '.Variables(24).Attributes = rrA1;']);
	eval(['nc' sprlNames{ix} '.Variables(24).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(24).FillValue = NaN;']);
	
	rrA10(1).Name = 'Units';
	rrA10(1).Value = '%';
	rrA10(2).Name = 'Description';
	rrA10(2).Value = 'Reject ratio - 10 sec average';
	eval(['nc' sprlNames{ix} '.Variables(25).Name = ''rjctRatio_10s'';']);
	eval(['nc' sprlNames{ix} '.Variables(25).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(25).Attributes = rrA10;']);
	eval(['nc' sprlNames{ix} '.Variables(25).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(25).FillValue = NaN;']);
	
	%% ML Removed TWC
	twcA10mlr(1).Name = 'Units';
	twcA10mlr(1).Value = 'g m-3';
	twcA10mlr(2).Name = 'Description';
	twcA10mlr(2).Value = 'Total water content - CIP obs + IGF fit - Melt Layer Data Removed';
	eval(['nc' sprlNames{ix} '.Variables(26).Name = ''cipTWC_hybrid_igf_mlr'';']);
	eval(['nc' sprlNames{ix} '.Variables(26).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(26).Attributes = twcA10mlr;']);
	eval(['nc' sprlNames{ix} '.Variables(26).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(26).FillValue = NaN;']);
    
    if bmxComp
        %% ML Removed BAMEX 29 June spiral 2 TWC
        twcA10B1mlr(1).Name = 'Units';
        twcA10B1mlr(1).Value = 'g m-3';
        twcA10B1mlr(2).Name = 'Description';
        twcA10B1mlr(2).Value = 'Total water content using BAMEX 29 Jun sprl2 m-D relat. - CIP obs + IGF fit - Melt Layer Data Removed';
        eval(['nc' sprlNames{ix} '.Variables(32).Name = ''cipTWCBMX29jun_hybrid_igf_mlr'';']);
        eval(['nc' sprlNames{ix} '.Variables(32).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
        eval(['nc' sprlNames{ix} '.Variables(32).Attributes = twcA10B1mlr;']);
        eval(['nc' sprlNames{ix} '.Variables(32).Datatype = ''double'';']);
        eval(['nc' sprlNames{ix} '.Variables(32).FillValue = NaN;']);
        %% ML Removed BAMEX 3 July spiral 1 TWC
        twcA10B2mlr(1).Name = 'Units';
        twcA10B2mlr(1).Value = 'g m-3';
        twcA10B2mlr(2).Name = 'Description';
        twcA10B2mlr(2).Value = 'Total water content using BAMEX 3 Jul sprl1 m-D relat. - CIP obs + IGF fit - Melt Layer Data Removed';
        eval(['nc' sprlNames{ix} '.Variables(33).Name = ''cipTWCBMX3jul_hybrid_igf_mlr'';']);
        eval(['nc' sprlNames{ix} '.Variables(33).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
        eval(['nc' sprlNames{ix} '.Variables(33).Attributes = twcA10B2mlr;']);
        eval(['nc' sprlNames{ix} '.Variables(33).Datatype = ''double'';']);
        eval(['nc' sprlNames{ix} '.Variables(33).FillValue = NaN;']);
        %% ML Removed BAMEX 6 July spiral 2 TWC
        twcA10B3mlr(1).Name = 'Units';
        twcA10B3mlr(1).Value = 'g m-3';
        twcA10B3mlr(2).Name = 'Description';
        twcA10B3mlr(2).Value = 'Total water content using BAMEX 6 Jul sprl2 m-D relat. - CIP obs + IGF fit - Melt Layer Data Removed';
        eval(['nc' sprlNames{ix} '.Variables(34).Name = ''cipTWCBMX6jul_hybrid_igf_mlr'';']);
        eval(['nc' sprlNames{ix} '.Variables(34).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
        eval(['nc' sprlNames{ix} '.Variables(34).Attributes = twcA10B3mlr;']);
        eval(['nc' sprlNames{ix} '.Variables(34).Datatype = ''double'';']);
        eval(['nc' sprlNames{ix} '.Variables(34).FillValue = NaN;']);
    end
	%% ML Removed Dmm
	dmmA10mlr(1).Name = 'Units';
	dmmA10mlr(1).Value = 'mm';
	dmmA10mlr(2).Name = 'Description';
	dmmA10mlr(2).Value = 'Median mass diameter - CIP obs + IGF fit - Melt Layer Data Removed';
	eval(['nc' sprlNames{ix} '.Variables(27).Name = ''cipDmm_hybrid_igf_mlr'';']);
	eval(['nc' sprlNames{ix} '.Variables(27).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(27).Attributes = dmmA10mlr;']);
	eval(['nc' sprlNames{ix} '.Variables(27).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(27).FillValue = NaN;']);
	%% ML Removed Effective Radius
	erA10mlr(1).Name = 'Units';
	erA10mlr(1).Value = 'um';
	erA10mlr(2).Name = 'Description';
	erA10mlr(2).Value = 'Effective radius - 10 sec average - Melt Layer Data Removed';
	eval(['nc' sprlNames{ix} '.Variables(28).Name = ''efctvRadius_10s_mlr'';']);
	eval(['nc' sprlNames{ix} '.Variables(28).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(28).Attributes = erA10mlr;']);
	eval(['nc' sprlNames{ix} '.Variables(28).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(28).FillValue = NaN;']);
    
    
    %% Temperature 60 sec Avg
	tempA60(1).Name = 'Units';
	tempA60(1).Value = 'deg C';
	tempA60(2).Name = 'Description';
	tempA60(2).Value = 'Flight-level temperature (w/Zipser wetting corr) - 60 sec average';
	eval(['nc' sprlNames{ix} '.Variables(29).Name = ''tempC_60s'';']);
	eval(['nc' sprlNames{ix} '.Variables(29).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(5);']);
	eval(['nc' sprlNames{ix} '.Variables(29).Attributes = tempA60;']);
	eval(['nc' sprlNames{ix} '.Variables(29).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(29).FillValue = NaN;']);
    
    %% Lambda 60 sec Avg
	lmdaA60(1).Name = 'Units';
	lmdaA60(1).Value = 'cm-1';
	lmdaA60(2).Name = 'Description';
	lmdaA60(2).Value = 'IGF Lambda - 60 sec average';
	eval(['nc' sprlNames{ix} '.Variables(30).Name = ''cipLmda_60s'';']);
	eval(['nc' sprlNames{ix} '.Variables(30).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(6);']);
	eval(['nc' sprlNames{ix} '.Variables(30).Attributes = lmdaA60;']);
	eval(['nc' sprlNames{ix} '.Variables(30).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(30).FillValue = NaN;']);
    
    %% Lambda 10 sec Avg
	lmdaA10(1).Name = 'Units';
	lmdaA10(1).Value = 'cm-1';
	lmdaA10(2).Name = 'Description';
	lmdaA10(2).Value = 'IGF Lambda - 10 sec average';
	eval(['nc' sprlNames{ix} '.Variables(31).Name = ''cipLmda_10s'';']);
	eval(['nc' sprlNames{ix} '.Variables(31).Dimensions(1) = nc' sprlNames{ix} '.Dimensions(2);']);
	eval(['nc' sprlNames{ix} '.Variables(31).Attributes = lmdaA10;']);
	eval(['nc' sprlNames{ix} '.Variables(31).Datatype = ''double'';']);
	eval(['nc' sprlNames{ix} '.Variables(31).FillValue = NaN;']);
	
	%% Apply metadata and write out data to file
	eval(['ncwriteschema(outfile, nc' sprlNames{ix} ');']);
	
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/cipTsec_1s'',time_secs_orig.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/cipTsec_10s'',time_secs_avg.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/flTsec_1s'',time_secsFL_orig.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/flTsec_10s'',time_secsFL_avg.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/cipNt_hybrid_igf'',cipNt_cm3_hybrid_igf.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/cipTWC_hybrid_igf'',cipTWC_gm3_hybrid_igf.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/cipDmm_hybrid_igf'',cipDmm_mm_hybrid_igf.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/tempC_1s'',tempC_orig.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/tempC_10s'',tempC_avg.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/rh_1s'',RH_orig.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/rh_10s'',RH_avg.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/area_1s'',permute(total_area_mm2cm4_orig.(sprlNames{ix}),[2 1]));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/area_10s'',permute(total_area_mm2cm4_avg.(sprlNames{ix}),[2 1]));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/area_calcd_1s'',permute(area_calcd_mm2cm4_orig.(sprlNames{ix}),[2 1]));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/area_calcd_10s'',permute(area_calcd_mm2cm4_avg.(sprlNames{ix}),[2 1]));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/areaRatio_1s'',area_ratio_orig.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/areaRatio_10s'',area_ratio_avg.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/efctvRadius_1s'',efct_rad_um_orig.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/efctvRadius_10s'',efct_rad_um_avg.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/meanAreaRatio_1s'',permute(mean_areaRatio_orig.(sprlNames{ix}),[2 1]));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/meanAreaRatio_10s'',permute(mean_areaRatio_avg.(sprlNames{ix}),[2 1]));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/meanPerim_1s'',permute(mean_perim_um_orig.(sprlNames{ix}),[2 1]));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/meanPerim_10s'',permute(mean_perim_um_avg.(sprlNames{ix}),[2 1]));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/rjctRatio_1s'',reject_ratio_orig.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/rjctRatio_10s'',reject_ratio_avg.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/cipTWC_hybrid_igf_mlr'',cipTWC_gm3_hybrid_igf_mlr.(sprlNames{ix}));']);
    if bmxComp
        eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/cipTWCBMX29jun_hybrid_igf_mlr'',cipTWCBMX29jun_gm3_hybrid_igf_mlr.(sprlNames{ix}));']);
        eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/cipTWCBMX3jul_hybrid_igf_mlr'',cipTWCBMX3jul_gm3_hybrid_igf_mlr.(sprlNames{ix}));']);
        eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/cipTWCBMX6jul_hybrid_igf_mlr'',cipTWCBMX6jul_gm3_hybrid_igf_mlr.(sprlNames{ix}));']);
    end
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/cipDmm_hybrid_igf_mlr'',cipDmm_mm_hybrid_igf_mlr.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/efctvRadius_10s_mlr'',efct_rad_um_avg_mlr.(sprlNames{ix}));']);
	eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/cipLmda_10s'',cipLmda_10s.(sprlNames{ix}));']);
    eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/cipLmda_60s'',cipLmda_60s.(sprlNames{ix}));']);
    eval(['ncwrite(outfile, ''/spiral_' num2str(ix) '/tempC_60s'',tempC_60s.(sprlNames{ix}));']);
end