% Compute 5- or 10-second averaged size distributions and flight-level
% data for each spiral of a particular flight

clearvars

flight = '20150709';

probe = 'CIP';

avgTime = 10; % sec

calcRej = 0;

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

startT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'startT');
endT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'endT');
mlBotTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlBotTime');

% outFileAppend = '_Dgt150um';
outFileAppend = '';

%% Define filenames
if calcRej
	sDistFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.' probe 'wRejSD.cdf'];
else
	sDistFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.' probe '.cdf'];
end
flFile = [dataPath 'FlightLevelData/Processed/' flight '_FltLvl_Processed.mat'];


%% Import all size distribution variables

timehhmmss = nc_varget(sDistFile,'time');

bin_min = nc_varget(sDistFile,'bin_min');
bin_max = nc_varget(sDistFile,'bin_max');
bin_mid = nc_varget(sDistFile,'bin_mid');
bin_size = nc_varget(sDistFile,'bin_dD');

conc_minR_all = (nc_varget(sDistFile,'conc_minR'))';
area_all = permute(nc_varget(sDistFile,'area'),[3 2 1]);
area_calcd_all = (nc_varget(sDistFile,'Calcd_area'))';
conc_AreaR_all = (nc_varget(sDistFile,'conc_AreaR'))';
n_all = nc_varget(sDistFile,'n');
total_area_all = (nc_varget(sDistFile,'total_area'))';
mass_ice_all = (nc_varget(sDistFile,'mass_ice'))';
mass_lw_all = (nc_varget(sDistFile,'mass_lw'))';
massBL_all = (nc_varget(sDistFile,'massBL'))';
sd_habit_all = permute(nc_varget(sDistFile,'habitsd'),[3 2 1]);
sdMass_habit_all = permute(nc_varget(sDistFile,'habitmsd'),[3 2 1]);
efct_rad_all = nc_varget(sDistFile,'re');
area_ratio_all = nc_varget(sDistFile,'ar');
reject_ratio_all = nc_varget(sDistFile,'Reject_ratio');
termVeloc_ice_all = (nc_varget(sDistFile,'vt_ice'))';
termVeloc_lw_all = (nc_varget(sDistFile,'vt_lw'))';
prec_rate_ice_all = (nc_varget(sDistFile,'Prec_rate_ice'))';
prec_rate_lw_all = (nc_varget(sDistFile,'Prec_rate_lw'))';
count_all = (nc_varget(sDistFile,'count'))';
sampleVol_all = (nc_varget(sDistFile,'sample_vol'))';

mean_aspectRatio_rect_all = (nc_varget(sDistFile,'mean_aspect_ratio_rectangle'))';
mean_aspectRatio_elps_all = (nc_varget(sDistFile,'mean_aspect_ratio_ellipse'))';
mean_areaRatio_all = (nc_varget(sDistFile,'mean_area_ratio'))';
mean_perim_all = (nc_varget(sDistFile,'mean_perimeter'))';

if calcRej
	rej_conc_minR_all = (nc_varget(sDistFile,'REJ_conc_minR'))';
	rej_area_all = permute(nc_varget(sDistFile,'REJ_area'),[3 2 1]);
	rej_area_calcd_all = (nc_varget(sDistFile,'REJ_Calcd_area'))';
	rej_conc_AreaR_all = (nc_varget(sDistFile,'REJ_conc_AreaR'))';
	rej_n_all = nc_varget(sDistFile,'REJ_n');
	rej_total_area_all = (nc_varget(sDistFile,'REJ_total_area'))';
	rej_mass_all = (nc_varget(sDistFile,'REJ_mass'))';
	rej_massBL_all = (nc_varget(sDistFile,'REJ_massBL'))';
	rej_sd_habit_all = permute(nc_varget(sDistFile,'REJ_habitsd'),[3 2 1]);
	rej_sdMass_habit_all = permute(nc_varget(sDistFile,'REJ_habitmsd'),[3 2 1]);
	rej_efct_rad_all = nc_varget(sDistFile,'REJ_re');
	rej_area_ratio_all = nc_varget(sDistFile,'REJ_ar');
	rej_termVeloc_all = (nc_varget(sDistFile,'REJ_vt'))';
	rej_prec_rate_all = (nc_varget(sDistFile,'REJ_Prec_rate'))';
	rej_count_all = (nc_varget(sDistFile,'REJ_count'))';
	
	rej_mean_aspectRatio_rect_all = (nc_varget(sDistFile,'REJ_mean_aspect_ratio_rectangle'))';
	rej_mean_aspectRatio_elps_all = (nc_varget(sDistFile,'REJ_mean_aspect_ratio_ellipse'))';
	rej_mean_areaRatio_all = (nc_varget(sDistFile,'REJ_mean_area_ratio'))';
	rej_mean_perim_all = (nc_varget(sDistFile,'REJ_mean_perimeter'))';
end

time_secs_all = hhmmss2insec(timehhmmss);

num_bins = length(bin_mid);


%% Import relevant flight-level data
importVars = {'time_secs_FL','TA','RH_hybrid','Alt'};

tempLoad = load(flFile,importVars{:});

time_secsFL_all = tempLoad.(importVars{1});
tempC_all = tempLoad.(importVars{2});
RH_all = tempLoad.(importVars{3});
alt_all = tempLoad.(importVars{4});


%% Compute any variables not available in the sDist output


for iii=1:length(time_secs_all)
	% Ice water content
	iwc_all_a(iii,:)=mass_ice_all(iii,:).*(bin_size')/10; % Multiply by bin size in cm [now g/cm3]
	
	% Liquid water content
	lwc_all_a(iii,:)=mass_lw_all(iii,:).*(bin_size')/10;
	
end

iwc_all = nansum(iwc_all_a,2)*1e6; % Sum all values in all bins and convert to g/m3

% Liquid water content
lwc_all = nansum(lwc_all_a,2)*1e6;

%% Determine periods over which to average

for ix=1:length(startT)
	sprlLocs = find(time_secs_all >= startT(ix) & time_secs_all <= endT(ix));
	sprlLocs_FL = find(time_secsFL_all >= startT(ix) & time_secsFL_all <= endT(ix));
	
	if (length(sprlLocs) ~= length(sprlLocs_FL))
		fprintf('Warning: Number of flight-level spiral locations does not equal number of sDist locations! Skipping (Spiral #%d)\n',ix);
		continue
	end

	spiralStr = ['sprl' num2str(ix)];
	
	% Determine the length of the averaged time dimension, and account for the remainder if necessary 
	if (rem(max(size(sprlLocs)),avgTime) ~= 0)
		neatLength = floor(max(size(sprlLocs))/avgTime)*avgTime; %Round length down to the nearest multiple of avgTime
		maxLength = floor(max(size(sprlLocs))/avgTime)+1; % Final length of time dimension for each variable (last value will be average of remainder)
		leftover = rem(max(size(sprlLocs)),avgTime);
		remain = 1;
	else
		neatLength = floor(max(size(sprlLocs))/avgTime)*avgTime;
		maxLength = neatLength/avgTime;
		remain = 0;
	end
	
	
	%% Create empty structs for each averaged variable with a new structure field for each spiral
	
	%%% Flight-level variables
	tempC_avg.(spiralStr) = zeros(maxLength,1)*NaN;
	RH_avg.(spiralStr) = zeros(maxLength,1)*NaN;
	alt_avg.(spiralStr) = zeros(maxLength,1)*NaN;
	time_secsFL_avg.(spiralStr) = zeros(maxLength,1)*NaN;
	
	%%% Accepted particles
	area_ratio_avg.(spiralStr) = zeros(maxLength,1)*NaN;
	reject_ratio_avg.(spiralStr) = zeros(maxLength,1)*NaN;
	efct_rad_avg.(spiralStr) = zeros(maxLength,1)*NaN;
	n_avg.(spiralStr) = zeros(maxLength,1)*NaN;
	
	iwc_avg.(spiralStr) = zeros(maxLength,1)*NaN;
	lwc_avg.(spiralStr) = zeros(maxLength,1)*NaN;
	twc_avg.(spiralStr) = zeros(maxLength,1)*NaN;
	
	conc_minR_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
	mean_perim_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
	mean_areaRatio_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
	mean_aspectRatio_elps_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
	mean_aspectRatio_rect_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
	count_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
	area_calcd_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
	prec_rate_ice_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
	prec_rate_lw_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
	termVeloc_ice_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
	termVeloc_lw_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
	massBL_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
	mass_ice_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
	mass_lw_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
	total_area_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
	conc_AreaR_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
	
	area_avg.(spiralStr) = zeros(maxLength,num_bins,10)*NaN;
	sd_habit_avg.(spiralStr) = zeros(maxLength,num_bins,10)*NaN;
	sdMass_habit_avg.(spiralStr) = zeros(maxLength,num_bins,10)*NaN;
	
	%%% Rejected particles
	if calcRej
		rej_area_ratio_avg.(spiralStr) = zeros(maxLength,1)*NaN;
		rej_efct_rad_avg.(spiralStr) = zeros(maxLength,1)*NaN;
		rej_n_avg.(spiralStr) = zeros(maxLength,1)*NaN;
		
		rej_conc_minR_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
		rej_mean_perim_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
		rej_mean_areaRatio_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
		rej_mean_aspectRatio_elps_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
		rej_mean_aspectRatio_rect_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
		rej_count_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
		rej_area_calcd_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
		rej_prec_rate_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
		rej_termVeloc_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
		rej_massBL_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
		rej_mass_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
		rej_total_area_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
		rej_conc_AreaR_avg.(spiralStr) = zeros(maxLength,num_bins)*NaN;
		
		rej_area_avg.(spiralStr) = zeros(maxLength,num_bins,10)*NaN;
		rej_sd_habit_avg.(spiralStr) = zeros(maxLength,num_bins,10)*NaN;
		rej_sdMass_habit_avg.(spiralStr) = zeros(maxLength,num_bins,10)*NaN;
	end
	
	
	
	%% Create variables containing original 1-sec data for each spiral
	
	%%% Flight-level variables
	tempC_orig.(spiralStr) = tempC_all(sprlLocs_FL);
	RH_orig.(spiralStr) = RH_all(sprlLocs_FL);
	alt_orig.(spiralStr) = alt_all(sprlLocs_FL);
	time_secsFL_orig.(spiralStr) = time_secsFL_all(sprlLocs_FL);
	
	%%% Accepted particles
	area_ratio_orig.(spiralStr) = area_ratio_all(sprlLocs);
	reject_ratio_orig.(spiralStr) = reject_ratio_all(sprlLocs);
	efct_rad_orig.(spiralStr) = efct_rad_all(sprlLocs);
	n_orig.(spiralStr) = n_all(sprlLocs);
	
	iwc_orig.(spiralStr) = iwc_all(sprlLocs);
	lwc_orig.(spiralStr) = lwc_all(sprlLocs);
	
	conc_minR_orig.(spiralStr) = conc_minR_all(sprlLocs,:);
	mean_perim_orig.(spiralStr) = mean_perim_all(sprlLocs,:);
	mean_areaRatio_orig.(spiralStr) = mean_areaRatio_all(sprlLocs,:);
	mean_aspectRatio_elps_orig.(spiralStr) = mean_aspectRatio_elps_all(sprlLocs,:);
	mean_aspectRatio_rect_orig.(spiralStr) = mean_aspectRatio_rect_all(sprlLocs,:);
	count_orig.(spiralStr) = count_all(sprlLocs,:);
	sampleVol_orig.(spiralStr) = sampleVol_all(sprlLocs,:);
	area_calcd_orig.(spiralStr) = area_calcd_all(sprlLocs,:);
	prec_rate_ice_orig.(spiralStr) = prec_rate_ice_all(sprlLocs,:);
	prec_rate_lw_orig.(spiralStr) = prec_rate_lw_all(sprlLocs,:);
	termVeloc_ice_orig.(spiralStr) = termVeloc_ice_all(sprlLocs,:);
	termVeloc_lw_orig.(spiralStr) = termVeloc_lw_all(sprlLocs,:);
	massBL_orig.(spiralStr) = massBL_all(sprlLocs,:);
	mass_ice_orig.(spiralStr) = mass_ice_all(sprlLocs,:);
	mass_lw_orig.(spiralStr) = mass_lw_all(sprlLocs,:);
	total_area_orig.(spiralStr) = total_area_all(sprlLocs,:);
	conc_AreaR_orig.(spiralStr) = conc_AreaR_all(sprlLocs,:);
	
	area_orig.(spiralStr) = area_all(sprlLocs,:,:);
	sd_habit_orig.(spiralStr) = sd_habit_all(sprlLocs,:,:);
	sdMass_habit_orig.(spiralStr) = sdMass_habit_all(sprlLocs,:,:);
	
	time_secs_orig.(spiralStr) = time_secs_all(sprlLocs);
	
	
	%%% Rejected particles
	if calcRej
		rej_area_ratio_orig.(spiralStr) = rej_area_ratio_all(sprlLocs);
		rej_efct_rad_orig.(spiralStr) = rej_efct_rad_all(sprlLocs);
		rej_n_orig.(spiralStr) = rej_n_all(sprlLocs);
		
		rej_conc_minR_orig.(spiralStr) = rej_conc_minR_all(sprlLocs,:);
		rej_mean_perim_orig.(spiralStr) = rej_mean_perim_all(sprlLocs,:);
		rej_mean_areaRatio_orig.(spiralStr) = rej_mean_areaRatio_all(sprlLocs,:);
		rej_mean_aspectRatio_elps_orig.(spiralStr) = rej_mean_aspectRatio_elps_all(sprlLocs,:);
		rej_mean_aspectRatio_rect_orig.(spiralStr) = rej_mean_aspectRatio_rect_all(sprlLocs,:);
		rej_count_orig.(spiralStr) = rej_count_all(sprlLocs,:);
		rej_area_calcd_orig.(spiralStr) = rej_area_calcd_all(sprlLocs,:);
		rej_prec_rate_orig.(spiralStr) = rej_prec_rate_all(sprlLocs,:);
		rej_termVeloc_orig.(spiralStr) = rej_termVeloc_all(sprlLocs,:);
		rej_massBL_orig.(spiralStr) = rej_massBL_all(sprlLocs,:);
		rej_mass_orig.(spiralStr) = rej_mass_all(sprlLocs,:);
		rej_total_area_orig.(spiralStr) = rej_total_area_all(sprlLocs,:);
		rej_conc_AreaR_orig.(spiralStr) = rej_conc_AreaR_all(sprlLocs,:);
		
		rej_area_orig.(spiralStr) = rej_area_all(sprlLocs,:,:);
		rej_sd_habit_orig.(spiralStr) = rej_sd_habit_all(sprlLocs,:,:);
		rej_sdMass_habit_orig.(spiralStr) = rej_sdMass_habit_all(sprlLocs,:,:);
	end
	
	%% Calculate total water content and median mass diameter 
	% Melting layer bottom (defined visually as the time when first evidence of
	% ice was observed) is used to differentiate where IWC and LWC are used to define
	% total water content and median mass diameter.
	% From bottom of spiral up to, but not including the melting layer bottom, define
	% TWC as LWC, and IWC everywhere else.
	
	for ii=1:length(mass_ice_orig.(spiralStr))
		massIceOrigTmp(ii,:)=mass_ice_orig.(spiralStr)(ii,:).*(bin_size')/10; % [g/cm3]
		massLwOrigTmp(ii,:)=mass_lw_orig.(spiralStr)(ii,:).*(bin_size')/10; % [g/cm3]
	end
	iwcOrigTmp = iwc_orig.(spiralStr)/1e6;
	lwcOrigTmp = lwc_orig.(spiralStr)/1e6;
	
	Dmm_ice_orig.(spiralStr) = calc_mmd(bin_mid,massIceOrigTmp,iwcOrigTmp); % [mm]
	Dmm_ice_orig.(spiralStr)(Dmm_ice_orig.(spiralStr) == 0) = NaN;
	
	Dmm_lw_orig.(spiralStr) = calc_mmd(bin_mid,massLwOrigTmp,lwcOrigTmp); % [mm]
	Dmm_lw_orig.(spiralStr)(Dmm_lw_orig.(spiralStr) == 0) = NaN;
	
	
	twc_orig.(spiralStr) = iwc_orig.(spiralStr);
	Dmm_twc_orig.(spiralStr) = Dmm_ice_orig.(spiralStr);
	mass_twc_orig.(spiralStr) = mass_ice_orig.(spiralStr);
	
	% If we have defined times/temps for the melting layer bottom of the current spiral,
	% we redefine TWC/Dmm/Mass_TWC to use liquid water values below the melting layer bottom
	% (NaN's are specified in the PECAN parameter file variables for periods where a given variable
	% is currently undefined)
	if ~isnan(mlBotTime(ix))
		botIx = find(time_secs_orig.(spiralStr) == mlBotTime(ix));
		if tempC_orig.(spiralStr)(1) < tempC_orig.(spiralStr)(end) % Spiral down
			twc_orig.(spiralStr)(botIx+1:end) = lwc_orig.(spiralStr)(botIx+1:end);
			Dmm_twc_orig.(spiralStr)(botIx+1:end) = Dmm_lw_orig.(spiralStr)(botIx+1:end);
			mass_twc_orig.(spiralStr)(botIx+1:end,:) = mass_lw_orig.(spiralStr)(botIx+1:end,:);
		else % Spiral up
			twc_orig.(spiralStr)(1:botIx) = lwc_orig.(spiralStr)(1:botIx);
			Dmm_twc_orig.(spiralStr)(1:botIx) = Dmm_lw_orig.(spiralStr)(1:botIx);
			mass_twc_orig.(spiralStr)(1:botIx,:) = mass_lw_orig.(spiralStr)(1:botIx,:);
		end
	end
	
	clearvars massIceOrigTmp massLwOrigTmp
	
	
	
	
	%% Calculate averages for neatLength period
	i = 1;
	for iz=1:avgTime:neatLength
		%%% Flight-level variables
		tempC_avg.(spiralStr)(i) = nansum(tempC_orig.(spiralStr)(iz:iz+(avgTime-1)),1)/avgTime;
		RH_avg.(spiralStr)(i) = nansum(RH_orig.(spiralStr)(iz:iz+(avgTime-1)),1)/avgTime;
		alt_avg.(spiralStr)(i) = nansum(alt_orig.(spiralStr)(iz:iz+(avgTime-1)),1)/avgTime;
		time_secsFL_avg.(spiralStr)(i) = nansum(time_secsFL_orig.(spiralStr)(iz:iz+(avgTime-1)),1)/avgTime;
		
		iwc_avg.(spiralStr)(i) = nansum(iwc_orig.(spiralStr)(iz:iz+(avgTime-1)),1)/avgTime;
		lwc_avg.(spiralStr)(i) = nansum(lwc_orig.(spiralStr)(iz:iz+(avgTime-1)),1)/avgTime;
		twc_avg.(spiralStr)(i) = nansum(twc_orig.(spiralStr)(iz:iz+(avgTime-1)),1)/avgTime;
		
		%%% Accepted particles
		area_ratio_avg.(spiralStr)(i) = nansum(area_ratio_orig.(spiralStr)(iz:iz+(avgTime-1)),1)/avgTime;
		reject_ratio_avg.(spiralStr)(i) = nansum(reject_ratio_orig.(spiralStr)(iz:iz+(avgTime-1)),1)/avgTime;
		efct_rad_avg.(spiralStr)(i) = nansum(efct_rad_orig.(spiralStr)(iz:iz+(avgTime-1)),1)/avgTime;
		n_avg.(spiralStr)(i) = nansum(n_orig.(spiralStr)(iz:iz+(avgTime-1)),1)/avgTime;
		
		conc_minR_avg.(spiralStr)(i,:) = nansum(conc_minR_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
		mean_perim_avg.(spiralStr)(i,:) = nansum(mean_perim_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
		mean_areaRatio_avg.(spiralStr)(i,:) = nansum(mean_areaRatio_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
		mean_aspectRatio_elps_avg.(spiralStr)(i,:) = nansum(mean_aspectRatio_elps_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
		mean_aspectRatio_rect_avg.(spiralStr)(i,:) = nansum(mean_aspectRatio_rect_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
		count_avg.(spiralStr)(i,:) = nansum(count_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
		sampleVol_avg.(spiralStr)(i,:) = nansum(sampleVol_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
		area_calcd_avg.(spiralStr)(i,:) = nansum(area_calcd_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
		prec_rate_ice_avg.(spiralStr)(i,:) = nansum(prec_rate_ice_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
		termVeloc_ice_avg.(spiralStr)(i,:) = nansum(termVeloc_ice_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
		massBL_avg.(spiralStr)(i,:) = nansum(massBL_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
		mass_ice_avg.(spiralStr)(i,:) = nansum(mass_ice_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
		total_area_avg.(spiralStr)(i,:) = nansum(total_area_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
		conc_AreaR_avg.(spiralStr)(i,:) = nansum(conc_AreaR_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
		
		area_avg.(spiralStr)(i,:,:) = nansum(area_orig.(spiralStr)(iz:iz+(avgTime-1),:,:),1)/avgTime;
		sd_habit_avg.(spiralStr)(i,:,:) = nansum(sd_habit_orig.(spiralStr)(iz:iz+(avgTime-1),:,:),1)/avgTime;
		sdMass_habit_avg.(spiralStr)(i,:,:) = nansum(sdMass_habit_orig.(spiralStr)(iz:iz+(avgTime-1),:,:),1)/avgTime;
		
		time_secs_avg.(spiralStr)(i) = nansum(time_secs_orig.(spiralStr)(iz:iz+(avgTime-1)),1)/avgTime;
		
		%%% Rejected particles
		if calcRej
			rej_area_ratio_avg.(spiralStr)(i) = nansum(rej_area_ratio_orig.(spiralStr)(iz:iz+(avgTime-1)),1)/avgTime;
			rej_efct_rad_avg.(spiralStr)(i) = nansum(rej_efct_rad_orig.(spiralStr)(iz:iz+(avgTime-1)),1)/avgTime;
			rej_n_avg.(spiralStr)(i) = nansum(rej_n_orig.(spiralStr)(iz:iz+(avgTime-1)),1)/avgTime;

			rej_conc_minR_avg.(spiralStr)(i,:) = nansum(rej_conc_minR_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
			rej_mean_perim_avg.(spiralStr)(i,:) = nansum(rej_mean_perim_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
			rej_mean_areaRatio_avg.(spiralStr)(i,:) = nansum(rej_mean_areaRatio_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
			rej_mean_aspectRatio_elps_avg.(spiralStr)(i,:) = nansum(rej_mean_aspectRatio_elps_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
			rej_mean_aspectRatio_rect_avg.(spiralStr)(i,:) = nansum(rej_mean_aspectRatio_rect_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
			rej_count_avg.(spiralStr)(i,:) = nansum(rej_count_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
			rej_area_calcd_avg.(spiralStr)(i,:) = nansum(rej_area_calcd_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
			rej_prec_rate_avg.(spiralStr)(i,:) = nansum(rej_prec_rate_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
			rej_termVeloc_avg.(spiralStr)(i,:) = nansum(rej_termVeloc_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
			rej_massBL_avg.(spiralStr)(i,:) = nansum(rej_massBL_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
			rej_mass_avg.(spiralStr)(i,:) = nansum(rej_mass_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
			rej_total_area_avg.(spiralStr)(i,:) = nansum(rej_total_area_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;
			rej_conc_AreaR_avg.(spiralStr)(i,:) = nansum(rej_conc_AreaR_orig.(spiralStr)(iz:iz+(avgTime-1),:),1)/avgTime;

			rej_area_avg.(spiralStr)(i,:,:) = nansum(rej_area_orig.(spiralStr)(iz:iz+(avgTime-1),:,:),1)/avgTime;
			rej_sd_habit_avg.(spiralStr)(i,:,:) = nansum(rej_sd_habit_orig.(spiralStr)(iz:iz+(avgTime-1),:,:),1)/avgTime;
			rej_sdMass_habit_avg.(spiralStr)(i,:,:) = nansum(rej_sdMass_habit_orig.(spiralStr)(iz:iz+(avgTime-1),:,:),1)/avgTime;
		end

		i=i+1;
	end
	
	%% Calculate averages for any times in the remainder period (if one exists)
	if remain
		%%% Flight-level variables
		tempC_avg.(spiralStr)(maxLength) = nansum(tempC_orig.(spiralStr)(neatLength+1:end),1)/leftover;
		RH_avg.(spiralStr)(maxLength) = nansum(RH_orig.(spiralStr)(neatLength+1:end),1)/leftover;
		alt_avg.(spiralStr)(maxLength) = nansum(alt_orig.(spiralStr)(neatLength+1:end),1)/leftover;
		time_secsFL_avg.(spiralStr)(maxLength) = nansum(time_secsFL_orig.(spiralStr)(neatLength+1:end),1)/leftover;
		
		iwc_avg.(spiralStr)(maxLength) = nansum(iwc_orig.(spiralStr)(neatLength+1:end),1)/leftover;
		lwc_avg.(spiralStr)(maxLength) = nansum(lwc_orig.(spiralStr)(neatLength+1:end),1)/leftover;
		twc_avg.(spiralStr)(maxLength) = nansum(twc_orig.(spiralStr)(neatLength+1:end),1)/leftover;
		
		%%% Accepted particles
		area_ratio_avg.(spiralStr)(maxLength) = nansum(area_ratio_orig.(spiralStr)(neatLength+1:end),1)/leftover;
		reject_ratio_avg.(spiralStr)(maxLength) = nansum(reject_ratio_orig.(spiralStr)(neatLength+1:end),1)/leftover;
		efct_rad_avg.(spiralStr)(maxLength) = nansum(efct_rad_orig.(spiralStr)(neatLength+1:end),1)/leftover;
		n_avg.(spiralStr)(maxLength) = nansum(n_orig.(spiralStr)(neatLength+1:end),1)/leftover;
		
		conc_minR_avg.(spiralStr)(maxLength,:) = nansum(conc_minR_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
		mean_perim_avg.(spiralStr)(maxLength,:) = nansum(mean_perim_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
		mean_areaRatio_avg.(spiralStr)(maxLength,:) = nansum(mean_areaRatio_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
		mean_aspectRatio_elps_avg.(spiralStr)(maxLength,:) = nansum(mean_aspectRatio_elps_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
		mean_aspectRatio_rect_avg.(spiralStr)(maxLength,:) = nansum(mean_aspectRatio_rect_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
		count_avg.(spiralStr)(maxLength,:) = nansum(count_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
		sampleVol_avg.(spiralStr)(maxLength,:) = nansum(sampleVol_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
		area_calcd_avg.(spiralStr)(maxLength,:) = nansum(area_calcd_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
		prec_rate_ice_avg.(spiralStr)(maxLength,:) = nansum(prec_rate_ice_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
		termVeloc_ice_avg.(spiralStr)(maxLength,:) = nansum(termVeloc_ice_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
		massBL_avg.(spiralStr)(maxLength,:) = nansum(massBL_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
		mass_ice_avg.(spiralStr)(maxLength,:) = nansum(mass_ice_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
		total_area_avg.(spiralStr)(maxLength,:) = nansum(total_area_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
		conc_AreaR_avg.(spiralStr)(maxLength,:) = nansum(conc_AreaR_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
		
		area_avg.(spiralStr)(maxLength,:,:) = nansum(area_orig.(spiralStr)(neatLength+1:end,:,:),1)/leftover;
		sd_habit_avg.(spiralStr)(maxLength,:,:) = nansum(sd_habit_orig.(spiralStr)(neatLength+1:end,:,:),1)/leftover;
		sdMass_habit_avg.(spiralStr)(maxLength,:,:) = nansum(sdMass_habit_orig.(spiralStr)(neatLength+1:end,:,:),1)/leftover;
		
		time_secs_avg.(spiralStr)(maxLength) = nansum(time_secs_orig.(spiralStr)(neatLength+1:end),1)/leftover;
		
		
		%%% Rejected particles
		if calcRej
			rej_area_ratio_avg.(spiralStr)(maxLength) = nansum(rej_area_ratio_orig.(spiralStr)(neatLength+1:end),1)/leftover;
			rej_efct_rad_avg.(spiralStr)(maxLength) = nansum(rej_efct_rad_orig.(spiralStr)(neatLength+1:end),1)/leftover;
			rej_n_avg.(spiralStr)(maxLength) = nansum(rej_n_orig.(spiralStr)(neatLength+1:end),1)/leftover;

			rej_conc_minR_avg.(spiralStr)(maxLength,:) = nansum(rej_conc_minR_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
			rej_mean_perim_avg.(spiralStr)(maxLength,:) = nansum(rej_mean_perim_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
			rej_mean_areaRatio_avg.(spiralStr)(maxLength,:) = nansum(rej_mean_areaRatio_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
			rej_mean_aspectRatio_elps_avg.(spiralStr)(maxLength,:) = nansum(rej_mean_aspectRatio_elps_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
			rej_mean_aspectRatio_rect_avg.(spiralStr)(maxLength,:) = nansum(rej_mean_aspectRatio_rect_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
			rej_count_avg.(spiralStr)(maxLength,:) = nansum(rej_count_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
			rej_area_calcd_avg.(spiralStr)(maxLength,:) = nansum(rej_area_calcd_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
			rej_prec_rate_avg.(spiralStr)(maxLength,:) = nansum(rej_prec_rate_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
			rej_termVeloc_avg.(spiralStr)(maxLength,:) = nansum(rej_termVeloc_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
			rej_massBL_avg.(spiralStr)(maxLength,:) = nansum(rej_massBL_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
			rej_mass_avg.(spiralStr)(maxLength,:) = nansum(rej_mass_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
			rej_total_area_avg.(spiralStr)(maxLength,:) = nansum(rej_total_area_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;
			rej_conc_AreaR_avg.(spiralStr)(maxLength,:) = nansum(rej_conc_AreaR_orig.(spiralStr)(neatLength+1:end,:),1)/leftover;

			rej_area_avg.(spiralStr)(maxLength,:,:) = nansum(rej_area_orig.(spiralStr)(neatLength+1:end,:,:),1)/leftover;
			rej_sd_habit_avg.(spiralStr)(maxLength,:,:) = nansum(rej_sd_habit_orig.(spiralStr)(neatLength+1:end,:,:),1)/leftover;
			rej_sdMass_habit_avg.(spiralStr)(maxLength,:,:) = nansum(rej_sdMass_habit_orig.(spiralStr)(neatLength+1:end,:,:),1)/leftover;
		end
	end
	
	
	for iiii=1:size(mass_ice_avg.(spiralStr),1)
		massIceAvgTmp(iiii,:)=mass_ice_avg.(spiralStr)(iiii,:).*(bin_size')/10; % [g/cm3]
		massLwAvgTmp(iiii,:)=mass_lw_avg.(spiralStr)(iiii,:).*(bin_size')/10; % [g/cm3]
	end
	iwcAvgTmp = iwc_avg.(spiralStr)/1e6;
	lwcAvgTmp = lwc_avg.(spiralStr)/1e6;
	
	Dmm_ice_avg.(spiralStr) = calc_mmd(bin_mid,massIceAvgTmp,iwcAvgTmp); % [mm]
	Dmm_ice_avg.(spiralStr)(Dmm_ice_avg.(spiralStr) == 0) = NaN;
	Dmm_lw_avg.(spiralStr) = calc_mmd(bin_mid,massLwAvgTmp,lwcAvgTmp); % [mm]
	Dmm_lw_avg.(spiralStr)(Dmm_lw_avg.(spiralStr) == 0) = NaN;
	
	
	twc_avg.(spiralStr) = iwc_avg.(spiralStr);
	Dmm_twc_avg.(spiralStr) = Dmm_ice_avg.(spiralStr);
	mass_twc_avg.(spiralStr) = mass_ice_avg.(spiralStr);
	
	% If we have defined times/temps for the melting layer bottom of the current spiral,
	% we redefine TWC/Dmm/Mass_TWC to use liquid water values below the melting layer bottom
	% (NaN's are specified in the PECAN parameter file variables for periods where a given variable
	% is currently undefined)
	if ~isnan(mlBotTime(ix))
		botIx = find(time_secs_avg.(spiralStr) == mlBotTime(ix));
		if tempC_avg.(spiralStr)(1) < tempC_avg.(spiralStr)(end) % Spiral down
			twc_avg.(spiralStr)(botIx+1:end) = lwc_avg.(spiralStr)(botIx+1:end);
			Dmm_twc_avg.(spiralStr)(botIx+1:end) = Dmm_lw_avg.(spiralStr)(botIx+1:end);
			mass_twc_avg.(spiralStr)(botIx+1:end,:) = mass_lw_avg.(spiralStr)(botIx+1:end,:);
		else % Spiral up
			twc_avg.(spiralStr)(1:botIx-1) = lwc_avg.(spiralStr)(1:botIx-1);
			Dmm_twc_avg.(spiralStr)(1:botIx-1) = Dmm_lw_avg.(spiralStr)(1:botIx-1);
			mass_twc_avg.(spiralStr)(1:botIx-1,:) = mass_lw_avg.(spiralStr)(1:botIx-1,:);
		end
	end
	
	clearvars massIceAvgTmp massLwAvgTmp
	
end


%% Save output file(s)

clearvars flFile FLstr hhFL i ii iii iiii ix iz jj leftover maxLength mmFL neatLength...
	num_bins remain sDistFile spiralStr sprlLocs ssFL calcRej botIx startT endT importVars...
	mlBotTime 
	
save([dataPath 'mp-data/' flight '/sDist/' 'sdistCI.' flight '.' probe '.' num2str(avgTime)	 'secAvg' outFileAppend '.mat']);
% save([dataPath 'mp-data/' flight '/sDist_matchBins/' 'sdistCI.' flight '.' probe '.' num2str(avgTime)	 'secAvg' outFileAppend '.mat']);
