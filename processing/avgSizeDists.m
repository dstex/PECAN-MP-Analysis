% Compute 5- or 10-second averaged size distributions and flight-level
% data for each spiral of a particular flight

clearvars

flight = '20150617';

probe = 'CIP';

avgTime = 10; % sec

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

startT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'startT');
endT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'endT');
mlBotTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlBotTime');

% outFileAppend = '_Dgt150um';
outFileAppend = '';

%% Define filenames
sDistFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.' probe '.cdf'];

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


time_secs_all = hhmmss2insec(timehhmmss);

num_bins = length(bin_mid);


%% Change 0's to NaNs in terminal velocity and effective radius
% Old version of sizeDist incorrectly set periods with no particles to 0
% for these vars, which does not make sense
efct_rad_all(efct_rad_all == 0) = NaN;
termVeloc_ice_all(termVeloc_ice_all == 0) = NaN;
termVeloc_lw_all(termVeloc_lw_all == 0) = NaN;



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

% Determine which times have all NaNs for binned IWC/LWC
nan_iwc = find(all(isnan(iwc_all_a),2));
nan_lwc = find(all(isnan(lwc_all_a),2));

iwc_all = nansum(iwc_all_a,2)*1e6; % Sum all values in all bins and convert to g/m3
iwc_all(nan_iwc) = NaN; % nansum incorrectly sets sums of all NaNs to 0 - fix this here

% Liquid water content
lwc_all = nansum(lwc_all_a,2)*1e6;
lwc_all(nan_lwc) = NaN; % nansum incorrectly sets sums of all NaNs to 0 - fix this here

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
	tempC_avg.(spiralStr) = NaN(maxLength,1);
	RH_avg.(spiralStr) = NaN(maxLength,1);
	alt_avg.(spiralStr) = NaN(maxLength,1);
	time_secsFL_avg.(spiralStr) = NaN(maxLength,1);
	
	%%% Accepted particles
	area_ratio_avg.(spiralStr) = NaN(maxLength,1);
	reject_ratio_avg.(spiralStr) = NaN(maxLength,1);
	efct_rad_avg.(spiralStr) = NaN(maxLength,1);
	n_avg.(spiralStr) = NaN(maxLength,1);
	
	iwc_avg.(spiralStr) = NaN(maxLength,1);
	lwc_avg.(spiralStr) = NaN(maxLength,1);
	twc_avg.(spiralStr) = NaN(maxLength,1);
	
	conc_minR_avg.(spiralStr) = NaN(maxLength,num_bins);
	mean_perim_avg.(spiralStr) = NaN(maxLength,num_bins);
	mean_areaRatio_avg.(spiralStr) = NaN(maxLength,num_bins);
	mean_aspectRatio_elps_avg.(spiralStr) = NaN(maxLength,num_bins);
	mean_aspectRatio_rect_avg.(spiralStr) = NaN(maxLength,num_bins);
	count_avg.(spiralStr) = NaN(maxLength,num_bins);
	area_calcd_avg.(spiralStr) = NaN(maxLength,num_bins);
	prec_rate_ice_avg.(spiralStr) = NaN(maxLength,num_bins);
	prec_rate_lw_avg.(spiralStr) = NaN(maxLength,num_bins);
	termVeloc_ice_avg.(spiralStr) = NaN(maxLength,num_bins);
	termVeloc_lw_avg.(spiralStr) = NaN(maxLength,num_bins);
	massBL_avg.(spiralStr) = NaN(maxLength,num_bins);
	mass_ice_avg.(spiralStr) = NaN(maxLength,num_bins);
	mass_lw_avg.(spiralStr) = NaN(maxLength,num_bins);
	total_area_avg.(spiralStr) = NaN(maxLength,num_bins);
	conc_AreaR_avg.(spiralStr) = NaN(maxLength,num_bins);
	
	area_avg.(spiralStr) = zeros(maxLength,num_bins,10)*NaN;
	sd_habit_avg.(spiralStr) = zeros(maxLength,num_bins,10)*NaN;
	sdMass_habit_avg.(spiralStr) = zeros(maxLength,num_bins,10)*NaN;
	
	
	
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
		% Check to see if all values in average period are NaN and ensure the averaged value is also NaN (not 0)
		if ~all(isnan(tempC_orig.(spiralStr)(iz:iz+(avgTime-1))))
			tempC_avg.(spiralStr)(i) = nanmean(tempC_orig.(spiralStr)(iz:iz+(avgTime-1)),1);
		end
		if ~all(isnan(RH_orig.(spiralStr)(iz:iz+(avgTime-1))))
			RH_avg.(spiralStr)(i) = nanmean(RH_orig.(spiralStr)(iz:iz+(avgTime-1)),1);
		end
		if ~all(isnan(alt_orig.(spiralStr)(iz:iz+(avgTime-1))))
			alt_avg.(spiralStr)(i) = nanmean(alt_orig.(spiralStr)(iz:iz+(avgTime-1)),1);
		end
		if ~all(isnan(time_secsFL_orig.(spiralStr)(iz:iz+(avgTime-1))))
			time_secsFL_avg.(spiralStr)(i) = nanmean(time_secsFL_orig.(spiralStr)(iz:iz+(avgTime-1)),1);
		end
		
		if ~all(isnan(iwc_orig.(spiralStr)(iz:iz+(avgTime-1))))
			iwc_avg.(spiralStr)(i) = nanmean(iwc_orig.(spiralStr)(iz:iz+(avgTime-1)),1);
		end
		if ~all(isnan(lwc_orig.(spiralStr)(iz:iz+(avgTime-1))))
			lwc_avg.(spiralStr)(i) = nanmean(lwc_orig.(spiralStr)(iz:iz+(avgTime-1)),1);
		end
		if ~all(isnan(twc_orig.(spiralStr)(iz:iz+(avgTime-1))))
			twc_avg.(spiralStr)(i) = nanmean(twc_orig.(spiralStr)(iz:iz+(avgTime-1)),1);
		end
		
		if ~all(isnan(area_ratio_orig.(spiralStr)(iz:iz+(avgTime-1))))
			area_ratio_avg.(spiralStr)(i) = nanmean(area_ratio_orig.(spiralStr)(iz:iz+(avgTime-1)),1);
		end
		if ~all(isnan(reject_ratio_orig.(spiralStr)(iz:iz+(avgTime-1))))
			reject_ratio_avg.(spiralStr)(i) = nanmean(reject_ratio_orig.(spiralStr)(iz:iz+(avgTime-1)),1);
		end
		if ~all(isnan(efct_rad_orig.(spiralStr)(iz:iz+(avgTime-1))))
			efct_rad_avg.(spiralStr)(i) = nanmean(efct_rad_orig.(spiralStr)(iz:iz+(avgTime-1)),1);
		end
		if ~all(isnan(n_orig.(spiralStr)(iz:iz+(avgTime-1))))
			n_avg.(spiralStr)(i) = nanmean(n_orig.(spiralStr)(iz:iz+(avgTime-1)),1);
		end
		if ~all(isnan(time_secs_orig.(spiralStr)(iz:iz+(avgTime-1))))
			time_secs_avg.(spiralStr)(i) = nanmean(time_secs_orig.(spiralStr)(iz:iz+(avgTime-1)),1);
		end
		
		for ib = 1:length(bin_mid)
			if ~all(isnan(conc_minR_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				conc_minR_avg.(spiralStr)(i,ib) = nanmean(conc_minR_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(mean_perim_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				mean_perim_avg.(spiralStr)(i,ib) = nanmean(mean_perim_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(mean_areaRatio_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				mean_areaRatio_avg.(spiralStr)(i,ib) = nanmean(mean_areaRatio_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(mean_aspectRatio_elps_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				mean_aspectRatio_elps_avg.(spiralStr)(i,ib) = nanmean(mean_aspectRatio_elps_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(mean_aspectRatio_rect_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				mean_aspectRatio_rect_avg.(spiralStr)(i,ib) = nanmean(mean_aspectRatio_rect_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(count_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				count_avg.(spiralStr)(i,ib) = nanmean(count_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(sampleVol_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				sampleVol_avg.(spiralStr)(i,ib) = nanmean(sampleVol_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(area_calcd_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				area_calcd_avg.(spiralStr)(i,ib) = nanmean(area_calcd_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(prec_rate_ice_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				prec_rate_ice_avg.(spiralStr)(i,ib) = nanmean(prec_rate_ice_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(termVeloc_ice_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				termVeloc_ice_avg.(spiralStr)(i,ib) = nanmean(termVeloc_ice_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(massBL_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				massBL_avg.(spiralStr)(i,ib) = nanmean(massBL_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(mass_ice_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				mass_ice_avg.(spiralStr)(i,ib) = nanmean(mass_ice_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(total_area_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				total_area_avg.(spiralStr)(i,ib) = nanmean(total_area_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(conc_AreaR_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				conc_AreaR_avg.(spiralStr)(i,ib) = nanmean(conc_AreaR_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end

			for ih = 1:10
				if ~all(isnan(area_orig.(spiralStr)(iz:iz+(avgTime-1),ib,ih)))
					area_avg.(spiralStr)(i,ib,ih) = nanmean(area_orig.(spiralStr)(iz:iz+(avgTime-1),ib,ih),1);
				end
				if ~all(isnan(sd_habit_orig.(spiralStr)(iz:iz+(avgTime-1),ib,ih)))
					sd_habit_avg.(spiralStr)(i,ib,ih) = nanmean(sd_habit_orig.(spiralStr)(iz:iz+(avgTime-1),ib,ih),1);
				end
				if ~all(isnan(sdMass_habit_orig.(spiralStr)(iz:iz+(avgTime-1),ib,ih)))
					sdMass_habit_avg.(spiralStr)(i,ib,ih) = nanmean(sdMass_habit_orig.(spiralStr)(iz:iz+(avgTime-1),ib,ih),1);
				end
			end
		end

		i=i+1;
	end
	
	%% Calculate averages for any times in the remainder period (if one exists)
	if remain
		%%% Flight-level variables
		if ~all(isnan(tempC_orig.(spiralStr)(neatLength+1:end)))
			tempC_avg.(spiralStr)(maxLength) = nanmean(tempC_orig.(spiralStr)(neatLength+1:end),1);
		end
		if ~all(isnan(RH_orig.(spiralStr)(neatLength+1:end)))
			RH_avg.(spiralStr)(maxLength) = nanmean(RH_orig.(spiralStr)(neatLength+1:end),1);
		end
		if ~all(isnan(alt_orig.(spiralStr)(neatLength+1:end)))
			alt_avg.(spiralStr)(maxLength) = nanmean(alt_orig.(spiralStr)(neatLength+1:end),1);
		end
		if ~all(isnan(time_secsFL_orig.(spiralStr)(neatLength+1:end)))
			time_secsFL_avg.(spiralStr)(maxLength) = nanmean(time_secsFL_orig.(spiralStr)(neatLength+1:end),1);
		end
		
		if ~all(isnan(iwc_orig.(spiralStr)(neatLength+1:end)))
			iwc_avg.(spiralStr)(maxLength) = nanmean(iwc_orig.(spiralStr)(neatLength+1:end),1);
		end
		if ~all(isnan(lwc_orig.(spiralStr)(neatLength+1:end)))
			lwc_avg.(spiralStr)(maxLength) = nanmean(lwc_orig.(spiralStr)(neatLength+1:end),1);
		end
		if ~all(isnan(twc_orig.(spiralStr)(neatLength+1:end)))
			twc_avg.(spiralStr)(maxLength) = nanmean(twc_orig.(spiralStr)(neatLength+1:end),1);
		end
		
		if ~all(isnan(area_ratio_orig.(spiralStr)(neatLength+1:end)))
			area_ratio_avg.(spiralStr)(maxLength) = nanmean(area_ratio_orig.(spiralStr)(neatLength+1:end),1);
		end
		if ~all(isnan(reject_ratio_orig.(spiralStr)(neatLength+1:end)))
			reject_ratio_avg.(spiralStr)(maxLength) = nanmean(reject_ratio_orig.(spiralStr)(neatLength+1:end),1);
		end
		if ~all(isnan(efct_rad_orig.(spiralStr)(neatLength+1:end)))
			efct_rad_avg.(spiralStr)(maxLength) = nanmean(efct_rad_orig.(spiralStr)(neatLength+1:end),1);
		end
		if ~all(isnan(n_orig.(spiralStr)(neatLength+1:end)))
			n_avg.(spiralStr)(maxLength) = nanmean(n_orig.(spiralStr)(neatLength+1:end),1);
		end
		if ~all(isnan(time_secs_orig.(spiralStr)(neatLength+1:end)))
			time_secs_avg.(spiralStr)(maxLength) = nanmean(time_secs_orig.(spiralStr)(neatLength+1:end),1);
		end
		
		for ib = 1:length(bin_mid)
			if ~all(isnan(conc_minR_orig.(spiralStr)(neatLength+1:end,ib)))
				conc_minR_avg.(spiralStr)(maxLength,ib) = nanmean(conc_minR_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(mean_perim_orig.(spiralStr)(neatLength+1:end,ib)))
				mean_perim_avg.(spiralStr)(maxLength,ib) = nanmean(mean_perim_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(mean_areaRatio_orig.(spiralStr)(neatLength+1:end,ib)))
				mean_areaRatio_avg.(spiralStr)(maxLength,ib) = nanmean(mean_areaRatio_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(mean_aspectRatio_elps_orig.(spiralStr)(neatLength+1:end,ib)))
				mean_aspectRatio_elps_avg.(spiralStr)(maxLength,ib) = nanmean(mean_aspectRatio_elps_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(mean_aspectRatio_rect_orig.(spiralStr)(neatLength+1:end,ib)))
				mean_aspectRatio_rect_avg.(spiralStr)(maxLength,ib) = nanmean(mean_aspectRatio_rect_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(count_orig.(spiralStr)(neatLength+1:end,ib)))
				count_avg.(spiralStr)(maxLength,ib) = nanmean(count_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(sampleVol_orig.(spiralStr)(neatLength+1:end,ib)))
				sampleVol_avg.(spiralStr)(maxLength,ib) = nanmean(sampleVol_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(area_calcd_orig.(spiralStr)(neatLength+1:end,ib)))
				area_calcd_avg.(spiralStr)(maxLength,ib) = nanmean(area_calcd_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(prec_rate_ice_orig.(spiralStr)(neatLength+1:end,ib)))
				prec_rate_ice_avg.(spiralStr)(maxLength,ib) = nanmean(prec_rate_ice_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(termVeloc_ice_orig.(spiralStr)(neatLength+1:end,ib)))
				termVeloc_ice_avg.(spiralStr)(maxLength,ib) = nanmean(termVeloc_ice_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(massBL_orig.(spiralStr)(neatLength+1:end,ib)))
				massBL_avg.(spiralStr)(maxLength,ib) = nanmean(massBL_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(mass_ice_orig.(spiralStr)(neatLength+1:end,ib)))
				mass_ice_avg.(spiralStr)(maxLength,ib) = nanmean(mass_ice_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(total_area_orig.(spiralStr)(neatLength+1:end,ib)))
				total_area_avg.(spiralStr)(maxLength,ib) = nanmean(total_area_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(conc_AreaR_orig.(spiralStr)(neatLength+1:end,ib)))
				conc_AreaR_avg.(spiralStr)(maxLength,ib) = nanmean(conc_AreaR_orig.(spiralStr)(neatLength+1:end,ib),1);
			end

			for ih = 1:10
				if ~all(isnan(area_orig.(spiralStr)(neatLength+1:end,ib,ih)))
					area_avg.(spiralStr)(maxLength,ib,ih) = nanmean(area_orig.(spiralStr)(neatLength+1:end,ib,ih),1);
				end
				if ~all(isnan(sd_habit_orig.(spiralStr)(neatLength+1:end,ib,ih)))
					sd_habit_avg.(spiralStr)(maxLength,ib,ih) = nanmean(sd_habit_orig.(spiralStr)(neatLength+1:end,ib,ih),1);
				end
				if ~all(isnan(sdMass_habit_orig.(spiralStr)(neatLength+1:end,ib,ih)))
					sdMass_habit_avg.(spiralStr)(maxLength,ib,ih) = nanmean(sdMass_habit_orig.(spiralStr)(neatLength+1:end,ib,ih),1);
				end
			end
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
	num_bins remain sDistFile spiralStr sprlLocs ssFL botIx startT endT importVars...
	mlBotTime 
	
save([dataPath 'mp-data/' flight '/sDist/' 'sdistCI.' flight '.' probe '.' num2str(avgTime)	 'secAvg' outFileAppend '.mat']);
% save([dataPath 'mp-data/' flight '/sDist_matchBins/' 'sdistCI.' flight '.' probe '.' num2str(avgTime)	 'secAvg' outFileAppend '.mat']);
