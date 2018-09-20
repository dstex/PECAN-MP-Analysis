% Compute 5- or 10-second averaged size distributions and flight-level
% data for each spiral of a particular flight

clearvars

flight = '20150620';

probe = 'CIP';

avgTime = 10; % sec

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';

startT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'startT');
endT = nc_varget([dataPath '/' flight '_PECANparams.nc'],'endT');
mlBotTime = nc_varget([dataPath '/' flight '_PECANparams.nc'],'mlBotTime');

outFileAppend = '';

%% Define filenames
sDistFile = [dataPath 'mp-data/' flight '/sDist/sdistCI.' flight '.' probe '.cdf'];

flFile = [dataPath 'FlightLevelData/Processed/' flight '_FltLvl_Processed.mat'];


%% Import all size distribution variables

timehhmmss = nc_varget(sDistFile,'time');

bin_min_mm = nc_varget(sDistFile,'bin_min'); % All bin_* in mm
bin_max_mm = nc_varget(sDistFile,'bin_max');
bin_mid_mm = nc_varget(sDistFile,'bin_mid');
bin_size_mm = nc_varget(sDistFile,'bin_dD');

% conc_minR_cm4_all = (nc_varget(sDistFile,'conc_minR'))'; % cm-4
% area_cm4_all = permute(nc_varget(sDistFile,'area'),[3 2 1]);
% area_calcd_mm2cm4_all = (nc_varget(sDistFile,'Calcd_area'))';
% conc_AreaR_cm4_all = (nc_varget(sDistFile,'conc_AreaR'))';
% n_cm3_all = nc_varget(sDistFile,'n');
% total_area_mm2cm4_all = (nc_varget(sDistFile,'total_area'))';
% mass_ice_gcm4_all = (nc_varget(sDistFile,'mass_ice'))'; % g cm-4
% mass_lw_gcm4_all = (nc_varget(sDistFile,'mass_lw'))'; % g cm-4
% massBL_gcm4_all = (nc_varget(sDistFile,'massBL'))';
% sd_habit_cm4_all = permute(nc_varget(sDistFile,'habitsd'),[3 2 1]);
% sdMass_habit_gcm4_all = permute(nc_varget(sDistFile,'habitmsd'),[3 2 1]);
% efct_rad_um_all = nc_varget(sDistFile,'re');
% area_ratio_all = nc_varget(sDistFile,'ar');
% reject_ratio_all = nc_varget(sDistFile,'Reject_ratio');
% termVeloc_ice_all = (nc_varget(sDistFile,'vt_ice'))';
% termVeloc_lw_all = (nc_varget(sDistFile,'vt_lw'))';
% prec_rate_ice_mmhr_all = (nc_varget(sDistFile,'Prec_rate_ice'))';
% prec_rate_lw_mmhr_all = (nc_varget(sDistFile,'Prec_rate_lw'))';
% count_all = (nc_varget(sDistFile,'count'))';
% sampleVol_cm3_all = (nc_varget(sDistFile,'sample_vol'))';
% 
% mean_aspectRatio_rect_all = (nc_varget(sDistFile,'mean_aspect_ratio_rectangle'))';
% mean_aspectRatio_elps_all = (nc_varget(sDistFile,'mean_aspect_ratio_ellipse'))';
% mean_areaRatio_all = (nc_varget(sDistFile,'mean_area_ratio'))';
% mean_perim_um_all = (nc_varget(sDistFile,'mean_perimeter'))';

% For some reason different versions of MATLAB handle the dimension order
% differently on import... Some require the above, others require this:
conc_minR_cm4_all = nc_varget(sDistFile,'conc_minR'); % cm-4
area_cm4_all = nc_varget(sDistFile,'area');
area_calcd_mm2cm4_all = nc_varget(sDistFile,'Calcd_area');
conc_AreaR_cm4_all = nc_varget(sDistFile,'conc_AreaR');
n_cm3_all = nc_varget(sDistFile,'n');
total_area_mm2cm4_all = nc_varget(sDistFile,'total_area');
mass_ice_gcm4_all = nc_varget(sDistFile,'mass_ice'); % g cm-4
mass_lw_gcm4_all = nc_varget(sDistFile,'mass_lw'); % g cm-4
massBL_gcm4_all = nc_varget(sDistFile,'massBL');
sd_habit_cm4_all = nc_varget(sDistFile,'habitsd');
sdMass_habit_gcm4_all = nc_varget(sDistFile,'habitmsd');
efct_rad_um_all = nc_varget(sDistFile,'re');
area_ratio_all = nc_varget(sDistFile,'ar');
reject_ratio_all = nc_varget(sDistFile,'Reject_ratio');
termVeloc_ice_all = nc_varget(sDistFile,'vt_ice');
termVeloc_lw_all = nc_varget(sDistFile,'vt_lw');
prec_rate_ice_mmhr_all = nc_varget(sDistFile,'Prec_rate_ice');
prec_rate_lw_mmhr_all = nc_varget(sDistFile,'Prec_rate_lw');
count_all = nc_varget(sDistFile,'count');
sampleVol_cm3_all = nc_varget(sDistFile,'sample_vol');

mean_aspectRatio_rect_all = nc_varget(sDistFile,'mean_aspect_ratio_rectangle');
mean_aspectRatio_elps_all = nc_varget(sDistFile,'mean_aspect_ratio_ellipse');
mean_areaRatio_all = nc_varget(sDistFile,'mean_area_ratio');
mean_perim_um_all = nc_varget(sDistFile,'mean_perimeter');


time_secs_all = hhmmss2insec(timehhmmss);

num_bins = length(bin_mid_mm);



%% Import relevant flight-level data
importVars = {'time_secs_FL','TA','RH_hybrid','Alt'};

tempLoad = load(flFile,importVars{:});

time_secsFL_all = tempLoad.(importVars{1});
tempC_all = tempLoad.(importVars{2});
RH_all = tempLoad.(importVars{3});
alt_all = tempLoad.(importVars{4});


%% Compute any variables not available in the sDist output

iwc_all_gcm3 = mass_ice_gcm4_all.*(bin_size_mm'/10); % Multiply by bin size in cm [now g/cm3]
lwc_all_gcm3 = mass_lw_gcm4_all.*(bin_size_mm'/10);

% Determine which times have all NaNs for binned IWC/LWC
nan_iwc = find(all(isnan(iwc_all_gcm3),2));
nan_lwc = find(all(isnan(lwc_all_gcm3),2));

iwc_gm3_all = nansum(iwc_all_gcm3,2)*1e6; % Sum all values in all bins and convert to g/m3
iwc_gm3_all(nan_iwc) = NaN; % nansum incorrectly sets sums of all NaNs to 0 - fix this here

% Liquid water content
lwc_gm3_all = nansum(lwc_all_gcm3,2)*1e6;
lwc_gm3_all(nan_lwc) = NaN;

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
	efct_rad_um_avg.(spiralStr) = NaN(maxLength,1);
	n_cm3_avg.(spiralStr) = NaN(maxLength,1);
	
	iwc_gm3_avg.(spiralStr) = NaN(maxLength,1);
	lwc_gm3_avg.(spiralStr) = NaN(maxLength,1);
	twc_gm3_avg.(spiralStr) = NaN(maxLength,1);
	
	conc_minR_cm4_avg.(spiralStr) = NaN(maxLength,num_bins);
	mean_perim_um_avg.(spiralStr) = NaN(maxLength,num_bins);
	mean_areaRatio_avg.(spiralStr) = NaN(maxLength,num_bins);
	mean_aspectRatio_elps_avg.(spiralStr) = NaN(maxLength,num_bins);
	mean_aspectRatio_rect_avg.(spiralStr) = NaN(maxLength,num_bins);
	count_avg.(spiralStr) = NaN(maxLength,num_bins);
	sampleVol_cm3_avg.(spiralStr) = NaN(maxLength,num_bins);
	area_calcd_mm2cm4_avg.(spiralStr) = NaN(maxLength,num_bins);
	prec_rate_ice_mmhr_avg.(spiralStr) = NaN(maxLength,num_bins);
	prec_rate_lw_mmhr_avg.(spiralStr) = NaN(maxLength,num_bins);
	termVeloc_ice_avg.(spiralStr) = NaN(maxLength,num_bins);
	termVeloc_lw_avg.(spiralStr) = NaN(maxLength,num_bins);
	massBL_gcm4_avg.(spiralStr) = NaN(maxLength,num_bins);
	mass_ice_gcm4_avg.(spiralStr) = NaN(maxLength,num_bins);
	mass_lw_gcm4_avg.(spiralStr) = NaN(maxLength,num_bins);
	total_area_mm2cm4_avg.(spiralStr) = NaN(maxLength,num_bins);
	conc_AreaR_cm4_avg.(spiralStr) = NaN(maxLength,num_bins);
	
	area_cm4_avg.(spiralStr) = zeros(maxLength,num_bins,10)*NaN;
	sd_habit_cm4_avg.(spiralStr) = zeros(maxLength,num_bins,10)*NaN;
	sdMass_habit_gcm4_avg.(spiralStr) = zeros(maxLength,num_bins,10)*NaN;
	
	
	
	%% Create variables containing original 1-sec data for each spiral
	
	%%% Flight-level variables
	tempC_orig.(spiralStr) = tempC_all(sprlLocs_FL);
	RH_orig.(spiralStr) = RH_all(sprlLocs_FL);
	alt_orig.(spiralStr) = alt_all(sprlLocs_FL);
	time_secsFL_orig.(spiralStr) = time_secsFL_all(sprlLocs_FL);
	
	%%% Accepted particles
	area_ratio_orig.(spiralStr) = area_ratio_all(sprlLocs);
	reject_ratio_orig.(spiralStr) = reject_ratio_all(sprlLocs);
	efct_rad_um_orig.(spiralStr) = efct_rad_um_all(sprlLocs);
	n_cm3_orig.(spiralStr) = n_cm3_all(sprlLocs);
	
	iwc_gm3_orig.(spiralStr) = iwc_gm3_all(sprlLocs); % [g m-3]
	
	conc_minR_cm4_orig.(spiralStr) = conc_minR_cm4_all(sprlLocs,:); % [cm-4]
	mean_perim_um_orig.(spiralStr) = mean_perim_um_all(sprlLocs,:);
	mean_areaRatio_orig.(spiralStr) = mean_areaRatio_all(sprlLocs,:);
	mean_aspectRatio_elps_orig.(spiralStr) = mean_aspectRatio_elps_all(sprlLocs,:);
	mean_aspectRatio_rect_orig.(spiralStr) = mean_aspectRatio_rect_all(sprlLocs,:);
	count_orig.(spiralStr) = count_all(sprlLocs,:);
	sampleVol_cm3_orig.(spiralStr) = sampleVol_cm3_all(sprlLocs,:);
	area_calcd_mm2cm4_orig.(spiralStr) = area_calcd_mm2cm4_all(sprlLocs,:);
	prec_rate_ice_mmhr_orig.(spiralStr) = prec_rate_ice_mmhr_all(sprlLocs,:);
	prec_rate_lw_mmhr_orig.(spiralStr) = prec_rate_lw_mmhr_all(sprlLocs,:);
	termVeloc_ice_orig.(spiralStr) = termVeloc_ice_all(sprlLocs,:);
	termVeloc_lw_orig.(spiralStr) = termVeloc_lw_all(sprlLocs,:);
	massBL_gcm4_orig.(spiralStr) = massBL_gcm4_all(sprlLocs,:);
	mass_ice_gcm4_orig.(spiralStr) = mass_ice_gcm4_all(sprlLocs,:); %[g cm-4]
	mass_lw_gcm4_orig.(spiralStr) = mass_lw_gcm4_all(sprlLocs,:);
	total_area_mm2cm4_orig.(spiralStr) = total_area_mm2cm4_all(sprlLocs,:);
	conc_AreaR_cm4_orig.(spiralStr) = conc_AreaR_cm4_all(sprlLocs,:);
	
	area_cm4_orig.(spiralStr) = area_cm4_all(sprlLocs,:,:);
	sd_habit_cm4_orig.(spiralStr) = sd_habit_cm4_all(sprlLocs,:,:);
	sdMass_habit_gcm4_orig.(spiralStr) = sdMass_habit_gcm4_all(sprlLocs,:,:);
	
	time_secs_orig.(spiralStr) = time_secs_all(sprlLocs);
	
	
	%% Calculate total water content and median mass diameter 
	% Melting layer bottom (defined visually as the time when first evidence of
	% ice was observed) is used to differentiate where IWC and LWC are used to define
	% total water content and median mass diameter.
	% From bottom of spiral up to, but not including the melting layer bottom, define
	% TWC as LWC, and IWC everywhere else.

	if ~isnan(mlBotTime(ix))
		[~, botIx] = min(abs(time_secs_orig.(spiralStr) - mlBotTime(ix)));
		if tempC_orig.(spiralStr)(1) < tempC_orig.(spiralStr)(end) % Spiral down
			sprlMeanAspR_lw_orig.(spiralStr) = nanmean(nanmean(mean_aspectRatio_elps_orig.(spiralStr)(botIx+1:end,:),1));
			mass_lw_gcm4_orig.(spiralStr)(botIx+1:end,:) = sprlMeanAspR_lw_orig.(spiralStr).*mass_lw_gcm4_orig.(spiralStr)(botIx+1:end,:);
		else % spiral up
			sprlMeanAspR_lw_orig.(spiralStr) = nanmean(nanmean(mean_aspectRatio_elps_orig.(spiralStr)(1:botIx-1,:),1));
			mass_lw_gcm4_orig.(spiralStr)(1:botIx-1,:) = sprlMeanAspR_lw_orig.(spiralStr).*mass_lw_gcm4_orig.(spiralStr)(1:botIx-1,:);
		end
	end
	
	nan_lwc_o = find(all(isnan(mass_lw_gcm4_orig.(spiralStr)),2));
	
	lwc_gm3_orig.(spiralStr) = nansum(mass_lw_gcm4_orig.(spiralStr).*(bin_size_mm'/10),2)*1e6; % [g m-3]
	lwc_gm3_orig.(spiralStr)(nan_lwc_o) = NaN;
	
	Dmm_ice_mm_orig.(spiralStr) = calc_mmd(bin_mid_mm,mass_ice_gcm4_orig.(spiralStr).*(bin_size_mm'./10),iwc_gm3_orig.(spiralStr)./1e6); % [mm]
	Dmm_ice_mm_orig.(spiralStr)(Dmm_ice_mm_orig.(spiralStr) == 0) = NaN;
	
	Dmm_lw_mm_orig.(spiralStr) = calc_mmd(bin_mid_mm,mass_lw_gcm4_orig.(spiralStr).*(bin_size_mm'./10),lwc_gm3_orig.(spiralStr)./1e6); % [mm]
	Dmm_lw_mm_orig.(spiralStr)(Dmm_lw_mm_orig.(spiralStr) == 0) = NaN;
	
	
	twc_gm3_orig.(spiralStr) = iwc_gm3_orig.(spiralStr); % [g m-3]
	Dmm_twc_mm_orig.(spiralStr) = Dmm_ice_mm_orig.(spiralStr); % [mm]
	mass_twc_gcm4_orig.(spiralStr) = mass_ice_gcm4_orig.(spiralStr); % [g cm-4]
	
	% If we have defined times/temps for the melting layer bottom of the current spiral,
	% we redefine TWC/Dmm/Mass_TWC to use liquid water values below the melting layer bottom
	% (NaN's are specified in the PECAN parameter file variables for periods where a given variable
	% is currently undefined)
	ice_flag_orig.(spiralStr) = ones(length(time_secs_orig.(spiralStr)),1); % Boolean array - 1=above/in ML; 0=below ML
	if ~isnan(mlBotTime(ix))
		if tempC_orig.(spiralStr)(1) < tempC_orig.(spiralStr)(end) % Spiral down
			twc_gm3_orig.(spiralStr)(botIx+1:end) = lwc_gm3_orig.(spiralStr)(botIx+1:end);
			Dmm_twc_mm_orig.(spiralStr)(botIx+1:end) = Dmm_lw_mm_orig.(spiralStr)(botIx+1:end);
			mass_twc_gcm4_orig.(spiralStr)(botIx+1:end,:) = mass_lw_gcm4_orig.(spiralStr)(botIx+1:end,:);
			ice_flag_orig.(spiralStr)(botIx+1:end) = 0;
		else % Spiral up
			twc_gm3_orig.(spiralStr)(1:botIx-1) = lwc_gm3_orig.(spiralStr)(1:botIx-1);
			Dmm_twc_mm_orig.(spiralStr)(1:botIx-1) = Dmm_lw_mm_orig.(spiralStr)(1:botIx-1);
			mass_twc_gcm4_orig.(spiralStr)(1:botIx-1,:) = mass_lw_gcm4_orig.(spiralStr)(1:botIx-1,:);
			ice_flag_orig.(spiralStr)(1:botIx-1) = 0;
		end
	end
	
	ZD_orig.(spiralStr) = (bin_mid_mm.^6)'.* ( (conc_minR_cm4_orig.(spiralStr).*(bin_size_mm'./10)).*1e6 ); %[mm^6] * [cm^-4 * cm] * [1e^6 cm^3 m^-3] --> [mm^6 m^-3]
	dBZ_orig.(spiralStr) = 10.*log10(nansum(ZD_orig.(spiralStr),2));
	
	
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
		
		if ~all(isnan(area_ratio_orig.(spiralStr)(iz:iz+(avgTime-1))))
			area_ratio_avg.(spiralStr)(i) = nanmean(area_ratio_orig.(spiralStr)(iz:iz+(avgTime-1)),1);
		end
		if ~all(isnan(reject_ratio_orig.(spiralStr)(iz:iz+(avgTime-1))))
			reject_ratio_avg.(spiralStr)(i) = nanmean(reject_ratio_orig.(spiralStr)(iz:iz+(avgTime-1)),1);
		end
		if ~all(isnan(efct_rad_um_orig.(spiralStr)(iz:iz+(avgTime-1))))
			efct_rad_um_avg.(spiralStr)(i) = nanmean(efct_rad_um_orig.(spiralStr)(iz:iz+(avgTime-1)),1);
		end

		if ~all(isnan(time_secs_orig.(spiralStr)(iz:iz+(avgTime-1))))
			time_secs_avg.(spiralStr)(i) = nanmean(time_secs_orig.(spiralStr)(iz:iz+(avgTime-1)),1);
		end
		
		for ib = 1:length(bin_mid_mm)
			if ~all(isnan(conc_minR_cm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				conc_minR_cm4_avg.(spiralStr)(i,ib) = nanmean(conc_minR_cm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(mean_perim_um_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				mean_perim_um_avg.(spiralStr)(i,ib) = nanmean(mean_perim_um_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
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
			if ~all(isnan(sampleVol_cm3_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				sampleVol_cm3_avg.(spiralStr)(i,ib) = nanmean(sampleVol_cm3_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(area_calcd_mm2cm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				area_calcd_mm2cm4_avg.(spiralStr)(i,ib) = nanmean(area_calcd_mm2cm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(prec_rate_ice_mmhr_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				prec_rate_ice_mmhr_avg.(spiralStr)(i,ib) = nanmean(prec_rate_ice_mmhr_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(termVeloc_ice_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				termVeloc_ice_avg.(spiralStr)(i,ib) = nanmean(termVeloc_ice_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(massBL_gcm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				massBL_gcm4_avg.(spiralStr)(i,ib) = nanmean(massBL_gcm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(mass_ice_gcm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				mass_ice_gcm4_avg.(spiralStr)(i,ib) = nanmean(mass_ice_gcm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(mass_lw_gcm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				mass_lw_gcm4_avg.(spiralStr)(i,ib) = nanmean(mass_lw_gcm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(total_area_mm2cm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				total_area_mm2cm4_avg.(spiralStr)(i,ib) = nanmean(total_area_mm2cm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end
			if ~all(isnan(conc_AreaR_cm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib)))
				conc_AreaR_cm4_avg.(spiralStr)(i,ib) = nanmean(conc_AreaR_cm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib),1);
			end

			for ih = 1:10
				if ~all(isnan(area_cm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib,ih)))
					area_cm4_avg.(spiralStr)(i,ib,ih) = nanmean(area_cm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib,ih),1);
				end
				if ~all(isnan(sd_habit_cm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib,ih)))
					sd_habit_cm4_avg.(spiralStr)(i,ib,ih) = nanmean(sd_habit_cm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib,ih),1);
				end
				if ~all(isnan(sdMass_habit_gcm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib,ih)))
					sdMass_habit_gcm4_avg.(spiralStr)(i,ib,ih) = nanmean(sdMass_habit_gcm4_orig.(spiralStr)(iz:iz+(avgTime-1),ib,ih),1);
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
		
		if ~all(isnan(area_ratio_orig.(spiralStr)(neatLength+1:end)))
			area_ratio_avg.(spiralStr)(maxLength) = nanmean(area_ratio_orig.(spiralStr)(neatLength+1:end),1);
		end
		if ~all(isnan(reject_ratio_orig.(spiralStr)(neatLength+1:end)))
			reject_ratio_avg.(spiralStr)(maxLength) = nanmean(reject_ratio_orig.(spiralStr)(neatLength+1:end),1);
		end
		if ~all(isnan(efct_rad_um_orig.(spiralStr)(neatLength+1:end)))
			efct_rad_um_avg.(spiralStr)(maxLength) = nanmean(efct_rad_um_orig.(spiralStr)(neatLength+1:end),1);
		end

		if ~all(isnan(time_secs_orig.(spiralStr)(neatLength+1:end)))
			time_secs_avg.(spiralStr)(maxLength) = nanmean(time_secs_orig.(spiralStr)(neatLength+1:end),1);
		end
		
		for ib = 1:length(bin_mid_mm)
			if ~all(isnan(conc_minR_cm4_orig.(spiralStr)(neatLength+1:end,ib)))
				conc_minR_cm4_avg.(spiralStr)(maxLength,ib) = nanmean(conc_minR_cm4_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(mean_perim_um_orig.(spiralStr)(neatLength+1:end,ib)))
				mean_perim_um_avg.(spiralStr)(maxLength,ib) = nanmean(mean_perim_um_orig.(spiralStr)(neatLength+1:end,ib),1);
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
			if ~all(isnan(sampleVol_cm3_orig.(spiralStr)(neatLength+1:end,ib)))
				sampleVol_cm3_avg.(spiralStr)(maxLength,ib) = nanmean(sampleVol_cm3_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(area_calcd_mm2cm4_orig.(spiralStr)(neatLength+1:end,ib)))
				area_calcd_mm2cm4_avg.(spiralStr)(maxLength,ib) = nanmean(area_calcd_mm2cm4_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(prec_rate_ice_mmhr_orig.(spiralStr)(neatLength+1:end,ib)))
				prec_rate_ice_mmhr_avg.(spiralStr)(maxLength,ib) = nanmean(prec_rate_ice_mmhr_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(termVeloc_ice_orig.(spiralStr)(neatLength+1:end,ib)))
				termVeloc_ice_avg.(spiralStr)(maxLength,ib) = nanmean(termVeloc_ice_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(massBL_gcm4_orig.(spiralStr)(neatLength+1:end,ib)))
				massBL_gcm4_avg.(spiralStr)(maxLength,ib) = nanmean(massBL_gcm4_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(mass_ice_gcm4_orig.(spiralStr)(neatLength+1:end,ib)))
				mass_ice_gcm4_avg.(spiralStr)(maxLength,ib) = nanmean(mass_ice_gcm4_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(mass_lw_gcm4_orig.(spiralStr)(neatLength+1:end,ib)))
				mass_lw_gcm4_avg.(spiralStr)(maxLength,ib) = nanmean(mass_lw_gcm4_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(total_area_mm2cm4_orig.(spiralStr)(neatLength+1:end,ib)))
				total_area_mm2cm4_avg.(spiralStr)(maxLength,ib) = nanmean(total_area_mm2cm4_orig.(spiralStr)(neatLength+1:end,ib),1);
			end
			if ~all(isnan(conc_AreaR_cm4_orig.(spiralStr)(neatLength+1:end,ib)))
				conc_AreaR_cm4_avg.(spiralStr)(maxLength,ib) = nanmean(conc_AreaR_cm4_orig.(spiralStr)(neatLength+1:end,ib),1);
			end

			for ih = 1:10
				if ~all(isnan(area_cm4_orig.(spiralStr)(neatLength+1:end,ib,ih)))
					area_cm4_avg.(spiralStr)(maxLength,ib,ih) = nanmean(area_cm4_orig.(spiralStr)(neatLength+1:end,ib,ih),1);
				end
				if ~all(isnan(sd_habit_cm4_orig.(spiralStr)(neatLength+1:end,ib,ih)))
					sd_habit_cm4_avg.(spiralStr)(maxLength,ib,ih) = nanmean(sd_habit_cm4_orig.(spiralStr)(neatLength+1:end,ib,ih),1);
				end
				if ~all(isnan(sdMass_habit_gcm4_orig.(spiralStr)(neatLength+1:end,ib,ih)))
					sdMass_habit_gcm4_avg.(spiralStr)(maxLength,ib,ih) = nanmean(sdMass_habit_gcm4_orig.(spiralStr)(neatLength+1:end,ib,ih),1);
				end
			end
		end
	end
	
	if ~isnan(mlBotTime(ix))
		[~, botIx] = min(abs(time_secs_avg.(spiralStr) - mlBotTime(ix)));
		if tempC_avg.(spiralStr)(1) < tempC_avg.(spiralStr)(end) % Spiral down
			sprlMeanAspR_lw_avg.(spiralStr) = nanmean(nanmean(mean_aspectRatio_elps_avg.(spiralStr)(botIx+1:end,:),1));
			mass_lw_gcm4_avg.(spiralStr)(botIx+1:end,:) = sprlMeanAspR_lw_avg.(spiralStr).*mass_lw_gcm4_avg.(spiralStr)(botIx+1:end,:);
		else % spiral up
			sprlMeanAspR_lw_avg.(spiralStr) = nanmean(nanmean(mean_aspectRatio_elps_avg.(spiralStr)(1:botIx-1,:),1));
			mass_lw_gcm4_avg.(spiralStr)(1:botIx-1,:) = sprlMeanAspR_lw_avg.(spiralStr).*mass_lw_gcm4_avg.(spiralStr)(1:botIx-1,:);	
		end
	end
	
	

	nan_iwc = find(all(isnan(mass_ice_gcm4_avg.(spiralStr)),2));
	nan_lwc = find(all(isnan(mass_lw_gcm4_avg.(spiralStr)),2));
	
	iwc_gm3_avg.(spiralStr) = nansum(mass_ice_gcm4_avg.(spiralStr).*(bin_size_mm'/10),2)*1e6; % [g m-3]
	iwc_gm3_avg.(spiralStr)(nan_iwc) = NaN;
	lwc_gm3_avg.(spiralStr) = nansum(mass_lw_gcm4_avg.(spiralStr).*(bin_size_mm'/10),2)*1e6; % [g m-3]
	lwc_gm3_avg.(spiralStr)(nan_lwc) = NaN;

	
	Dmm_ice_mm_avg.(spiralStr) = calc_mmd(bin_mid_mm,mass_ice_gcm4_avg.(spiralStr).*(bin_size_mm'/10),iwc_gm3_avg.(spiralStr)./1e6); % [mm]
	Dmm_ice_mm_avg.(spiralStr)(Dmm_ice_mm_avg.(spiralStr) == 0) = NaN;
	Dmm_lw_mm_avg.(spiralStr) = calc_mmd(bin_mid_mm,mass_lw_gcm4_avg.(spiralStr).*(bin_size_mm'/10),lwc_gm3_avg.(spiralStr)./1e6); % [mm]
	Dmm_lw_mm_avg.(spiralStr)(Dmm_lw_mm_avg.(spiralStr) == 0) = NaN;
	
	
	twc_gm3_avg.(spiralStr) = iwc_gm3_avg.(spiralStr); % [g m-3]
	Dmm_twc_mm_avg.(spiralStr) = Dmm_ice_mm_avg.(spiralStr); % [mm]
	mass_twc_gcm4_avg.(spiralStr) = mass_ice_gcm4_avg.(spiralStr); % [g cm-4]
	
	% If we have defined times/temps for the melting layer bottom of the current spiral,
	% we redefine TWC/Dmm/Mass_TWC to use liquid water values below the melting layer bottom
	% (NaN's are specified in the PECAN parameter file variables for periods where a given variable
	% is currently undefined)
	ice_flag_avg.(spiralStr) = ones(length(time_secs_avg.(spiralStr)),1); % Boolean array - 1=above/in ML; 0=below ML
	if ~isnan(mlBotTime(ix))
		if tempC_avg.(spiralStr)(1) < tempC_avg.(spiralStr)(end) % Spiral down
			twc_gm3_avg.(spiralStr)(botIx+1:end) = lwc_gm3_avg.(spiralStr)(botIx+1:end);
			Dmm_twc_mm_avg.(spiralStr)(botIx+1:end) = Dmm_lw_mm_avg.(spiralStr)(botIx+1:end);
			mass_twc_gcm4_avg.(spiralStr)(botIx+1:end,:) = mass_lw_gcm4_avg.(spiralStr)(botIx+1:end,:);
			ice_flag_avg.(spiralStr)(botIx+1:end) = 0;
		else % Spiral up
			twc_gm3_avg.(spiralStr)(1:botIx-1) = lwc_gm3_avg.(spiralStr)(1:botIx-1);
			Dmm_twc_mm_avg.(spiralStr)(1:botIx-1) = Dmm_lw_mm_avg.(spiralStr)(1:botIx-1);
			mass_twc_gcm4_avg.(spiralStr)(1:botIx-1,:) = mass_lw_gcm4_avg.(spiralStr)(1:botIx-1,:);
			ice_flag_avg.(spiralStr)(1:botIx-1) = 0;
		end
	end
	
	ZD_avg.(spiralStr) = (bin_mid_mm.^6)'.* ( (conc_minR_cm4_avg.(spiralStr).*(bin_size_mm'./10)).*1e6 ); %[mm^6] * [cm^-4 * cm] * [1e^6 cm^3 m^-3] --> [mm^6 m^-3]
	dBZ_avg.(spiralStr) = 10.*log10(nansum(ZD_avg.(spiralStr),2));
	
	n_cm3_avg.(spiralStr) = nansum(conc_minR_cm4_avg.(spiralStr).*(bin_size_mm'/10),2); % cm-3
end


%% Save output file(s)

clearvars botIx i ib ih importVars ix iz leftover iwc_all_a lwc_all_a maxLength mlBotTime...
	nan_iwc nan_lwc nan_lwc_o neatLength num_bins remain sDistFile spiralStr sprlLocs sprlLocs_FL...
	startT endT tempLoad

save(sprintf('%smp-data/%s/sDist/sdistCI.%s.%s.%dsecAvg%s.mat',dataPath,flight,flight,probe,avgTime,outFileAppend),'-regexp',...
			'^(?!dataPath$|flFile$|avgTime$|outFileAppend$|probe$)\w');