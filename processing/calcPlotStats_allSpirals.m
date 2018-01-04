clearvars; close all;

%% Specify various plotting/calculation parameters
flights = {'20150617','20150620','20150701','20150702','20150706','20150709'};


avgTime = 10;

fileIdStr = ['_Fit-CIP_' num2str(avgTime) 'secAvg_1.2cm'];

dataPath = '/Users/danstechman/GoogleDrive/PECAN-Data/';


%% Specify plotting parameters

plotNtSprdF_G1	= 1;
plotNtSprdF_G2	= 0;
plotNtSprdF_All = 0;

cntrLine = 'Median'; % Do we plot median or mean on spread plots?
% cntrLine = 'Mean';

saveFigs	= 0;
noDisp		= 0;
Ftype		= '-dpdf';
% Ftype		= '-dpng';

if strcmp(Ftype,'-dpdf')
	Fres = '-painters'; % Use this for fully vectorized files
else
	Fres = '-r0'; % Use this for smaller files - saves figure with same resolution/size as displayed on screen
end


savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';

%% Create directories to save plots in if they don't already exist
if saveFigs
    saveDir = [savePath flight '/CIP-Fitting'];
	if (exist(saveDir, 'dir') ~= 7)
        mkdir(saveDir)
	end

	if ((plotNtSprdF_G1 || plotNtSprdF_G2 || plotNtSprdF_All) && exist([saveDir '/CIP-NtTempSprd_' num2str(avgTime) 's'], 'dir') ~= 7)
		mkdir([saveDir '/CIP-NtTempSprd_' num2str(avgTime) 's'])
	end
end

NtSprdLim = [-7 0];

%% Load in data from the original averaged PSD files and the extended distribution files
% Data will be assigned to structs, with substructs containing the data for each flight
for iFlt = 1:length(flights)
	fltStr = ['f' flights{iFlt}];
	sDistF.(fltStr) = load([dataPath 'mp-data/' flights{iFlt} '/sDist/sdistCI.' flights{iFlt} '.CIP.' num2str(avgTime) 'secAvg.mat']);
	extSDistF.(fltStr) = load([dataPath 'mp-data/' flights{iFlt} '/sDist/' flights{iFlt} fileIdStr '.mat']);
	
	startT.(fltStr) = nc_varget([dataPath '/' flights{iFlt} '_PECANparams.nc'],'startT');
	sprlZone.(fltStr) = nc_varget([dataPath '/' flights{iFlt} '_PECANparams.nc'],'sprlZone');
	mcsType.(fltStr) = nc_attget([dataPath '/' flights{iFlt} '_PECANparams.nc'],nc_global,'mcsGroup');
	
	tzSprls = find(sprlZone.(fltStr) == 'T');
	srSprls = find(sprlZone.(fltStr) == 'S');
	raSprls = find(sprlZone.(fltStr) == 'A');
	
	sprlNames = fieldnames(sDistF.(fltStr).time_secs_avg);

	
	
	if ~isempty(tzSprls)
		tzTWC_tmp = extSDistF.(fltStr).cipTWC_hybrid_igf.(sprlNames{tzSprls(1)}); % g m-3
		tzNt_tmp = extSDistF.(fltStr).cipNt_hybrid_igf.(sprlNames{tzSprls(1)}); % cm-3
		tzDmm_tmp = extSDistF.(fltStr).cipDmm_hybrid_igf.(sprlNames{tzSprls(1)}); % cm
		tzTemp_tmp = sDistF.(fltStr).tempC_avg.(sprlNames{tzSprls(1)});

		% Concatenate additional profiles from the same spiral zone
		if length(tzSprls) > 1
			for ix = 2:length(tzSprls)
				tzTWC_tmp = vertcat(tzTWC_tmp, extSDistF.(fltStr).cipTWC_hybrid_igf.(sprlNames{tzSprls(ix)}));
				tzNt_tmp = vertcat(tzNt_tmp, extSDistF.(fltStr).cipNt_hybrid_igf.(sprlNames{tzSprls(ix)}));
				tzDmm_tmp = vertcat(tzDmm_tmp, extSDistF.(fltStr).cipDmm_hybrid_igf.(sprlNames{tzSprls(ix)}));
				tzTemp_tmp = vertcat(tzTemp_tmp, sDistF.(fltStr).tempC_avg.(sprlNames{tzSprls(ix)}));
			end
		end
		
		tzTWC.(fltStr) = tzTWC_tmp;
		tzNt.(fltStr) = tzNt_tmp;
		tzDmm.(fltStr) = tzDmm_tmp;
		tzTemp.(fltStr) = tzTemp_tmp;
	end
	
	
	
	if ~isempty(srSprls)
		srTWC_tmp = extSDistF.(fltStr).cipTWC_hybrid_igf.(sprlNames{srSprls(1)}); % g m-3
		srNt_tmp = extSDistF.(fltStr).cipNt_hybrid_igf.(sprlNames{srSprls(1)}); % cm-3
		srDmm_tmp = extSDistF.(fltStr).cipDmm_hybrid_igf.(sprlNames{srSprls(1)}); % cm
		srTemp_tmp = sDistF.(fltStr).tempC_avg.(sprlNames{srSprls(1)});

		% Concatenate additional profiles from the same spiral zone
		if length(srSprls) > 1
			for ix = 2:length(srSprls)
				srTWC_tmp = vertcat(srTWC_tmp, extSDistF.(fltStr).cipTWC_hybrid_igf.(sprlNames{srSprls(ix)}));
				srNt_tmp = vertcat(srNt_tmp, extSDistF.(fltStr).cipNt_hybrid_igf.(sprlNames{srSprls(ix)}));
				srDmm_tmp = vertcat(srDmm_tmp, extSDistF.(fltStr).cipDmm_hybrid_igf.(sprlNames{srSprls(ix)}));
				srTemp_tmp = vertcat(srTemp_tmp, sDistF.(fltStr).tempC_avg.(sprlNames{srSprls(ix)}));
			end
		end
		
		srTWC.(fltStr) = srTWC_tmp;
		srNt.(fltStr) = srNt_tmp;
		srDmm.(fltStr) = srDmm_tmp;
		srTemp.(fltStr) = srTemp_tmp;
	end
	
	
	
	if ~isempty(raSprls)
		raTWC_tmp = extSDistF.(fltStr).cipTWC_hybrid_igf.(sprlNames{raSprls(1)}); % g m-3
		raNt_tmp = extSDistF.(fltStr).cipNt_hybrid_igf.(sprlNames{raSprls(1)}); % cm-3
		raDmm_tmp = extSDistF.(fltStr).cipDmm_hybrid_igf.(sprlNames{raSprls(1)}); % cm
		raTemp_tmp = sDistF.(fltStr).tempC_avg.(sprlNames{raSprls(1)});

		% Concatenate additional profiles from the same spiral zone
		if length(raSprls) > 1
			for ix = 2:length(raSprls)
				raTWC_tmp = vertcat(raTWC_tmp, extSDistF.(fltStr).cipTWC_hybrid_igf.(sprlNames{raSprls(ix)}));
				raNt_tmp = vertcat(raNt_tmp, extSDistF.(fltStr).cipNt_hybrid_igf.(sprlNames{raSprls(ix)}));
				raDmm_tmp = vertcat(raDmm_tmp, extSDistF.(fltStr).cipDmm_hybrid_igf.(sprlNames{raSprls(ix)}));
				raTemp_tmp = vertcat(raTemp_tmp, sDistF.(fltStr).tempC_avg.(sprlNames{raSprls(ix)}));
			end
		end
		
		raTWC.(fltStr) = raTWC_tmp;
		raNt.(fltStr) = raNt_tmp;
		raDmm.(fltStr) = raDmm_tmp;
		raTemp.(fltStr) = raTemp_tmp;
	end


end

fltNames = fieldnames(sDistF);


%% Make new arrays with concatenated data from each group and spiral region
% Create empty arrays
[tzTWC_g1,tzNt_g1,tzDmm_g1,tzTemp_g1,...
	tzTWC_g2,tzNt_g2,tzDmm_g2,tzTemp_g2,...
	tzTWC_all,tzNt_all,tzDmm_all,tzTemp_all] = deal([]);
[srTWC_g1,srNt_g1,srDmm_g1,srTemp_g1,...
	srTWC_g2,srNt_g2,srDmm_g2,srTemp_g2,...
	srTWC_all,srNt_all,srDmm_all,srTemp_all] = deal([]);
[raTWC_g1,raNt_g1,raDmm_g1,raTemp_g1,...
	raTWC_g2,raNt_g2,raDmm_g2,raTemp_g2,...
	raTWC_all,raNt_all,raDmm_all,raTemp_all] = deal([]);

for iFlt = 1:length(flights)
	fltStr = ['f' flights{iFlt}];
	
	tzSprls = find(sprlZone.(fltStr) == 'T');
	srSprls = find(sprlZone.(fltStr) == 'S');
	raSprls = find(sprlZone.(fltStr) == 'A');
	
	% Concatenate spiral zone variables from each flight in group 1 (TS)
	if (strcmp(fltStr,'f20150617') || strcmp(fltStr,'f20150620') || strcmp(fltStr,'f20150706'))
		if ~isempty(tzSprls)
			for ix = 1:length(tzSprls)
				tzTWC_g1 = vertcat(tzTWC_g1, extSDistF.(fltStr).cipTWC_hybrid_igf.(sprlNames{tzSprls(ix)})); % g m-3
				tzNt_g1 = vertcat(tzNt_g1, extSDistF.(fltStr).cipNt_hybrid_igf.(sprlNames{tzSprls(ix)})); % cm-3
				tzDmm_g1 = vertcat(tzDmm_g1, extSDistF.(fltStr).cipDmm_hybrid_igf.(sprlNames{tzSprls(ix)})); % cm
				tzTemp_g1 = vertcat(tzTemp_g1, sDistF.(fltStr).tempC_avg.(sprlNames{tzSprls(ix)}));
			end
			
			tzTWC.grp1 = tzTWC_g1;
			tzNt.grp1 = tzNt_g1;
			tzDmm.grp1 = tzDmm_g1;
			tzTemp.grp1 = tzTemp_g1;
		end
		
		if ~isempty(srSprls)
			for ix = 1:length(srSprls)
				srTWC_g1 = vertcat(srTWC_g1, extSDistF.(fltStr).cipTWC_hybrid_igf.(sprlNames{srSprls(ix)})); % g m-3
				srNt_g1 = vertcat(srNt_g1, extSDistF.(fltStr).cipNt_hybrid_igf.(sprlNames{srSprls(ix)})); % cm-3
				srDmm_g1 = vertcat(srDmm_g1, extSDistF.(fltStr).cipDmm_hybrid_igf.(sprlNames{srSprls(ix)})); % cm
				srTemp_g1 = vertcat(srTemp_g1, sDistF.(fltStr).tempC_avg.(sprlNames{srSprls(ix)}));
			end
			
			srTWC.grp1 = srTWC_g1;
			srNt.grp1 = srNt_g1;
			srDmm.grp1 = srDmm_g1;
			srTemp.grp1 = srTemp_g1;
		end
		
		if ~isempty(raSprls)
			for ix = 1:length(raSprls)
				raTWC_g1 = vertcat(raTWC_g1, extSDistF.(fltStr).cipTWC_hybrid_igf.(sprlNames{raSprls(ix)})); % g m-3
				raNt_g1 = vertcat(raNt_g1, extSDistF.(fltStr).cipNt_hybrid_igf.(sprlNames{raSprls(ix)})); % cm-3
				raDmm_g1 = vertcat(raDmm_g1, extSDistF.(fltStr).cipDmm_hybrid_igf.(sprlNames{raSprls(ix)})); % cm
				raTemp_g1 = vertcat(raTemp_g1, sDistF.(fltStr).tempC_avg.(sprlNames{raSprls(ix)}));
			end
			
			raTWC.grp1 = raTWC_g1;
			raNt.grp1 = raNt_g1;
			raDmm.grp1 = raDmm_g1;
			raTemp.grp1 = raTemp_g1;
		end
	end
	
	% Concatenate spiral zone variables from each flight in group 2 (PS or LS)
	if (strcmp(fltStr,'f20150701') || strcmp(fltStr,'f20150702') || strcmp(fltStr,'f20150709'))
		if ~isempty(tzSprls)
			for ix = 1:length(tzSprls)
				tzTWC_g2 = vertcat(tzTWC_g2, extSDistF.(fltStr).cipTWC_hybrid_igf.(sprlNames{tzSprls(ix)})); % g m-3
				tzNt_g2 = vertcat(tzNt_g2, extSDistF.(fltStr).cipNt_hybrid_igf.(sprlNames{tzSprls(ix)})); % cm-3
				tzDmm_g2 = vertcat(tzDmm_g2, extSDistF.(fltStr).cipDmm_hybrid_igf.(sprlNames{tzSprls(ix)})); % cm
				tzTemp_g2 = vertcat(tzTemp_g2, sDistF.(fltStr).tempC_avg.(sprlNames{tzSprls(ix)}));
			end
			
			tzTWC.grp2 = tzTWC_g2;
			tzNt.grp2 = tzNt_g2;
			tzDmm.grp2 = tzDmm_g2;
			tzTemp.grp2 = tzTemp_g2;
		end
		
		if ~isempty(srSprls)
			for ix = 1:length(srSprls)
				srTWC_g2 = vertcat(srTWC_g2, extSDistF.(fltStr).cipTWC_hybrid_igf.(sprlNames{srSprls(ix)})); % g m-3
				srNt_g2 = vertcat(srNt_g2, extSDistF.(fltStr).cipNt_hybrid_igf.(sprlNames{srSprls(ix)})); % cm-3
				srDmm_g2 = vertcat(srDmm_g2, extSDistF.(fltStr).cipDmm_hybrid_igf.(sprlNames{srSprls(ix)})); % cm
				srTemp_g2 = vertcat(srTemp_g2, sDistF.(fltStr).tempC_avg.(sprlNames{srSprls(ix)}));
			end
			
			srTWC.grp2 = srTWC_g2;
			srNt.grp2 = srNt_g2;
			srDmm.grp2 = srDmm_g2;
			srTemp.grp2 = srTemp_g2;
		end
		
		if ~isempty(raSprls)
			for ix = 1:length(raSprls)
				raTWC_g2 = vertcat(raTWC_g2, extSDistF.(fltStr).cipTWC_hybrid_igf.(sprlNames{raSprls(ix)})); % g m-3
				raNt_g2 = vertcat(raNt_g2, extSDistF.(fltStr).cipNt_hybrid_igf.(sprlNames{raSprls(ix)})); % cm-3
				raDmm_g2 = vertcat(raDmm_g2, extSDistF.(fltStr).cipDmm_hybrid_igf.(sprlNames{raSprls(ix)})); % cm
				raTemp_g2 = vertcat(raTemp_g2, sDistF.(fltStr).tempC_avg.(sprlNames{raSprls(ix)}));
			end
			
			raTWC.grp2 = raTWC_g2;
			raNt.grp2 = raNt_g2;
			raDmm.grp2 = raDmm_g2;
			raTemp.grp2 = raTemp_g2;
		end
	end
	
	
	% Combine all zones from all flights
	if ~isempty(tzSprls)
		for ix = 1:length(tzSprls)
			tzTWC_all = vertcat(tzTWC_all, extSDistF.(fltStr).cipTWC_hybrid_igf.(sprlNames{tzSprls(ix)})); % g m-3
			tzNt_all = vertcat(tzNt_all, extSDistF.(fltStr).cipNt_hybrid_igf.(sprlNames{tzSprls(ix)})); % cm-3
			tzDmm_all = vertcat(tzDmm_all, extSDistF.(fltStr).cipDmm_hybrid_igf.(sprlNames{tzSprls(ix)})); % cm
			tzTemp_all = vertcat(tzTemp_all, sDistF.(fltStr).tempC_avg.(sprlNames{tzSprls(ix)}));
		end
		
		tzTWC.all = tzTWC_all;
		tzNt.all = tzNt_all;
		tzDmm.all = tzDmm_all;
		tzTemp.all = tzTemp_all;
	end
	
	if ~isempty(srSprls)
		for ix = 1:length(srSprls)
			srTWC_all = vertcat(srTWC_all, extSDistF.(fltStr).cipTWC_hybrid_igf.(sprlNames{srSprls(ix)})); % g m-3
			srNt_all = vertcat(srNt_all, extSDistF.(fltStr).cipNt_hybrid_igf.(sprlNames{srSprls(ix)})); % cm-3
			srDmm_all = vertcat(srDmm_all, extSDistF.(fltStr).cipDmm_hybrid_igf.(sprlNames{srSprls(ix)})); % cm
			srTemp_all = vertcat(srTemp_all, sDistF.(fltStr).tempC_avg.(sprlNames{srSprls(ix)}));
		end
		
		srTWC.all = srTWC_all;
		srNt.all = srNt_all;
		srDmm.all = srDmm_all;
		srTemp.all = srTemp_all;
	end
	
	if ~isempty(raSprls)
		for ix = 1:length(raSprls)
			raTWC_all = vertcat(raTWC_all, extSDistF.(fltStr).cipTWC_hybrid_igf.(sprlNames{raSprls(ix)})); % g m-3
			raNt_all = vertcat(raNt_all, extSDistF.(fltStr).cipNt_hybrid_igf.(sprlNames{raSprls(ix)})); % cm-3
			raDmm_all = vertcat(raDmm_all, extSDistF.(fltStr).cipDmm_hybrid_igf.(sprlNames{raSprls(ix)})); % cm
			raTemp_all = vertcat(raTemp_all, sDistF.(fltStr).tempC_avg.(sprlNames{raSprls(ix)}));
		end
		
		raTWC.all = raTWC_all;
		raNt.all = raNt_all;
		raDmm.all = raDmm_all;
		raTemp.all = raTemp_all;
	end
	
end

clearvars('*_tmp','*_all','*_g1','*_g2','raSprls','srSprls','tzSprls','iFlt','ix','fltStr');

%% Create plots

if plotNtSprdF_G1
	edgesTemp = (-19.25:.5:20.25);
% 	edgesTemp = (-15.125:.25:20.125);
	numBins = length(edgesTemp)-1;

	
	if saveFigs && noDisp
		figure('visible','off','Position', [10,10,1500,1500]);
	else
		figure('Position', [10,10,1500,1500]);
	end
	hold on
	
	
	if isfield(tzTemp,'grp1') % Check to see if there were any TZ spirals in group 1 (tzTemp,tzTWC,... all will work to check)
		tzNt.grp1Edit = tzNt.grp1;
		tzNt.grp1Edit(tzNt.grp1Edit == 0) = NaN;%1e-12; % Set zeros to a very small value to allow plotting/fill to work properly
% 		tzNt.grp1Edit(isnan(tzNt.grp1Edit)) = 1e-12;
		
		tzNtBin_mid = (edgesTemp(1:end-1) + edgesTemp(2:end))/2;

		tzNtTempBinIx = discretize(tzTemp.grp1,edgesTemp);

		[tzNtFlgRmvBin,tzNtBinMean,tzNtBinMedian,...
		tzNtBinMax,tzNtBinMin,tzNtBin25p,tzNtBin75p] = deal([]);
		

		for ix=1:numBins
			tzNtBinMem = (tzNtTempBinIx == ix);
			tzBinNt = tzNt.grp1Edit(tzNtBinMem);
			tzBinNt_orig = tzNt.grp1(tzNtBinMem);
			
			if ~isempty(tzBinNt)
				tzNtBinMean = horzcat(tzNtBinMean,nanmean(tzBinNt_orig));
				tzNtBinMedian = horzcat(tzNtBinMedian,nanmedian(tzBinNt));
				tzNtBinMax = horzcat(tzNtBinMax,max(tzBinNt));
				tzNtBinMin = horzcat(tzNtBinMin,min(tzBinNt));
				tzNtBin25p = horzcat(tzNtBin25p,prctile(tzBinNt,25));
				tzNtBin75p = horzcat(tzNtBin75p,prctile(tzBinNt,75));
				
			else
				tzNtFlgRmvBin = horzcat(tzNtFlgRmvBin,ix);
			end
		end

		tzNtBin_mid(tzNtFlgRmvBin) = [];

		tzNtBin_midFlip = [tzNtBin_mid,fliplr(tzNtBin_mid)];
% 		tzNtSpread = [log10(tzNtBinMin), fliplr(log10(tzNtBinMax))];
		tzNtSpread = [log10(tzNtBin25p), fliplr(log10(tzNtBin75p))];

		f1 = fill(tzNtBin_midFlip,tzNtSpread,'g','FaceAlpha', 0.5);

		if strcmp(cntrLine,'Mean')
			m1 = plot(tzNtBin_mid,log10(tzNtBinMean),'Color',[27/255 100/255 27/255],'LineWidth', 3);
		elseif strcmp(cntrLine,'Median')
			m1 = plot(tzNtBin_mid,log10(tzNtBinMedian),'Color',[27/255 100/255 27/255],'LineWidth', 3);
		end
	end
	
	
	if isfield(srTemp,'grp1') % Check to see if there were any sr spirals in group 1 (srTemp,srTWC,... all will work to check)
		srNt.grp1Edit = srNt.grp1;
		srNt.grp1Edit(srNt.grp1Edit == 0) = NaN;%1e-12;
% 		srNt.grp1Edit(isnan(srNt.grp1Edit)) = 1e-12;
		
		srNtBin_mid = (edgesTemp(1:end-1) + edgesTemp(2:end))/2;

		srNtTempBinIx = discretize(srTemp.grp1,edgesTemp);

		[srNtFlgRmvBin,srNtBinMean,srNtBinMedian,...
		srNtBinMax,srNtBinMin,srNtBin25p,srNtBin75p] = deal([]);
		

		for ix=1:numBins
			srNtBinMem = (srNtTempBinIx == ix);
			srBinNt = srNt.grp1Edit(srNtBinMem);
			srBinNt_orig = srNt.grp1(srNtBinMem);
			
			if ~isempty(srBinNt)
				srNtBinMean = horzcat(srNtBinMean,nanmean(srBinNt_orig));
				srNtBinMedian = horzcat(srNtBinMedian,nanmedian(srBinNt));
				srNtBinMax = horzcat(srNtBinMax,max(srBinNt));
				srNtBinMin = horzcat(srNtBinMin,min(srBinNt));
				srNtBin25p = horzcat(srNtBin25p,prctile(srBinNt,25));
				srNtBin75p = horzcat(srNtBin75p,prctile(srBinNt,75));
				
			else
				srNtFlgRmvBin = horzcat(srNtFlgRmvBin,ix);
			end
		end

		srNtBin_mid(srNtFlgRmvBin) = [];

		srNtBin_midFlip = [srNtBin_mid,fliplr(srNtBin_mid)];
% 		srNtSpread = [log10(srNtBinMin), fliplr(log10(srNtBinMax))];
		srNtSpread = [log10(srNtBin25p), fliplr(log10(srNtBin75p))];

		f2 = fill(srNtBin_midFlip,srNtSpread,'r','FaceAlpha', 0.5);

		if strcmp(cntrLine,'Mean')
			m2 = plot(srNtBin_mid,log10(srNtBinMean),'Color','y','LineWidth', 3);
		elseif strcmp(cntrLine,'Median')
			m2 = plot(srNtBin_mid,log10(srNtBinMedian),'Color','y','LineWidth', 3);
		end
	end
	
	
	if isfield(raTemp,'grp1') % Check to see if there were any ra spirals in group 1 (raTemp,raTWC,... all will work to check)
		raNt.grp1Edit = raNt.grp1;
		raNt.grp1Edit(raNt.grp1Edit == 0) = NaN;%1e-12; % Set zeros to a very small value to allow plotting/fill to work properly
% 		raNt.grp1Edit(isnan(raNt.grp1Edit)) = 1e-12;
		
		
		raNtBin_mid = (edgesTemp(1:end-1) + edgesTemp(2:end))/2;

		raNtTempBinIx = discretize(raTemp.grp1,edgesTemp);

		[raNtFlgRmvBin,raNtBinMean,raNtBinMedian,...
		raNtBinMax,raNtBinMin,raNtBin25p,raNtBin75p] = deal([]);
		

		for ix=1:numBins
			raNtBinMem = (raNtTempBinIx == ix);
			raBinNt = raNt.grp1Edit(raNtBinMem);
			raBinNt_orig = raNt.grp1(raNtBinMem);
			
			if ~isempty(raBinNt)
				raNtBinMean = horzcat(raNtBinMean,nanmean(raBinNt_orig));
				raNtBinMedian = horzcat(raNtBinMedian,nanmedian(raBinNt));
				raNtBinMax = horzcat(raNtBinMax,max(raBinNt));
				raNtBinMin = horzcat(raNtBinMin,min(raBinNt));
				raNtBin25p = horzcat(raNtBin25p,prctile(raBinNt,25));
				raNtBin75p = horzcat(raNtBin75p,prctile(raBinNt,75));
				
			else
				raNtFlgRmvBin = horzcat(raNtFlgRmvBin,ix);
			end
		end

		raNtBin_mid(raNtFlgRmvBin) = [];

		raNtBin_midFlip = [raNtBin_mid,fliplr(raNtBin_mid)];
% 		raNtSpread = [log10(raNtBinMin), fliplr(log10(raNtBinMax))];
		raNtSpread = [log10(raNtBin25p), fliplr(log10(raNtBin75p))];

		f3 = fill(raNtBin_midFlip,raNtSpread,'b','FaceAlpha', 0.5);

		if strcmp(cntrLine,'Mean')
			m3 = plot(raNtBin_mid,log10(raNtBinMean),'Color',[27/255 27/255 100/255],'LineWidth', 3);
		elseif strcmp(cntrLine,'Median')
			m3 = plot(raNtBin_mid,log10(raNtBinMedian),'Color',[27/255 27/255 100/255],'LineWidth', 3);
		end
	end
	
	
	title(['Nt vs. Temp - ' cntrLine ' and Spread - CIP - All TS Spirals (Group 1)']);

	set(gca,'XMinorTick','on','YMinorTick','on','XDir','reverse');
% 	if ~isempty(tempRangeAll)
% 		xlim(tempRangeAll);
% 	end
	if ~isempty(NtSprdLim)
		ylim(NtSprdLim)
	end
	ylabel('log_{10}(N_t) [cm^{-3}]')
	xlabel(sprintf('Temp (%cC)', char(176)));
	set(findall(gcf,'-property','FontSize'),'FontSize',28);
	view([90 -90])
	
	legend([f1 f2 f3],{'Transition Zone','Stratiform Region','Rear Anvil'},'Location','eastoutside');
	
	if saveFigs
		set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		print([saveDir '/CIP-NtTempSprd_' num2str(avgTime) 's/CIP_Nt-Temp_' cntrLine '-Spread-Fill_' num2str(avgTime) 's_Grp1'],Ftype,Fres)
	end
end

if plotNtSprdF_G2
	edgesTemp = (-19.25:.5:20.25);
% 	edgesTemp = (-15.125:.25:20.125);
	numBins = length(edgesTemp)-1;

	
	if saveFigs && noDisp
		figure('visible','off','Position', [10,10,1500,1500]);
	else
		figure('Position', [10,10,1500,1500]);
	end
	hold on
	
	
	if isfield(tzTemp,'grp2') % Check to see if there were any TZ spirals in group 1 (tzTemp,tzTWC,... all will work to check)
		tzNt.grp2Edit = tzNt.grp2;
		tzNt.grp2Edit(tzNt.grp2Edit == 0) = 1e-12; % Set zeros to a very small value to allow plotting/fill to work properly
		tzNt.grp2Edit(isnan(tzNt.grp2Edit)) = 1e-12;
		
		
		tzNtBin_mid = (edgesTemp(1:end-1) + edgesTemp(2:end))/2;

		tzNtTempBinIx = discretize(tzTemp.grp2,edgesTemp);

		[tzNtFlgRmvBin,tzNtBinMean,tzNtBinMedian,...
		tzNtBinMax,tzNtBinMin,tzNtBin25p,tzNtBin75p] = deal([]);
		

		for ix=1:numBins
			tzNtBinMem = (tzNtTempBinIx == ix);
			tzBinNt = tzNt.grp2Edit(tzNtBinMem);
			tzBinNt_orig = tzNt.grp2(tzNtBinMem);
			
			if ~isempty(tzBinNt)
				tzNtBinMean = horzcat(tzNtBinMean,nanmean(tzBinNt_orig));
				tzNtBinMedian = horzcat(tzNtBinMedian,nanmedian(tzBinNt));
				tzNtBinMax = horzcat(tzNtBinMax,max(tzBinNt));
				tzNtBinMin = horzcat(tzNtBinMin,min(tzBinNt));
				tzNtBin25p = horzcat(tzNtBin25p,prctile(tzBinNt,25));
				tzNtBin75p = horzcat(tzNtBin75p,prctile(tzBinNt,75));
				
			else
				tzNtFlgRmvBin = horzcat(tzNtFlgRmvBin,ix);
			end
		end

		tzNtBin_mid(tzNtFlgRmvBin) = [];

		tzNtBin_midFlip = [tzNtBin_mid,fliplr(tzNtBin_mid)];
% 		tzNtSpread = [log10(tzNtBinMin), fliplr(log10(tzNtBinMax))];
		tzNtSpread = [log10(tzNtBin25p), fliplr(log10(tzNtBin75p))];

		f1 = fill(tzNtBin_midFlip,tzNtSpread,'g','FaceAlpha', 0.5);

		if strcmp(cntrLine,'Mean')
			m1 = plot(tzNtBin_mid,log10(tzNtBinMean),'Color',[27/255 100/255 27/255],'LineWidth', 3);
		elseif strcmp(cntrLine,'Median')
			m1 = plot(tzNtBin_mid,log10(tzNtBinMedian),'Color',[27/255 100/255 27/255],'LineWidth', 3);
		end
	end
	
	
	if isfield(srTemp,'grp2') % Check to see if there were any sr spirals in group 1 (srTemp,srTWC,... all will work to check)
		srNt.grp2Edit = srNt.grp2;
		srNt.grp2Edit(srNt.grp2Edit == 0) = 1e-12;
		srNt.grp2Edit(isnan(srNt.grp2Edit)) = 1e-12;
		
		srNtBin_mid = (edgesTemp(1:end-1) + edgesTemp(2:end))/2;

		srNtTempBinIx = discretize(srTemp.grp2,edgesTemp);

		[srNtFlgRmvBin,srNtBinMean,srNtBinMedian,...
		srNtBinMax,srNtBinMin,srNtBin25p,srNtBin75p] = deal([]);
		

		for ix=1:numBins
			srNtBinMem = (srNtTempBinIx == ix);
			srBinNt = srNt.grp2Edit(srNtBinMem);
			srBinNt_orig = srNt.grp2(srNtBinMem);
			
			if ~isempty(srBinNt)
				srNtBinMean = horzcat(srNtBinMean,nanmean(srBinNt_orig));
				srNtBinMedian = horzcat(srNtBinMedian,nanmedian(srBinNt));
				srNtBinMax = horzcat(srNtBinMax,max(srBinNt));
				srNtBinMin = horzcat(srNtBinMin,min(srBinNt));
				srNtBin25p = horzcat(srNtBin25p,prctile(srBinNt,25));
				srNtBin75p = horzcat(srNtBin75p,prctile(srBinNt,75));
				
			else
				srNtFlgRmvBin = horzcat(srNtFlgRmvBin,ix);
			end
		end

		srNtBin_mid(srNtFlgRmvBin) = [];

		srNtBin_midFlip = [srNtBin_mid,fliplr(srNtBin_mid)];
% 		srNtSpread = [log10(srNtBinMin), fliplr(log10(srNtBinMax))];
		srNtSpread = [log10(srNtBin25p), fliplr(log10(srNtBin75p))];

		f2 = fill(srNtBin_midFlip,srNtSpread,'r','FaceAlpha', 0.5);

		if strcmp(cntrLine,'Mean')
			m2 = plot(srNtBin_mid,log10(srNtBinMean),'Color','y','LineWidth', 3);
		elseif strcmp(cntrLine,'Median')
			m2 = plot(srNtBin_mid,log10(srNtBinMedian),'Color','y','LineWidth', 3);
		end
	end
	
	
	if isfield(raTemp,'grp2') % Check to see if there were any ra spirals in group 1 (raTemp,raTWC,... all will work to check)
		raNt.grp2Edit = raNt.grp2;
		raNt.grp2Edit(raNt.grp2Edit == 0) = 1e-12; % Set zeros to a very small value to allow plotting/fill to work properly
		raNt.grp2Edit(isnan(raNt.grp2Edit)) = 1e-12;
		
		raNtBin_mid = (edgesTemp(1:end-1) + edgesTemp(2:end))/2;

		raNtTempBinIx = discretize(raTemp.grp2,edgesTemp);

		[raNtFlgRmvBin,raNtBinMean,raNtBinMedian,...
		raNtBinMax,raNtBinMin,raNtBin25p,raNtBin75p] = deal([]);
		

		for ix=1:numBins
			raNtBinMem = (raNtTempBinIx == ix);
			raBinNt = raNt.grp2Edit(raNtBinMem);
			raBinNt_orig = raNt.grp2(raNtBinMem);
			
			if ~isempty(raBinNt)
				raNtBinMean = horzcat(raNtBinMean,nanmean(raBinNt_orig));
				raNtBinMedian = horzcat(raNtBinMedian,nanmedian(raBinNt));
				raNtBinMax = horzcat(raNtBinMax,max(raBinNt));
				raNtBinMin = horzcat(raNtBinMin,min(raBinNt));
				raNtBin25p = horzcat(raNtBin25p,prctile(raBinNt,25));
				raNtBin75p = horzcat(raNtBin75p,prctile(raBinNt,75));
				
			else
				raNtFlgRmvBin = horzcat(raNtFlgRmvBin,ix);
			end
		end

		raNtBin_mid(raNtFlgRmvBin) = [];

		raNtBin_midFlip = [raNtBin_mid,fliplr(raNtBin_mid)];
% 		raNtSpread = [log10(raNtBinMin), fliplr(log10(raNtBinMax))];
		raNtSpread = [log10(raNtBin25p), fliplr(log10(raNtBin75p))];

		f3 = fill(raNtBin_midFlip,raNtSpread,'b','FaceAlpha', 0.5);

		if strcmp(cntrLine,'Mean')
			m3 = plot(raNtBin_mid,log10(raNtBinMean),'Color',[27/255 27/255 100/255],'LineWidth', 3);
		elseif strcmp(cntrLine,'Median')
			m3 = plot(raNtBin_mid,log10(raNtBinMedian),'Color',[27/255 27/255 100/255],'LineWidth', 3);
		end
	end
	
	
	title(['Nt vs. Temp - ' cntrLine ' and Spread - CIP - All PS/LS Spirals (Group 2)']);

	set(gca,'XMinorTick','on','YMinorTick','on','XDir','reverse');
% 	if ~isempty(tempRangeAll)
% 		xlim(tempRangeAll);
% 	end
	if ~isempty(NtSprdLim)
		ylim(NtSprdLim)
	end
	ylabel('log_{10}(N_t) [cm^{-3}]')
	xlabel(sprintf('Temp (%cC)', char(176)));
	set(findall(gcf,'-property','FontSize'),'FontSize',28);
	view([90 -90])
	
	legend([f1 f2 f3],{'Transition Zone','Stratiform Region','Rear Anvil'},'Location','eastoutside');
	
	if saveFigs
		set(gcf,'Units','Inches');
		pos = get(gcf,'Position');
		set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
		print([saveDir '/CIP-NtTempSprd_' num2str(avgTime) 's/CIP_Nt-Temp_' cntrLine '-Spread-Fill_' num2str(avgTime) 's_Grp2'],Ftype,Fres)
	end
end