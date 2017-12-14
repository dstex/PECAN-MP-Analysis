clearvars; close all;


savePath = '/Users/danstechman/GoogleDrive/School/Research/PECAN/Microphysics/plots/';

flights = {'20150617','20150620','20150701','20150702','20150706','20150709'};

% Exponential Fit Params
sse_exp = [];
rsquare_exp = [];
adjrsquare_exp = [];
rmse_exp = [];
n0_exp = [];
lmda_exp = [];

% IGF params
chisquare_igf = [];
n0_igf = [];
mu_igf = [];
lmda_igf = [];

for iflt = 1:length(flights)
	flight = flights{iflt};
	
	mFiles = dir([savePath flight '/CIP-Fitting/' flight '_FitComparison-CIP_Whole*.mat']);
	if length(mFiles) > 1
		fprintf('More than one fit comparison m-file in directory');
		return
	end
	dataTmp = load([mFiles.folder '/' mFiles.name]);
	
	sprlNames = dataTmp.sprlNames;
	
	for ix = 1:length(sprlNames)
		gof = dataTmp.gof_expLMWhl.(sprlNames{ix});
		sse_exp = [sse_exp; gof.sse];
		rsquare_exp = [rsquare_exp; gof.rsquare];
		adjrsquare_exp = [adjrsquare_exp; gof.adjrsquare];
		rmse_exp = [rmse_exp; gof.rmse];
		nl_exp = dataTmp.cip_expLM_nlWhl.(sprlNames{ix});
		n0_exp = [n0_exp; nl_exp(1)];
		lmda_exp = [lmda_exp; nl_exp(2)];
		
		chisquare_igf = [chisquare_igf; dataTmp.cip_igf_chisquareWhl.(sprlNames{ix})];
		nml_igf = dataTmp.cip_igf_nmlWhl.(sprlNames{ix});
		n0_igf = [n0_igf; nml_igf(1)];
		mu_igf = [mu_igf; nml_igf(2)];
		lmda_igf = [lmda_igf; nml_igf(3)];
	end
	
end
clearvars('-except','savePath','flights','sse_exp','rsquare_exp','adjrsquare_exp','rmse_exp','n0_exp','lmda_exp',...
	'chisquare_igf','nml_igf','n0_igf','mu_igf','lmda_igf');