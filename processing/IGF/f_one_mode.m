function [nml, mobs, mfit, chisquare, psdfit, flag, output] = f_one_mode(psd, bin_div, moments)
	% IGF fitting routine.
	% Created by Shichu 10/29/2016
	% Input:
	%   psd -> particle size distribution (preferrably in 1/L/um)
	%   bin_div -> EDGES of each size bin. Length should be length(psd)+1
	%       i.e. [bin1_min, bin1_max, bin2_max, ... , binLast_min, binLast_max]
	%       unit in um
	%   moments -> specify at least 3 moments to fit to. eg. [0 2 3]
	% Output:
	%   nml -> [ log10(N0), mu, lambda ]
	%       unit N0 in 1/L/um, mu unitless, lambda 1/um
	%       Fitted psd(D) = N0 x D^mu x exp(-lambda D)
	%   mobs -> the moments of the input psd
	%   mfit -> the moments of the fitted gamma psd
	%   chisquare -> best fit (nml) chisquare, usually very close to 0
	%   psdfit -> fitted psd at the same points of the original input psd
	%   flag -> the exit signal of matlab fit function
	%   output -> detailed fitting output from matlab fit function
	%   For flag and output, see matlab fminsearch function docu for details
	%   https://www.mathworks.com/help/matlab/ref/fminsearch.html
	% Last edited by Shichu 10/29/2016
	
	% To modify this for the exponential fit, we only need N0 and lambda params
	% Will also need to modify the my_gamma function to be the actual exponential 
	
	%% Parameter settings
	MaxFunEvals = 6e4;
	MaxIter = 6e4;
	TolFun = 1e-12;
	TolX = 1e-5;
	%TolFun=1e-3;
	%TolX=1e-3;
	%%% A general initial paras to start with
	x00 = [log10(300) -1 0.0014];
	
	bin_diff = diff(bin_div);   % bin width
	bin_mid = (bin_div(1:end-1)+bin_div(2:end))/2;  % bin mid size
	
	mobs = zeros(size(moments));
	mfit = zeros(size(moments));
	
	for szmoment = 1:size(moments,2)
		mobs(szmoment) = sum( psd.*bin_diff.*bin_mid.^moments(szmoment) );
	end
	
	option1 = optimset('MaxFunEvals',MaxFunEvals, 'MaxIter',MaxIter);
	option2 = optimset('MaxFunEvals',MaxFunEvals, 'MaxIter',MaxIter, 'TolFun',TolFun, 'TolX',TolX);
	
	%%% upper/lower bounds of parameters nml
	lb = [-Inf, -Inf, -Inf];
	ub = [Inf, Inf, Inf];
	
	%%% Find where psd>0 just for least square fitting of log(psd)
	I = find(psd>0);
	[nml, ~] = lsqcurvefit( @(nml, bin_mid)log10(f_mygamma(nml, bin_mid)),...
		x00, bin_mid(I), log10(psd(I)), lb, ub, option1 );
	
	%%% IGF moment fitting based on the previous initial guess
	x00 = nml;
	f_fit = @(nml)f_sum_chisquare(nml, moments, bin_mid, bin_diff, mobs);
	% default tolerences are both 1e-4
	[nml, chisquare, flag, output] = fminsearch(f_fit, x00, option2);
	psdfit = f_mygamma(nml,bin_mid);
	for szmoment=1:size(mfit,2)
		mfit(szmoment) = sum(psdfit.*bin_diff.*bin_mid.^(moments(szmoment)));
	end
	
end


function [ chisquare ] = f_sum_chisquare( nml, moments, bin_mid, bin_diff, mobs)
	% Calculate the chisquare value of distribution
	%   nml: log10(N0), mu, lambda
	
	if any(isnan(f_mygamma(nml, bin_mid))) || any(isinf(f_mygamma(nml, bin_mid)))
		chisquare = Inf;
		return;
	end
	
	mfit = zeros(size(mobs));
	
	for szmoment = 1:size(mfit,2)
		mfit(szmoment) = sum( f_mygamma(nml, bin_mid).*bin_diff.*bin_mid.^(moments(szmoment)) );
		%    int_f=@(x) f_mygamma(nml,x).*x.^(moments(szmoment));
		%    mfit(szmoment)=integral(int_f,xmin,xmax);
	end
	
	chisquare = sum( (mfit-mobs).^2./mobs./mfit );
	
end


function [ psd ] = f_mygamma( nml, x )
	% Calculate the continuous psd for the distribution
	%   nml: log10(N0), mu, lambda
	
	%N1=10^nml(1);
	%m1=nml(2);
	%l1=nml(3);
	
	psd = 10^nml(1).*x.^nml(2).*exp(-nml(3).*x);
	
end
