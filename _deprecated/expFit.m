function [nl, mobs, mfit, chisquare, psdfit, flag, output] = expFit(psd, bin_div, moments)
	% Exponential fitting routine.
	% Adapted by Dan Stechman 12/20/2016 from IGF routine by Shichu (10/29/2016)
	% Input:
	%   psd -> particle size distribution (preferrably in 1/L/um)
	%   bin_div -> EDGES of each size bin. Length should be length(psd)+1
	%       i.e. [bin1_min, bin1_max, bin2_max, ... , binLast_min, binLast_max]
	%       unit in um
	%   moments -> specify at least 2 moments to fit to. eg. [0 2]
	% Output:
	%   nl -> [ log10(N0),lambda ]
	%       unit N0 in 1/L/um, lambda 1/um
	%       Fitted psd(D) = N0 x exp(-lambda D)
	%   mobs -> the moments of the input psd
	%   mfit -> the moments of the fitted exponential psd
	%   chisquare -> best fit (nl) chisquare, usually very close to 0
	%   psdfit -> fitted psd at the same points of the original input psd
	%   flag -> the exit signal of matlab fit function
	%   output -> detailed fitting output from matlab fit function
	%   For flag and output, see matlab fminsearch function docu for details
	%   https://www.mathworks.com/help/matlab/ref/fminsearch.html
	% Last edited by Dan Stechman 12/20/2016
	
	%% Parameter settings
	MaxFunEvals = 6e4;
	MaxIter = 6e4;
	TolFun = 1e-12;
	TolX = 1e-5;
	%TolFun=1e-3;
	%TolX=1e-3;
	%%% A general initial paras to start with
	x00 = [log10(300) 0.0014];
	
	bin_diff = diff(bin_div);   % bin width
	bin_mid = (bin_div(1:end-1)+bin_div(2:end))/2;  % bin mid size
	
	mobs = zeros(size(moments));
	mfit = zeros(size(moments));
	
	for szmoment = 1:size(moments,2)
		mobs(szmoment) = sum( psd.*bin_diff.*bin_mid.^moments(szmoment) );
	end
	
	option1 = optimset('MaxFunEvals',MaxFunEvals, 'MaxIter',MaxIter,'Display','off');
	option2 = optimset('MaxFunEvals',MaxFunEvals, 'MaxIter',MaxIter, 'TolFun',TolFun, 'TolX',TolX);
	
	%%% upper/lower bounds of parameters nl
	lb = [-Inf, -Inf];
	ub = [Inf, Inf];
	
	%%% Find where psd>0 just for least square fitting of log(psd)
	I = find(psd>0);
	[nl, ~] = lsqcurvefit( @(nl, bin_mid)log10(f_myExp(nl, bin_mid)),...
		x00, bin_mid(I), log10(psd(I)), lb, ub, option1 );

	%%% Exponential moment fitting based on the previous initial guess
	x00 = nl;
	f_fit = @(nl)f_sum_chisquare(nl, moments, bin_mid, bin_diff, mobs);
	% default tolerences are both 1e-4
	[nl, chisquare, flag, output] = fminsearch(f_fit, x00, option2);
	psdfit = f_myExp(nl,bin_mid);
	for szmoment=1:size(mfit,2)
		mfit(szmoment) = sum(psdfit.*bin_diff.*bin_mid.^(moments(szmoment)));
	end
	
end


function [ chisquare ] = f_sum_chisquare( nl, moments, bin_mid, bin_diff, mobs)
	% Calculate the chisquare value of distribution
	%   nl: log10(N0), lambda
	
	if any(isnan(f_myExp(nl, bin_mid))) || any(isinf(f_myExp(nl, bin_mid)))
		chisquare = Inf;
		return;
	end
	
	mfit = zeros(size(mobs));
	
	for szmoment = 1:size(mfit,2)
		mfit(szmoment) = sum( f_myExp(nl, bin_mid).*bin_diff.*bin_mid.^(moments(szmoment)) );
		%    int_f=@(x) f_myExp(nl,x).*x.^(moments(szmoment));
		%    mfit(szmoment)=integral(int_f,xmin,xmax);
	end
	
	chisquare = sum( (mfit-mobs).^2./mobs./mfit );
	
end


function [ psd ] = f_myExp( nl, x )
	% Calculate the continuous psd for the distribution
	%   nl: log10(N0), lambda
	
	%N1=10^nl(1);
	%l1=nl(2);
	
	psd = 10^nl(1).*exp(-nl(2).*x);
	
end
