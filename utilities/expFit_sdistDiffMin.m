%% Minimizing chisquare difference of NDobs and NDfit
function [nl, psdFit, chisquare] = expFit_sdistDiffMin(bin_mid, psdObs)

	opts = optimset('Display','off');
	lb = [-Inf,-Inf];
	ub = [Inf, Inf];
	
	% General initial params to start with
	x00 = [0.3 16];
	
	
	% Get initial best guess for n0/lambda
	% Find where psdObs>0 just for least square fitting
	I = find(psdObs>0);
	[nl, ~] = lsqcurvefit( @(nl, bin_mid)f_myExp(nl, bin_mid),...
		x00, bin_mid(I), psdObs(I),lb,ub,opts);

	% Exponential moment fitting based on the previous initial guess
	
	f_fit = @(nl)f_sum_chisquare(nl, bin_mid, psdObs);
	
% 	x00 = nl;
	[nl, chisquare, flag, output] = fminsearch(f_fit, x00);
	psdFit = f_myExp(nl,bin_mid);
end

function [ chisquare ] = f_sum_chisquare( nl, bin_mid, psdObs)
	% Calculate the chisquare value of distribution
	
	if any(isnan(f_myExp(nl, bin_mid))) || any(isinf(f_myExp(nl, bin_mid)))
		chisquare = Inf;
		return;
	end
	
	psdFit = f_myExp(nl, bin_mid);
	
	chisquare = sum( (psdObs-psdFit).^2 );
	
end


function [ psdFit ] = f_myExp( nl, x )
	% Calculate the continuous psd for the distribution
	
	psdFit = nl(1).*exp(-nl(2).*x);
	
end