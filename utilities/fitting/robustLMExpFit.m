%% cftool generated fitting algorithm 
function [nl, obsFit, fitresult, gof] = robustLMExpFit(binMid, obsPSD)
	%CREATEFIT(CIP_BINMID,PSD3)
	%  Create a fit.
	%
	%  Data for 'ExpFit' fit:
	%      X Input : cip_binMid
	%      Y Output: PSD3
	%  Output:
	%      fitresult : a fit object representing the fit.
	%      gof : structure with goodness-of fit info.
	
	
	%% Fit: 'ExpFit'.
	[xData, yData] = prepareCurveData( binMid, obsPSD );
	
	% Set up fittype and options.
% 	ft = fittype( 'a*exp(-b*x)', 'independent', 'x', 'dependent', 'y' );
% 	opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% 	opts.Display = 'Off';
% 	opts.Robust = 'LAR';
	ft = fittype( 'a*exp(-b*x)', 'independent', 'x', 'dependent', 'y' );
	opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
	opts.Algorithm = 'Levenberg-Marquardt';
	opts.Display = 'Off';
	opts.Robust = 'Bisquare';
	opts.StartPoint = [0.278498218867048 0.546881519204984];
	
	% Fit model to data.
	[fitresult, gof] = fit( xData, yData, ft, opts );
	
	nl = coeffvalues(fitresult);
	
	obsFit = f_myexp(nl,xData);
	
end

function [ obsFit ] = f_myexp( nl, x )
	% Calculate the continuous psd for the distribution
	
	obsFit = nl(1).*exp(-nl(2).*x);
	
end
	
