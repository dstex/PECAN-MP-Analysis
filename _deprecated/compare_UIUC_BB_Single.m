clearvars;
readBobBlack_SDs;

tempUIUC_cipConc = uiucCIPSDall.conc_minR_all(21244:21249,:);
tempUIUC_pipConc = uiucPIPSDall.conc_minR_all(21244:21249,:);
tempBB_cipConc = bb_cip_conc_centerIn_all(1615,:);
tempBB_pipConc = bb_pip_conc_centerIn_all(1170,:);

tempUIUC_cipSmplVol = uiucCIPSDall.sampleVol_all(21244:21249,:); % cm3
tempUIUC_pipSmplVol = uiucPIPSDall.sampleVol_all(21244:21249,:);
tempBB_cipSmplVol = (bb_cip_smplVol_L_all(1615))*1000; %cm3
tempBB_pipSmplVol = (bb_pip_smplVol_L_all(1170))*1000;

tempUIUC_cipNt = nansum(tempUIUC_cipConc'.*repmat(uiuc_cip_binwidth/10,[1 size(tempUIUC_cipConc,1)]),1); % Units of cm-3
tempUIUC_pipNt = nansum(tempUIUC_pipConc'.*repmat(uiuc_pip_binwidth/10,[1 size(tempUIUC_pipConc,1)]),1);
tempBB_cipNt = nansum(tempBB_cipConc*0.025)/1000; % Units of cm-3 (Converted from L-1 0.025 mm-1)
tempBB_pipNt = nansum(tempBB_pipConc*0.1)/1000; % Units of cm-3 (Converted from L-1 0.1 mm-1)

% Brown and Francis 1995 coefficients
a = 0.00294; % g cm-3
b = 1.9;

tempUIUC_cipBF = a*(uiuc_cip_binMid/10).^b; % Units of grams(?)
tempUIUC_pipBF = a*(uiuc_pip_binMid/10).^b;
tempBB_cipBF = a*(bb_cip_binMid/10).^b;
tempBB_pipBF = a*(bb_pip_binMid/10).^b;


for ix=1:size(tempUIUC_cipConc,1)
	tempUIUC_cip_massIce(ix,:) = ((tempUIUC_cipConc(ix,:).*(uiuc_cip_binwidth/10)')*1e6).*tempUIUC_cipBF'; % Units of g m-3
	tempUIUC_pip_massIce(ix,:) = ((tempUIUC_pipConc(ix,:).*(uiuc_pip_binwidth/10)')*1e6).*tempUIUC_pipBF';
end

% Following gives same result as above steps, but without loop
% tempUIUC_cip_massIce = ((tempUIUC_cipConc'.*repmat(uiuc_cip_binwidth/10,[1 size(tempUIUC_cipConc,1)]))*1e6).*repmat(a*(uiuc_cip_binMid/10).^b,[1 size(tempUIUC_cipConc,1)]);
% tempUIUC_pip_massIce = ((tempUIUC_pipConc'.*repmat(uiuc_pip_binwidth/10,[1 size(tempUIUC_pipConc,1)]))*1e6).*repmat(a*(uiuc_pip_binMid/10).^b,[1 size(tempUIUC_pipConc,1)]);

tempBB_cip_massIce = ((tempBB_cipConc*0.025)*1000).*tempBB_cipBF; % Units of g m-3 (Converted to m-3 from L-1 0.025 mm-1)
tempBB_pip_massIce = ((tempBB_pipConc*0.1)*1000).*tempBB_pipBF; % Units of g m-3 (Converted to m-3 from L-1 0.1 mm-1)

tempUIUC_cipIWC = nansum(tempUIUC_cip_massIce,2);
tempUIUC_pipIWC = nansum(tempUIUC_pip_massIce,2);
tempBB_cipIWC = nansum(tempBB_cip_massIce);
tempBB_pipIWC = nansum(tempBB_pip_massIce);