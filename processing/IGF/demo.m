clear;
load data.mat

[nml, mobs, mfit, chisquare, psdfit, flag, output] = f_one_mode(PSD, bin_div, [0 2 3]);

stairs(bin_mid, PSD,'r');
hold on
stairs(bin_mid, psdfit,'b');

legend('Observed','Fitted');
xlim([bin_div(1),bin_div(end)]);
ylim([min(PSD(PSD~=0)),max(PSD)]);
xlabel('D (um)');
ylabel('N(D) 1/um/L');
set(gca,'yscale','log','xscale','log');
grid minor
hold off