%% Rabi Oscillations and Rabi Frequency

clear, clc, clf, close all

% Initialize variables to store data from all txt files
allData = table();  % Initialize an empty table

% Define the folder where the txt files are located
folder = 'C:\Users\mbdssase\Dropbox (The University of Manchester)\Azzah\VO_Zr_Hf_project\VO_Hf_sample\VO_Hf_Qband_rabi\results\1191.4';

% Get a list of all txt files in the folder
filelist = dir(fullfile(folder, '*FreqData*.txt'));

% Loop through each file and read the data
for i = 1:numel(filelist)
    % Read the data from the current file
    data = readtable(fullfile(folder, filelist(i).name));
    magnetic_field = regexp(filelist(1).name, '(?<=\D)(\d+)(?=G)', 'match'); %Temperature
    magnetic_field = str2double( magnetic_field)/10;
    % Append the data to the allData table
    allData = [allData; data];
end

% Extract columns from the merged table to MATLAB workspace
rabi_freq = allData.Rabi_frequency_MHz;
B1 = allData.B1;
figure();

subplot(2,2,4)
plot(B1, rabi_freq, 'o', 'color', 'k')
hold on
xlabel('Relative B1 (arb. unit)'); xlim([min(B1)-0.5 max(B1)+0.5]); xticks([0:1:10])
ylabel('\Omega_R (MHz)');
ylim([0, rabi_freq(1)+1.5]);
xlim([-0.1, B1(1)+0.5]);


% Linear fit
% 
P = polyfit(B1,rabi_freq,1); % fits a polynomial of order 1
yfit = polyval(P,B1); % evaluates the polynomial
plot(B1,yfit, 'red')


% Exporting Figure

filename = sprintf('B1_rabifreq_Hf.txt');
dlmwrite(filename, [B1,rabi_freq], 'delimiter', '\t', 'precision', 6);

%%%%  plotting Intensityvs time traces



% Define the folder where the txt files are located
% Get a list of all txt files in the folder
filelist2 = dir(fullfile(folder, '*Int*.txt'));

for i = 1:numel(filelist2)
    % Read the data from the current file
    data2 = readtable(fullfile(folder, filelist2(i).name));
    
    % Extract time and intensity data from the table
    time = data2.Nutation_pulse_length_ns;
    intensity = data2.Echo_Intensity_arb_unit;
    intensity=(intensity-min(intensity))/max(intensity);
    subplot(2,2,[1,3])
    % Plot the data for the current sample
    plot(time, intensity+i*0.5);
    
    hold on; % To keep all plots in the same figure
    title("Rabi Oscillations at " + magnetic_field + " mT", 'Interpreter', 'none');
    xlabel('Nutation pulse length (ns)');
    ylabel('Normalised echo intensity (arb. unit offset)'); yticks([]);
    set(gcf, 'Position', [100, 100, 1000, 800]);  % Adjust figure size if needed
        
end

filelist3 = dir(fullfile(folder, '*FT*.txt'));

for i = 1:numel(filelist3)
    % Read the data from the current file
    data3 = readtable(fullfile(folder, filelist3(i).name));
    
    % Extract time and intensity data from the table
    Freq = data3.Rabi_frequency_MHz;
    Rabi_FT= data3.FT_arb_unit;
    Rabi_FT=(Rabi_FT-min(Rabi_FT))/max(Rabi_FT);
    subplot(2,2,2)
    % Plot the data for the current sample
    plot(Freq, Rabi_FT+i*0.5);
    
    hold on; % To keep all plots in the same figure
     title("Rabi Frequency", 'Interpreter', 'none');
     xlabel('\Omega_R (MHz)');
    ylabel('Fourier Transform (arb. unit offset)'); yticks([]);
     xlim([-0.1, 200]);   
end

figFilename = sprintf('Field_%dM_Rabi_Hf.fig',magnetic_field);
    savefig(figFilename);

% Save the figure in .png format
    pngFilename = sprintf('Field_%dM_Rabi_Hf.png',magnetic_field);
    saveas(gcf, pngFilename, 'png'); 




%  rabi_freq(1) = 33;
%  rabi_freq(2) = 32;
%  rabi_freq(3) = 29;
%  rabi_freq(4) = 33;
%  rabi_freq(16) = 6;
%  rabi_freq(20) = 3;




% Plotting Rabi
%     colorMap = parula(nnn);
%    
%     xrange = [0 600]; xlim([min(xrange) max(xrange)]); xticks(min(xrange):100:max(xrange))
%     %colours = {[0 0 0], [1 0 0], [0 .8 0], [0 0 1], [1 0.5 0], [1 0 1], [0 0.5 1], [.5 0 .5], [0.5 0.5 0]}; % define the colours of the plots
%     
%     plot(time{i}, int{i}, 'LineWidth', 1.5, 'color', colorMap(i, :))
%     text(max(xrange), mean(int{i}), attenuations(i),'color', colorMap(i, :), 'FontSize',12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom') %writes the legend next to the curves
%     
%     hold on
%     
%     title("Rabi Oscillations at " + magnetic_field(i) + " mT", 'Interpreter', 'none');
%     xlabel('Nutation pulse length (ns)');
%     ylabel('Normalised echo intensity (arb. unit offset)'); yticks([]);
%     set(gcf, 'Position', [-1000 100 600 600]);
% 
%      subplot(2,2,2)
%     xrange = [0 60]; xlim([min(xrange) max(xrange)]); xticks([min(xrange):10:max(xrange)]);
%     
%     plot(f{i},P1{i}, 'color', colorMap(i, :), 'LineWidth', 1.5)
%     %     text(max(xrange), mean(P1{i}), attenuations(i),'color', colours{i}, 'FontSize',12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom') %writes the legend next to the curves
%     hold on
%     
%     title("Rabi Frequency", 'Interpreter', 'none');
%     xlabel('\Omega_R (MHz)');
%     ylabel('Fourier Transform (arb. unit offset)'); yticks([]);




% subplot(2,2,4)
% plot(B1, rabi_freq, 'o', 'color', 'k')
% hold on
% xlabel('Relative B1 (arb. unit)'); xlim([min(B1)-0.5 max(B1)+0.5]); xticks([0:1:10])
% ylabel('\Omega_R (MHz)');
% 
% % Linear fit
% % 
% P = polyfit(B1,rabi_freq,1); % fits a polynomial of order 1
% yfit = polyval(P,B1); % evaluates the polynomial
% plot(B1,yfit, 'red')
% 
% % Exporting Figure
% 
% filename = sprintf('B1_rabifreq_Hf.txt');
% dlmwrite(filename, [B1',rabi_freq'], 'delimiter', '\t', 'precision', 6);
% 
% saveas(gcf, fullfile(output,string(magnetic_field(i)),"Rabi Oscillations at " + magnetic_field(i) + " mT.png"))
% saveas(gcf, fullfile(output,string(magnetic_field(i)),"Rabi Oscillations at " + magnetic_field(i) + " mT.fig"))
% writematrix(P,'omega_b1_param_Hf.xlsx')