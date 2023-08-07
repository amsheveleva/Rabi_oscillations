%% Rabi Oscillations and Rabi Frequency

clear, clc, clf, close all


% % Define the folder path 
 folder = 'the folder path';
% 
% % Get a list of files with the extension .DSC and T1
 filelist = dir(fullfile(folder, '*rabi*.DSC'));

output = fullfile(folder,"results"); % sets the output folder
nnn=length(filelist);
for i = 1:nnn
    
    [time{i}, int{i}, par{i}, filename{i}] = eprload(fullfile(folder,filelist(i).name));
    
    % Takes the magnetic field and attenuations to check against the title.
    
    magnetic_field(i) = (par{i}.B0VL*1e3);
    attenuations{i} = par{i}.PowerAtten;

    a(i) = str2num(attenuations{i}(1:end-3)); % transforms the string into double
    
    B1(i) = sqrt((10^(-0.1*a(i)))/(10^(-1.6))); % calculates relative B1
    
    % Takes real part of data, normalises and offsets.
    
    int{i} = real(int{i});
    int{i} = int{i}/max(int{i});
    int{i} = int{i} + 2*i;
      
    % Fourier transform
    
    Fs = par{i}.XPTS;                  % Sampling frequency is the number of points
    T = 1/Fs;                          % Sampling period
    L{i} = length(int{i});             % Length of signal
    
    Y{i} = fft(int{i});        % Compute the Fourier transform of the signal.
    
    P2{i} = abs(Y{i}/L{i});     % Compute the two-sided spectrum P2.
    P1{i} = P2{i}(1:L{i}/2+1);          % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
    P1{i}(2:end-1) = 2*P1{i}(2:end-1);
    P1{i}(1:4) = 0; % removes the intense signal at x = 0
    
    P1{i} = P1{i}/max(P1{i}); % normalises by the maximim
    P1{i} = P1{i} + 1.5*i; % offsets the signal
    
    f{i} = Fs*(0:(L{i}/2))/L{i}; %Define the frequency domain f and plot the single-sided amplitude spectrum P1.
    
    
    %     Export graphs
    rabi_freq(i) = f{i}(P1{i} == max(P1{i}));
    rabi_freq(16) = 6;
    rabi_freq(20) = 3;
    rabi_freq(1) = 33;
    rabi_freq(2) = 32;
    rabi_freq(3) = 30;
    rabi_freq(4) = 30;

    headings = ["Nutation_pulse_length_ns", "Echo_Intensity_arb_unit"];
    graph = table(time{i}, int{i}, 'VariableNames', headings);
    
    if isfolder(fullfile(output,string(magnetic_field(i)))) == 1
        writetable(graph, fullfile(output,string(magnetic_field(i)),"Int_" + filelist(i).name + ".txt"))
    else
        mkdir(fullfile(output,string(magnetic_field(i))))
        writetable(graph, fullfile(output,string(magnetic_field(i)),"Int_" + filelist(i).name + ".txt"))
    end
    
    headings = ["Rabi_frequency_MHz", "FT_arb_unit"];
    graph = table(f{i}', P1{i}, 'VariableNames', headings);
    writetable(graph, fullfile(output,string(magnetic_field(i)), "FT_" + filelist(i).name + ".txt"))
    
    % Save rabi_freq(i) and B1(i) to separate txt files
    headings_freq = ["Rabi_frequency_MHz", "B1"];
    data_freq = table(rabi_freq(i), B1(i), 'VariableNames', headings_freq);
    writetable(data_freq, fullfile(output, string(magnetic_field(i)), "FreqData_" + filelist(i).name + ".txt"))
    
    % Finding maxima
    
   
    
end

