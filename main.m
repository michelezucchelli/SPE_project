%% 
% Start
clc;
clear;
close all;
format long;
disp('Simulation started:');
disp(datetime("now"));
tic

%% 
% Read config file

try    
    config = load_config('config_file');
catch
    error('Specified config file does not exists or exists with bad value(s)');
end

techniques = textscan(config{'type'},'%d','Delimiter',',');
techniques = permute( techniques{1}, [2,1] );
T = str2num(config{'num_Tx'});
R = str2num(config{'num_Rx'});
snr_min = str2num(config{'snr_min'});
snr_max = str2num(config{'snr_max'});
points = str2num(config{'points'});
M = str2num(config{'M'});
N_tot = str2num(config{'num_symb'});
seeds = str2num(config{'seeds'});
confidence_level = str2num(config{'confidence'});
x_lim = str2num(config{'xlim'});
%% 
% Init variables

% set random number generator seed to obtain reproducible results
num_loop = numel(seeds);

% Number of symbols that each tx antenna transmit
N = N_tot / T;

% Define SNR vector
SNR_dB = linspace(snr_min, snr_max, points);
lenSNR = length(SNR_dB);

% Init 3D matrix for received symbols corrupted by noise. These are lenSNR
% matrix of dimension RxN
received_symbols_noisy = zeros(lenSNR, R, N);

% define BER null vector for the number of techniques used
bers = zeros(length(find(techniques==1)), 1, lenSNR);

% Creeate a folder for the current configuration where to store all i/o files
dirname = strcat(config{'num_Tx'}, '_', config{'num_Rx'}, '_', ...
            config{'M'}, '_', config{'num_symb'}, '_', ...
            config{'snr_min'}, '_', config{'snr_max'},  '_', ...
            config{'points'});
% dirname = strcat(config{'num_Tx'}, '_', config{'num_Rx'}, '_', ...
%         config{'points'},  '_', config{'M'}, '_', config{'num_symb'});
[status, msg, msgID] = mkdir(dirname);

% Check if the folder is correctly created
if status ~= 1
    error("Error! Folder not created")
end

% define the filename for the BER files
filenameZF = strcat(dirname, filesep,'ber_zf.txt');
filenameMM = strcat(dirname, filesep,'ber_mmse.txt');
filenameZVB = strcat(dirname, filesep,'ber_zf_vblast.txt');
filenameMVB = strcat(dirname, filesep,'ber_mmse_vblast.txt');

%% 
% % Real program begins here

for m=1:num_loop
    % set random number generator seed to obtain reproducible results
    rng(seeds(m));

    % Generate random symbols to transmit
    tx_symbols = randi([0 M-1], T, N);

    % Modulate the simbols to transmit
    modulated_symbols = qammod(tx_symbols, M, 'UnitAveragePower', true);
    % generate the channel matrix and save it on file
    H = (randn(R, T) + 1i * randn(R, T)) / sqrt(2); % MIMO channel mat. (H)

    % compute the theoretical received symbols (no noise) and save them on
    % files
    received_symbols = H * modulated_symbols;
    % % end

    % ---------------- Received Noisy ---------------- % 
    for f=1:lenSNR
        SNR = 10^(SNR_dB(f) / 10);
        noise_variance = 1/SNR;
        noise = sqrt(noise_variance/2) * (randn(R, N) + 1i * randn(R, N));
        received_symbols_noisy(f,:,:) = received_symbols + noise;
    end
    % ------------------------------------------------ % 

    % ---------------- TECHNIQUES ---------------- %
    if techniques(1) % Zero Forcing
        for f=1:lenSNR
            % received symbols noisy in 2 dimension (current SNR value)
            rsn_2d(:,:) = received_symbols_noisy(f,:,:);
            % perform zf
            zf_symbols = zf(H,rsn_2d);
            % demodulate
            demodulated_symbols_zf = qamdemod(zf_symbols, M, 'UnitAveragePower', true);
            % compute SER and convert to BER
            [~, ser] = symerr(tx_symbols, demodulated_symbols_zf);
            bers(1,f) = ser / log2(M);
        end
        try
            write_config(filenameZF, seeds(m), bers(1,:));
        catch ME
            disp('Not writing on file since BER file for this technique already exists ')
        end
    end


    if techniques(2) %MMSE
        for f=1:lenSNR
            SNR = 10^(SNR_dB(f) / 10);
            % received symbols noisy in 2 dimension (current SNR value)
            rsn_2d(:,:) = received_symbols_noisy(f,:,:);
            % perform mmse
            mmse_symbols = mmse(T,H,rsn_2d,SNR);
            % demodulate
            demodulated_symbols_mmse = qamdemod(mmse_symbols, M, 'UnitAveragePower', true);
            % compute SER and convert to BER
            [~, ser] = symerr(tx_symbols, demodulated_symbols_mmse);
            bers(2,f) = ser / log2(M);
        end
        try
            write_config(filenameMM, seeds(m), bers(2,:));
        catch ME
            disp('Not writing on file since BER file for this technique already exists ')
        end
    end

    if techniques(3) % VBLAST
        for f=1:lenSNR
            % received symbols noisy in 2 dimension (current SNR value)
            rsn_2d(:,:) = received_symbols_noisy(f,:,:);
            % perform vblast
            vblast_symbols = zfvblast(T,R,M,N,H,rsn_2d); 
            % contrary to zf and mmse, the demodulation operation
            % is already applied within the vblast function
            % compute SER and convert to BER
            [~, ser] = symerr(tx_symbols, vblast_symbols);
            bers(3,f) = ser / log2(M);
        end
        try
            write_config(filenameZVB, seeds(m), bers(3,:));
        catch ME
            disp('Not writing on file since BER file for this technique already exists ')
        end
    end

    if techniques(4) % MMSE-VBLAST
        for f=1:lenSNR
            SNR = 10^(SNR_dB(f) / 10);
            % received symbols noisy in 2 dimension (current SNR value)
            rsn_2d(:,:) = received_symbols_noisy(f,:,:);
            % perform vblast
            mmse_vblast_symbols = mmse_vblast(T,R,M,N,H,rsn_2d,SNR); 
            % contrary to zf and mmse, the demodulation operation
            % is already applied within the vblast function
            % compute SER and convert to BER
            [~, ser] = symerr(tx_symbols, mmse_vblast_symbols);
            bers(4,f) = ser / log2(M);
        end
        try
            write_config(filenameMVB, seeds(m), bers(4,:));
        catch ME
            disp('Not writing on file since BER file for this technique already exists ')
        end
    end
end

fclose('all');

%% Read bers data from file
% The reason why we read the BERs from file is because one could run a
% simulation with ten seed, and then another simulation with ten different
% seeds, and plot the results considering all 20 seeds instead of just the
% ten of the current simulation run

% Read Zero forcing data
if isfile(filenameZF) && techniques(1)
    mat_ber_zf = load_config(filenameZF);
    ber_zf = zeros(numEntries(mat_ber_zf), points);
    fclose('all');
    for k=1:numEntries(mat_ber_zf)
        ber_zf(k,:) = str2num(mat_ber_zf{seeds(k)});
    end
end

% Read MMSE data
if isfile(filenameMM) && techniques(2)
    mat_ber_mmse = load_config(filenameMM);
    ber_mmse = zeros(numEntries(mat_ber_mmse), points);
    fclose('all');
    for k=1:numEntries(mat_ber_mmse)
        ber_mmse(k,:) = str2num(mat_ber_mmse{seeds(k)});
    end
end

% Read vblast data
if isfile(filenameZVB) && techniques(3)
    mat_ber_vblast = load_config(filenameZVB);
    ber_vblast = zeros(numEntries(mat_ber_vblast), points);
    fclose('all');
    for k=1:numEntries(mat_ber_vblast)
        ber_vblast(k,:) = str2num(mat_ber_vblast{seeds(k)});
    end
end

% Read mmse-vblast data
if isfile(filenameMVB) && techniques(4)
    mat_ber_mmse_vblast = load_config(filenameMVB);
    ber_mmse_vblast = zeros(numEntries(mat_ber_mmse_vblast), points);
    fclose('all');
    for k=1:numEntries(mat_ber_mmse_vblast)
        ber_mmse_vblast(k,:) = str2num(mat_ber_mmse_vblast{seeds(k)});
    end
end
%% Compute CI and Plot

% Legend array and index
leg = strings([1,length(techniques)*2]);
leg_index = 1;

% fig = figure('Position', get(0, 'Screensize'));
fig = figure;
% ---- ZF ---- %
if exist("ber_zf", "var")
    mean_zf = mean(ber_zf);
    snr = linspace(snr_min, snr_max, length(mean_zf));
    % Plotta mean values with error bar for the CI
    semilogy(snr, mean_zf, 'm');
    leg(1,leg_index) = "ZF";
    leg_index = leg_index + 1;
    hold on;    

    % Compute CI - You can choose between bootsrap or asymptotic method,
    % at least 30 seeds are needed for the asymptotic method
    % % Compute CI using bootstrap method
    % for c=1:points
    %     ci_zf(c,:) = bootstrap_ci(ber_zf(:,c), @mean, confidence_level, 999);
    % end
    % errorbar(snr, mean_zf, mean_zf - ci_zf(:, 1)', ci_zf(:, 2)' - mean_zf, 'm');
    % leg(1,leg_index) = "C.I. ZF";
    % leg_index = leg_index + 1;
    % hold on;
    % % End bootstrap method
    % % Compute CI using asymptotic method
    if length(ber_zf) > 30
        ci_zf = zeros(points, 2);
        for c=1:points
            ci_zf(c,:) = asymptotic_ci(ber_zf(:,c), confidence_level);
        end 
        if ~isempty(ci_zf)
            errorbar(snr, mean_zf, mean_zf - ci_zf(:, 1)', ci_zf(:, 2)' - mean_zf, 'm');
            leg(1,leg_index) = "C.I. ZF";
            leg_index = leg_index + 1;
            hold on;
        end
    end    
    % % End asymptotic method
    
end

% ---- MMSE ---- %
if exist("ber_mmse", "var")
    mean_mmse = mean(ber_mmse);

    snr = linspace(snr_min, snr_max, length(mean_mmse));
    % Plotta mean values with error bar for the CI
    semilogy(snr, mean_mmse, 'b--');
    leg(1,leg_index) = "MMSE";
    leg_index = leg_index + 1;
    hold on;    
    
    % Compute CI - You can choose between bootsrap or asymptotic method,
    % at least 30 seeds are needed for the asymptotic method
    % % Compute CI using bootstrap method
    % for c=1:points
    %     ci_mmse(c,:) = bootstrap_ci(ber_mmse(:,c), @mean, confidence_level, 999);
    % end
    % errorbar(snr, mean_mmse, mean_mmse - ci_mmse(:, 1)', ci_mmse(:, 2)' - mean_mmse, 'b--');
    % leg(1,leg_index) = "C.I. MMSE";
    % leg_index = leg_index + 1;
    % hold on;
    % % End bootstrap method
    % % Compute CI using asymptotic method
    if length(ber_mmse) > 30
        ci_mmse = zeros(points, 2);
        for c=1:points
            ci_mmse(c,:) = asymptotic_ci(ber_mmse(:,c), confidence_level);
        end 
        if ~isempty(ci_mmse)
            errorbar(snr, mean_mmse, mean_mmse - ci_mmse(:, 1)', ci_mmse(:, 2)' - mean_mmse, 'b--');
            leg(1,leg_index) = "C.I. MMSE";
            leg_index = leg_index + 1;
            hold on;
        end
    end    
    % % End asymptotic method
end

% ---- VBLAST ---- %
if exist("ber_vblast", "var")
    % Compute CI
    mean_vblast = mean(ber_vblast);


    snr = linspace(snr_min, snr_max, length(mean_vblast));
    % Plotta mean values with error bar for the CI
    semilogy(snr, mean_vblast, 'r');
    leg(1,leg_index) = "ZF-VBLAST";
    leg_index = leg_index + 1;
    hold on;

    % Compute CI - You can choose between bootsrap or asymptotic method,
    % at least 30 seeds are needed for the asymptotic method
    % % Compute CI using bootstrap method
    % for c=1:points
    %     ci_vblast(c,:) = bootstrap_ci(ber_vblast(:,c), @mean, confidence_level, 999);
    % end
    % errorbar(snr, mean_vblast, mean_vblast - ci_vblast(:, 1)', ci_vblast(:, 2)' - mean_vblast, 'r');
    % leg(1,leg_index) = "C.I. ZF-VBLAST";
    % leg_index = leg_index + 1;
    % hold on;
    % % End bootstrap method
    % % Compute CI using asymptotic method
    if length(ber_vblast) > 30
        ci_vblast = zeros(points, 2);
        for c=1:points
            ci_vblast(c,:) = asymptotic_ci(ber_vblast(:,c), confidence_level);
        end 
        if ~isempty(ci_vblast)
            errorbar(snr, mean_vblast, mean_vblast - ci_vblast(:, 1)', ci_vblast(:, 2)' - mean_vblast, 'r');
            leg(1,leg_index) = "C.I. ZF-VBLAST";
            leg_index = leg_index + 1;
            hold on;
        end
    end    
    % % End asymptotic method    
end

% ---- MMSE-VBLAST ---- %
if exist("ber_mmse_vblast", "var")
    % Compute CI
    mean_mmse_vblast = mean(ber_mmse_vblast);

    snr = linspace(snr_min, snr_max, length(mean_mmse_vblast));
    % Plotta mean values with error bar for the CI
    semilogy(snr, mean_mmse_vblast, 'g--');
    leg(1,leg_index) = "MMSE-VBLAST";
    leg_index = leg_index + 1;
    hold on;
    
    % Compute CI - You can choose between bootsrap or asymptotic method,
    % at least 30 seeds are needed for the asymptotic method
    % % Compute CI using bootstrap method
    % for c=1:points
    %     ci_mmse_vblast(c,:) = bootstrap_ci(ber_mmse_vblast(:,c), @mean, confidence_level, 999);
    % end
    % errorbar(snr, mean_mmse_vblast, mean_mmse_vblast - ci_mmse_vblast(:, 1)', ci_mmse_vblast(:, 2)' - mean_mmse_vblast, 'g--');
    % leg(1,leg_index) = "C.I. MMSE-VBLAST";
    % leg_index = leg_index + 1;
    % hold on;
    % % End bootstrap method
    % % Compute CI using asymptotic method
    if length(ber_mmse_vblast) > 30
        ci_mmse_vblast = zeros(points, 2);
        for c=1:points
            ci_mmse_vblast(c,:) = asymptotic_ci(ber_mmse_vblast(:,c), confidence_level);
        end 
        if ~isempty(ci_mmse_vblast)
            errorbar(snr, mean_mmse_vblast, mean_mmse_vblast - ci_mmse_vblast(:, 1)', ci_mmse_vblast(:, 2)' - mean_mmse_vblast, 'g--');
            leg(1,leg_index) = "C.I. MMSE-VBLAST";
            leg_index = leg_index + 1;
            hold on;
        end
    end    
    % % End asymptotic method
end

% ZF          --> MAGENTA
% MMSE        --> BLUE
% VBLAST      --> RED
% MMSE-VBLAST --> GREEN
hold off;
legend(leg)
title({ strcat('BER for: ', num2str(T), 'x' , num2str(R), ' Antennas - M=', num2str(M))})
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
xlim(x_lim)
grid on
saveas(fig, strcat(config{'num_Tx'}, '_', config{'num_Rx'}, '_', ...
        config{'M'}, '_', config{'num_symb'}, '_', config{'snr_min'}, ...
        '_', config{'snr_max'}, '_', config{'points'}, '.jpg'));

%% 
% Plot distribution of data using QQ plot
% % Plot zf if exists
% if exist("ber_zf", "var")
%     figure;
%     qqplot(ber_zf(:,1));    
% end
% 
% % Plot mmse if exists
% if exist("ber_mmse", "var")    
%     figure;
%     qqplot(ber_mmse(:,1));
% end
% 
% % Plot vblast if exists
% if exist("ber_vblast", "var")
%     figure;
%     qqplot(ber_vblast(:,1));    
% end
% 
% % Plot vblast if exists
% if exist("ber_mmse_vblast", "var")    
%     figure;
%     qqplot(ber_mmse_vblast(:,1));
% end

%% 
% Close all files and stop

fclose('all');
disp('Simulation ended:');
disp(datetime("now"));
disp('Simulation time: ')
disp(toc)
beep