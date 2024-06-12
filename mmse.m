function mmse_symbols = mmse(T, H, received_symbols_noisy, SNR)
    % check that the matrices dimensions are ok for matrix multiplication
    if size(pinv(H),2) == size(received_symbols_noisy,1)
        % apply mmse
        W = inv(H' * H + (eye(T) / SNR)) * H';
        mmse_symbols = W * received_symbols_noisy;
    else
        error('[ERROR] matrices size are not adequate for multiplication')
    end
end

%% Function Description
% Compute the Minimum Meand Squared Error technique on a matrix of received 
% symbols.
% Input: 
%   - T = number of transmit antennas
%   - H = channel matrix. Dimension: RxT, where R is the number of
%   receive antennas
%   - received_symbols_noisy = matrix of received symbols corrupted by
%   noise. Dimension: RxN, where N is the number of symbols transmitted by
%   each T antennas
%   - SNR = Signal to Noise Ratio value, computed as 10^(SNR_dB / 10)
% Output:
%   - mmse_symbols = modified version of the received symbols where mmse 
%   technique is applied. Dimensionion: TxN