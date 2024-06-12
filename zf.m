function zf_symbols = zf(H, received_symbols_noisy)
    % check that the matrices dimensions are ok for multiplication
    if size(pinv(H),2) == size(received_symbols_noisy,1)
        % apply zf
        zf_symbols = pinv(H) * received_symbols_noisy;
    else
        error('[ERROR] matrices size are not adequate for multiplication')
    end
end

%% Function Description
% Compute the Zero Forcing technique on a matrix of received symbols.
% Input: 
%   - H = channel matrix. Dimension: RxT, where R is the number of
%   receive antennas ant T the number of transmit antennas
%   - received_symbols_noisy = matrix of received symbols corrupted by
%   noise. Dimension: RxN, where N is the number of symbols transmitted by
%   every T antennas
% Output:
%   - zf_symbols = modified version of the received symbols where zf 
%   technique is applied. Dimensionion: TxN