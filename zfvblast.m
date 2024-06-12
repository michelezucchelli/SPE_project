function a = zfvblast(T, R, M, N, H, r)    
    % init variables
    Hv = H;
    k = zeros(1,T);
    y = zeros(T,N);
    a = zeros(T,N);
    w = zeros(T,R);

    % counter used to not select the same k(i) more than once
    r_counter = zeros(1,T); 
    
    for i=1:T
        % equation (9b) and (9h)
        G = pinv(Hv);
        % ------------ 

        % equation (9c) and (9i)
        norms_squared = sum(G.^2, 2);
        norms_squared = norms_squared + reshape(r_counter, T, 1);
        [~, k(i)] = min(norms_squared);
        % don't select it again in the next iterations
        r_counter(k(i)) = 100; 
        % ------------
        
        % equation (9d)
        w(k(i),:) = G(k(i),:);
        % usato per confermare delta
        % disp(w(k(i),:)*Hv(:,k(i))) % should be equal to 1 + 0i        
        % ------------

        % equation (9e)
        y(k(i),:) = w(k(i),:) * r;
        % ------------

        % equation (9f)
        a(k(i),:) = qamdemod(y(k(i),:), M, 'UnitAveragePower', true);
        % ------------

        % equation (9g)
        temp = qammod(a(k(i),:),M, 'UnitAveragePower', true) .* H(:,k(i));
        % temp = H(:,k(i)) * qammod(a(k(i),:),M); % uguale a sopra
        r = r - temp;
        % ------------
        
        % equation (9h)
        Hv(:,k(i)) = 0; 
        % ------------        
    end
end

%% Function Description
% Reference: V-BLAST: An Architecture for Realizing Very High Data Rates
% Over the Rich-Scattering Wireless Channel
% P. W. Wolniansky, G. J. Foschini, G. D. Golden, R. A. Valenzuela
%
% Compute the V-BLAST technique on a matrix of received symbols.
% Input: 
%   - T = number of transmit antennas
%   - R = number of receive antennas
%   - M = modulation order
%   - N = number of symbols transmitted by each T antennas
%   - H = channel matrix. Dimension: RxT
%   - r = matrix of received symbols corrupted by noise. Dimension: RxN
% Output:
%   - a = modified version of the received symbols where vblast 
%   technique is applied. Dimensionion: TxN
% 
% Dimensions:
%  G: TxR
%  w(k(i)): 1xR
%  y(k(i),:): 1xN
%  a(k(i),:): 1xN
