function[INVD] = invdist(TM)
%% CODE
N = size(TM,2);
B = (TM-eye(N));
B(1:N,N) = ones(N,1);
%
o = zeros(1,N);
o(N) = 1;
INVD = o*inv(B); % Invariant Distribution
%