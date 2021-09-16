function [iniIFfit] = IF_param(IF,orderIF,iii)
% Intrinsic Chirp Component Decomposition(ICCD)
% the code is only for complex-valued data analysis
%%%%%%%%%%%%%%%%%%%%%%%  input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sig£ºmeasured signal,a row vector
% SampFreq: sampling frequency
% iniIFset: the instantaneous frequency series
% orderIF: the order of the Fourier model used for fitting the instantaneous frequency of the signal
% orderamp£ºFourier order for characterizing signal amplitudes
% alpha£ºTikhonov regularization parameter for ICCD.
% iii Samples that have been correctly sampled
%%%%%%%%%%%%%%%%%%%%%%%  output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extr_Sig: the reconstructed components
% ampmatrix: the estimated amplitudes
%
%%%%%%%%%%%%%%%%%%%%Construction of the Fourier matrix%%%%%%%%%%%%%%%%%%%%%%%
SampFreq=64;
N=length(IF);
f0 = SampFreq/2/N;%%Base Frequency
l_t = 2*orderIF + 1;%%the number of the Fourier coefficients
tmatrix = zeros(N,l_t);
tmatrix(:,1) = ones(N,1);
dt = [0:N-1]/SampFreq;
for j = 2:l_t
        tmatrix(:,j) = cos(2*pi*f0*(j-1)*dt); 
        if j >(l_t+1)/2
            tmatrix(:,j) = sin(2*pi*f0*(j-((l_t+1)/2))*dt);
        end
end

kmatrix =  tmatrix;
kmatrix1=kmatrix(iii,:);

%% regularized least-squares solution

size(IF(iii));
size(kmatrix1);
theta =( kmatrix1'*kmatrix1)\(kmatrix1'*IF(iii));%coefficient vector
iniIFfit=kmatrix*theta;

end
    


