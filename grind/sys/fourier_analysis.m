%program that read data and make Fourier transform and red noise significance spectra of the data
%data are in sigle files with only one column which is species
%concentration ( without zeros, detrended , normalized)

function spectra = fourier_analysis(X, method, ~)
if nargin < 2
   method = '';
end
%read data
%X=load(filename); %data analysis
N = length(X); % number of data points

% Sampling frequency
fs = 1;

%Calculate power spectra of data - Elisa's method
% !! fft() calculates with unit frequency, on interval [0,1/2]
Y = fft(X, N);  % make Fourier transform
Pyy = Y .*  conj(Y) / N;
normPyy = Pyy / (var(X));
% transform from unit frequency to real frequency
freq = (0:floor(N / 50)) / N * fs;
time = (freq.^(-1));
time = transpose(time);
freq = transpose(freq);
spectra = Y;

plot(freq, normPyy(1:length(freq)));
if 1
   %Calculate autocorrelation coefficient data
   alpha = corrcoef(X(1:N - 1), X(2:N));
   % Calculate red noise spectra
   P=cell(1,floor(N / 2) + 1);
   for k = 1:(floor(N / 2) + 1)
      P{k} = (1 - (alpha(1, 2)^2)) / (1 + alpha(1, 2)^2 - 2 * alpha(1, 2) * cos(2 * pi * (k - 1) / N));
   end
   %impose a level of significance
   Chi22 =  5.99;    % 95th percentile of Chisquare distribution with 2 degree of freedom
   SignP=cell(1,floor(N / 2) + 1);
   for k = 1:(floor(N / 2) + 1)
      SignP{k} = (P{k} .* Chi22) / 2;
   end
   SignP = transpose(SignP);
   P = transpose(P);
   
   
   % Calculate Periodogram
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   switch method
    case 'ham'
      %Simple Hamming window does not do a good job
      w1 = window('hamming', N);
      [Pss, ff] = periodogram(X, w1, N);
      lg = 'estimator: hamming window';
      
    case 'wel'
      %Welch: 8 windows overlapping by 50%
      [Pss, ff] = pwelch(X, [], [], N);
      lg = 'estimator: Welch - 8 overlapping windows';
      
    case 'mtm'
      %MTM: p data taper
      [Pss, ff] = pmtm(X, (param + 1) / 2, N);
      lg = 'estimator: Multitaper method - 5 taper';
      
    case 'yul'
      %Autocorr: par= AR coefficients (p=2 is identical to red noise, test it!)
      [Pss, ff] = pyulear(X, param, N);
      lg = 'estimator: autocorrelation (Y-W) - p=50';
      
    case 'bur'
      if (2 * param > N)
         param = 2; %only red spectrum in case of too view data
      end
      
      [Pss, ff] = pburg(X, param, N);   % this method yields slightly higher peaks
      lg = 'e';
      %stimator: autocorrelation (Burg) - p = 50';
      
    case 'cov'
      [Pss, ff] = pmcov(X, param, N);
      lg = 'estimator: autocorrelation (covar) - p=50';
      
    otherwise
      % raw periodogram - identical to fft-periodogram
      % !! Periodogram() using unit frequency, on interval [0,pi] !!
      [Pss, ff] = periodogram(X, [], N);
      lg = 'estimator: raw periodogram';
   end
   
   disp(lg);
   % transform frequencies from output interval [0,pi] to [0, fs/2]
   ff  = ff * fs / (2 * pi);
   tt  = (ff).^(-1);
   % transform PSD from unit frequencies to fs
   % Multiply with PI: has sometrhing to do with unit-range [0,pi]
   % Tested with nino-sst dataset and Figure 3 of Torrence and Compo, 1998
   Pss = Pss / var(X) * pi;
   
   %Return spectra to main call
   spectra = [ff, tt, Pss, P, SignP];
   
   %time plot
   subplot(2, 1, 1);
   semilogx(tt,Pss,'k',tt,P,'g',tt,SignP,'r');
   xlabel('Time');
   ylabel('Power');
   xlim([5 250]);
   ylim([0 25]);
   
   subplot(2, 1, 2);
   semilogx(time,normPyy(1:floor(N/2)+1),'k',time,P,'g',time,SignP,'r');
   xlabel('Time');
   ylabel('Power');
   xlim([5 250]);
   ylim([0 25]);
   legend('raw periodogram');
end
