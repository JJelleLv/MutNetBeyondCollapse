function [periods,frequencies]=fourier_spect(times,X)
tend=times(end);
periods=1./((0:length(times)-1)/tend);
frequencies=abs(fft(X));

