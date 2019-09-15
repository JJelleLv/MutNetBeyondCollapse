function [PDFvalue,CDFvalue]=func_VectAngleProb(Angle,n)

%% PDF taken from: 
%% Cai, T., Fan, J., & Jiang, T. (2013). Distributions of angles in random packing on spheres. The Journal of Machine Learning Research, 14(1), 1837-1864.

if isnan(Angle)==1
    PDFvalue=NaN;
    CDFvalue=NaN;
    return
end

%% convert degrees to radians
RADIAN=(Angle/180).*pi;

%% Probability Distribution Function
PDFvalue=(1/(pi.^0.5)).*(gamma(n./2))/(gamma((n-1)./2)).*sin(RADIAN).^(n-2);

%% determine cumulative distribution from 0 to radian
NRsteps=1e3;
NULLtoRADIAN=[RADIAN/NRsteps:RADIAN/NRsteps:RADIAN]';
for STEPNR=1:NRsteps
    testRADIAN=NULLtoRADIAN(STEPNR,1);
    PDF_NULLtoRADIAN(STEPNR,1)=dot((1/(pi.^0.5)).*(gamma(n./2))/(gamma((n-1)./2)),sin(testRADIAN).^(n-2));
end

%% cumulative distribution function
CDF=cumtrapz(NULLtoRADIAN,PDF_NULLtoRADIAN);

CDFvalue=CDF(NRsteps,1);
