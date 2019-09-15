function x=blocksignal(t,tperiod,tduration)
tper=round(t./tperiod).*tperiod;
x= (t>tper)&(t<tper+tduration);
