%% recover Baseband signal
function [TrainLow]=BasebandRecovery(WinData,CenterFreq,curser,BaseGuard)
t=((curser+1-BaseGuard:(curser+length(WinData)-BaseGuard))/44100)';
H=hilbert(WinData,10*length(WinData));
H=resample(H,1,10);
TrainLow=(WinData+1i*H).*exp(-2*pi*1i*CenterFreq*t);
FirBase = fir1(64,0.1);
TrainLow=2*conv(TrainLow,FirBase,'same');
TrainLow=TrainLow(BaseGuard+1:end-BaseGuard);
