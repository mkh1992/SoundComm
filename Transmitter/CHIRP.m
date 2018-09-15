function ESS_sig=CHIRP(CenterFrequency,SymbolRate)
%% ESS signal and inverse
BW=SymbolRate;
duration = 0.2;
f1=CenterFrequency-BW/2;
f2=CenterFrequency+BW/2;
n=(0:(1/44100):(duration-1/44100)).'; % designed for sampling rate 44100
k=exp(n*log(f2/f1)/duration);
ESS_sig=sin(2*pi*f1*duration/log(f2/f1)*(k-1));
% ESS_sig_inv=ESS_sig(end:-1:1)./k;