function ESS_sig=Chirp(CenterFreq,SymbolRate)

%% ESS signal Creating ...
BW=SymbolRate;
duration = 2;
f1=CenterFreq-BW/2;
f2=CenterFreq+BW/2;
n=(0:(1/44100):(duration-1/44100)).'; % designed for sampling rate 44100
k=exp(n*log(f2/f1)/duration);
ESS_sig=sin(2*pi*f1*duration/log(f2/f1)*(k-1));