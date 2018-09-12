%% transmitter main
[HEX,CenterFreq,SymbolRate,FEC]=TRANSMITER_START();
[MsgSize,~]=size(HEX);
HEX=[HEX;repmat('20',18*ceil(MsgSize/18)-MsgSize,1)];
NumberOfBurtsts=ceil(MsgSize/18);
NumBurst=1;  % Number Of Bursts;
ModStCnt=[]; % Modulator State Control (Modulator PhaseOffset Should be zero at the beginning of each Burst)
while NumBurst<=NumberOfBurtsts
    Burst=HEX((18*(NumBurst-1)+1):18*NumBurst,:); % Each burst contains 18 Bytes
    CodedData=MessageCoder(Burst,FEC);
    ModStCnt=mod(2*sum(CodedData)-length(CodedData),4)/2;
    BasebandBits(:,NumBurst)=[ModStCnt;1;0;0;1;0;1;0;CodedData;0;1;0;0;1;0;0;1]; % padding zeros at two sides of bursts
    NumBurst=NumBurst+1;
end
%% Channel estimate train signal
ESS=CHIRP(CenterFreq,SymbolRate);
Train=[1,1,1,0,0,0,1,0,1,0,1,0,0,0,1,0,0,1,1,1,1,0,1,1,0,1,0,0,1,...
    1,1,0,1,0,0,0,0,1,0,1,0,1,0,1,1,1,0,0,0,1,0,1,1,1,0,0,0,1,0,1,1,1,0,0]';
AllBits=[];
for j=1:NumberOfBurtsts
    if mod(j,2)
        AllBits=[AllBits;Train;BasebandBits(:,j)];
    else
        AllBits=[AllBits;BasebandBits(:,j)];
    end
end
%% modulate BaseBandBits
SigBase=gmsk_mod(double(AllBits),SymbolRate);
%% Adding train data for channel estimation
t=(0:1:(length(SigBase)-1))'./44100;
MsgTx=SigBase.*exp(1i*(2*pi*CenterFreq*t));
PhisycalSignal=[zeros(22050,1);ESS;real(MsgTx);zeros(22050,1)];
audiowrite('5K_200s_3_4.wav',PhisycalSignal,44100);