clear;
clc;
clear DfeEqualizer;
clear viterbi_GMSK_demod
% Reciever Main
%% GUI of reciever for Setiings
[CenterFreq,SymbolRate,FEC,path,DemodMod,nfwdweights,nfbkweights,cons1,cons2]=RECIEVER_GUI();
%% initial variables for other functions
ESS_sig=Chirp(CenterFreq,SymbolRate);
%% loading signal and set variables
Audio=audioread(path);
Audio=Audio(:,1);
if DemodMod==2;Audio=Audio./max(Audio);end
SamplePerSymbol=floor(44100/SymbolRate);
BaseGuard=SamplePerSymbol;
%% Filter data
Filter=fir1(128,[CenterFreq-0.75*SymbolRate,CenterFreq+0.75*SymbolRate]*2/44100);
Audio=conv(Audio,Filter,'same');
%% Syncronization Part1: find begining of symbols.
[acor,lag] = xcorr(Audio,ESS_sig);
[~,I] = max(acor);
lagDiff = lag(I);
%% Creating Train modulated signal for using in Equalizer
Train=[1,1,1,0,0,0,1,0,1,0,1,0,0,0,1,0,0,1,1,1,1,0,1,1,0,1,0,0,1,...
    1,1,0,1,0,0,0,0,1,0,1,0,1,0,1,1,1,0,0,0,1,0,1,1,1,0,0,0,1,0,1,1,1,0,0]';
TrainSec=gmsk_mod(Train,SymbolRate);
curser=lagDiff+8820;
%% Extract data Frames
TrainCheck=0;
TXT=[];
SymbolsInFrame=64;
while curser+(SymbolsInFrame*SamplePerSymbol)<=length(Audio)
    if mod(TrainCheck,3)==0
        TrainFlag=1;
        TrainMod=TrainSec;
        SymbolsInFrame=64;
    else
        TrainFlag=0;
        TrainMod=[];
        switch FEC
           case '1/2'
              SymbolsInFrame=466;
           case '3/4'
              SymbolsInFrame=316;
           case '5/6'
              SymbolsInFrame=286;
        end
    end
    WinData=Audio((curser+1-BaseGuard):(curser+(SymbolsInFrame*SamplePerSymbol)+BaseGuard));
    [FrameBase]=BasebandRecovery(WinData,CenterFreq,curser,BaseGuard);
    switch DemodMod
        case 1
            [DemodRes,MSE]=DfeEqualizer(FrameBase,...
               SymbolRate,nfwdweights,nfbkweights,TrainMod,cons1,cons2);
        case 2
            [DemodRes,MSE]=DfeEqualizer(FrameBase,...
               SymbolRate,nfwdweights,nfbkweights,TrainMod,cons1,cons2);
        case 3
            if TrainFlag==1  % PHASE OF Carrier
            CompensatePhi=TrainMod.'*conj(FrameBase);
            CompensatePhi=CompensatePhi./norm(CompensatePhi);
            end
            FrameBase=FrameBase*CompensatePhi;
            DemodRes=viterbi_GMSK_demod(FrameBase,SymbolRate);
    end
    %% reliability of frame
    if TrainFlag==0 && sum(DemodRes([2:8,(SymbolsInFrame-7):SymbolsInFrame])==[1;0;0;1;0;1;0;0;1;0;0;1;0;0;1])>= 10
        DecData=MessageDecoder(DemodRes(9:end-8),FEC);
        TEXTDEC=[128 64 32 16 8 4 2 1]*double(DecData');
        TXT=[TXT,char(TEXTDEC)];
        disp(['Frame # ',num2str(TrainCheck-floor(TrainCheck/3)),' Recieved successfully']);
    elseif TrainFlag==1
        QOS=(sum(DemodRes(3:62)==Train(3:62))/60)*100; % should be near 100.  50 means Equalizer doesn't work properly
        disp(['Quality of Signal in Train # ' num2str(1+floor(TrainCheck/3)),' :',num2str(QOS)]);
    else
        disp('Failed to decode Frame')
    end
    curser=curser+(SymbolsInFrame*SamplePerSymbol);
    TrainCheck=TrainCheck+1;
end
TXT=strtrim(TXT)
            