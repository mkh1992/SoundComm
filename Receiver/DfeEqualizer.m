% DFE-RLS Equalizer With Viterbi Decision Device for GMSK Modulation
% This is a DFE Equalizer with symbolwise delays
% This Equaliser uses RLS algorithm for setting Weights. also the 
% Deciesion devise in this Equalizer is a Viterbi Trellis search

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HELP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function inputs:
%UnEqSig: Un Equalized Signal
%trellis: for Decision Device
%SymbolRate: Transmitter SymbolRate
%nfwdweights: Number of Forward Registers
%nfbkweights: number of backward Registers
%Weights: Initial Weights for ALL registers Weights=[ForwardWeights;BachwardWeights]
%TrainMod: Modulated Training signal, importing this signal to Function means that we Are in Training mod so FeedBackMod will be set to 1
%CovMatrix: Covariance matrics of RLS Algorithm importing this matrix to function means that you want use RLS Method
function [DemodRes,MSE]=DfeEqualizer(InputSignal,SymbolRate,...
    nfwdweights,nfbkweights,TrainMod,lambda,Constant2)
%% internal settings
AdaptiveStep=1;
SamplePerSymbol=floor(44100./SymbolRate);
NumNewSymbols=length(InputSignal)/SamplePerSymbol;
e=zeros(SamplePerSymbol/AdaptiveStep,1);
%% set equalizer in correct modes
if isempty(TrainMod) && isempty(Constant2)
    FeedBackMod=2;
    AdaptiveMod=2;
    mu=lambda;
end
if isempty(TrainMod) && ~isempty(Constant2)
    FeedBackMod=2;
    AdaptiveMod=1;
end 
if ~isempty(TrainMod) && isempty(Constant2)
    FeedBackMod=1;
    AdaptiveMod=2;
    mu=lambda;
    TrainPointer=1:SamplePerSymbol;
end
if ~isempty(TrainMod) && ~isempty(Constant2)
    FeedBackMod=1;
    AdaptiveMod=1;
    TrainPointer=1:SamplePerSymbol;
end
%% persistant variavles
persistent trellis Weights FwReg BwReg CovMatrix symbol UnEqSig SurvivalPath
if isempty(trellis)
        trellis=struct('states',[],'path',[],'MetricFunc',[],'PathMeter',[],'SurvivalPath',[],'Index',5);
    trellis.states=[0,1,1;0,1,-1;0,-1,1;0,-1,-1
        pi/2,1,1;pi/2,1,-1;pi/2,-1,1;pi/2,-1,-1
        pi,1,1;pi,1,-1;pi,-1,1;pi,-1,-1
     3*pi/2,1,1;3*pi/2,1,-1;3*pi/2,-1,1;3*pi/2,-1,-1];
    trellis.path=[5,6;7,8;13,14;15,16;9,10
        11,12;1,2;3,4;13,14;15,16;5,6
        7,8;1,2;3,4;9,10;11,12]; % first coloumn is for input 1 and second for input -1
    % calculate trellis metric function
    g=@(t) 1/2*(qfunc(2*pi*0.3*(t-1/2)/sqrt(log(2)))-qfunc(2*pi*0.3*(t+1/2)/sqrt(log(2))));
    q=@(t) integral(g,-2,t);
    % piece wise Q
    t=linspace(0,1,SamplePerSymbol+1);
    t=t(1:end-1);
    t=[t,1+t,1*2+t]-1.5*1;
    q_piece=zeros(length(t),1);
    for i=1:length(t)
        q_piece(i)=q(t(i));
    end
    % cunstruct metric_functions
    trellis.Index=zeros(16,2);
    trellis.PathMeter=zeros(16,1);
    for i=1:16
        hist_phase=repelem(trellis.states(i,1),SamplePerSymbol)';
        bpre=trellis.states(i,2).*q_piece((2*SamplePerSymbol+1):end);
        pre=trellis.states(i,3).*q_piece((SamplePerSymbol+1):(2*SamplePerSymbol));
        Iself=1*q_piece(1:SamplePerSymbol);
        IIself=-1*q_piece(1:SamplePerSymbol);
        trellis.MetricFunc.I(i,:)=exp(1i*(hist_phase+(bpre+pre+Iself)*pi));
        trellis.MetricFunc.II(i,:)=exp(1i*(hist_phase+(bpre+pre+IIself)*pi));
        trellis.Index(i,:)=find(trellis.path==(ceil(i/2)*2-1))';
    end
end
if isempty(Weights)
    Weights=[1;zeros(nfwdweights+nfbkweights-1,1)];  % initial State for Weights
end
if isempty(FwReg)
    FwReg=zeros(SamplePerSymbol,nfwdweights);
    BwReg=zeros(SamplePerSymbol,nfbkweights);
end
if isempty(CovMatrix)
    CovMatrix=Constant2*eye(length(Weights));
end
if isempty(symbol)
   symbol=1;
end
if isempty(UnEqSig)
    UnEqSig=InputSignal;
else
    UnEqSig=[UnEqSig;InputSignal];
end
L=length(UnEqSig)/SamplePerSymbol;
if isempty(SurvivalPath)
   SurvivalPath=zeros(16,L);
else
    SurvivalPath=[SurvivalPath,zeros(16,NumNewSymbols)];
end
while symbol<=L
    %% loading data to Fwd registers And Shift others
    Pointer=((symbol-1)*SamplePerSymbol +1):(symbol*SamplePerSymbol );
    FwReg=[UnEqSig(Pointer),FwReg(:,1:end-1)];
    %% Decision Device-Viterbi
    index=trellis.Index;
    metrics_ind=index-16*(mod((1:16)',2)-1);
    R=conj([FwReg,BwReg]*conj(Weights)); % Equalizer output to Decision Device
    metrics=real([trellis.MetricFunc.I*R,trellis.MetricFunc.II*R]);
    ind=trellis.PathMeter(index)+metrics(metrics_ind);
    ind=2-(ind(:,1)>ind(:,2));
    SurvivalPath(:,symbol)=index((ind-1)*16+(1:16)'); % SurvivalPath
    trellis.PathMeter=trellis.PathMeter(SurvivalPath(:,symbol))+metrics(SurvivalPath(:,symbol)-16*(mod((1:16)',2)-1));
    %% desired signal (this part has 2 mods Training and Decision FeedBack)
    switch FeedBackMod
        case 1
            %Training
            d=TrainMod(TrainPointer);
            TrainPointer=TrainPointer+SamplePerSymbol;
        case 2
            % Decision FeedBack
            [~,FinState]=max(trellis.PathMeter);
            InitState=SurvivalPath(FinState,symbol);
            if mod(FinState,2)==1
                d=trellis.MetricFunc.I(InitState,:).';
            else
                d=trellis.MetricFunc.II(InitState,:).';
            end
    end
    %% updating Weights RLS and LMS algorithms (this part has 2 mods)
    switch AdaptiveMod
        case 1  %With RLS Algorithm
            for i=1:AdaptiveStep:(SamplePerSymbol)
                    u=[FwReg(i,:).';BwReg(i,:).'];
                    k=lambda^(-1)*CovMatrix*u./(1+lambda^(-1)*u'*CovMatrix*u);
                    e(i)=d(i)-Weights'*u;
                    Weights=Weights+k*e(i)';
                    CovMatrix=lambda^(-1)*(eye(length(Weights))-k*u')*CovMatrix;
            end
        case 2  %With LMS Algorithm
            for i=1:AdaptiveStep:(SamplePerSymbol)
                 u=[FwReg(i,:).';BwReg(i,:).'];
                 e(i)=d(i)-Weights'*u;
                 Weights=Weights+mu*e(i)'*u;
            end
    end
    %% loading data to Bwd registers
    switch FeedBackMod
        case 1
            BwReg=[d,BwReg(:,1:end-1)];
        case 2
            BKRegUpdate=zeros(min(nfbkweights,symbol)+1,1);
            BwReg(:,length(BKRegUpdate):end)=BwReg(:,length(BKRegUpdate)-1:end-1);
            [~,BKRegUpdate(end)]=max(trellis.PathMeter);
            for i=(length(BKRegUpdate)-1):-1:1
                BKRegUpdate(i)=SurvivalPath(BKRegUpdate(i+1),symbol+i+1-length(BKRegUpdate));
                if mod(BKRegUpdate(i+1),2)==1
                   BwReg(:,length(BKRegUpdate)-i)=trellis.MetricFunc.I(BKRegUpdate(i),:).';
                else
                   BwReg(:,length(BKRegUpdate)-i)=trellis.MetricFunc.I(BKRegUpdate(i),:).';
                end
            end
    end
    symbol=symbol+1;
end
[~,MAX_METER]=max(trellis.PathMeter);
DemodRes=zeros(NumNewSymbols,1);
DemodRes(end)=MAX_METER;
for i=(L-1):-1:(L-NumNewSymbols+1)
    DemodRes(i-L+NumNewSymbols)=SurvivalPath(DemodRes(i+1-L+NumNewSymbols),i+1);
end
DemodRes=mod(DemodRes,2);
MSE=abs(e);