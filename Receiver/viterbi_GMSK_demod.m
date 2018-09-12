function DemodRes=viterbi_GMSK_demod(Sig_LP,SymbolRate)
% Viterbi_GMSK_demodulator
SamplePerSymbol=floor(44100./SymbolRate);
persistent trellis
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
%% loading Low pass signal
symbol=1;
L=length(Sig_LP)/SamplePerSymbol;
index=trellis.Index;
metrics_ind=index-16*(mod((1:16)',2)-1);
SurvivalPath=zeros(16,L);
R=reshape(Sig_LP,SamplePerSymbol,L);
R=conj(R);
metricI=trellis.MetricFunc.I*R;
metricII=trellis.MetricFunc.II*R;
metrics_ALL=metricII(:,[1;1]*(1:size(metricII,2)));
metrics_ALL(:,1:2:end) = metricI;
metrics_ALL=real(metrics_ALL);
while(symbol<=L)
        metrics=metrics_ALL(:,(symbol*2-1):symbol*2);
        ind=trellis.PathMeter(index)+metrics(metrics_ind);
        ind=2-(ind(:,1)>ind(:,2));
        SurvivalPath(:,symbol)=index((ind-1)*16+(1:16)'); % SurvivalPath
        trellis.PathMeter=trellis.PathMeter(SurvivalPath(:,symbol))+metrics(SurvivalPath(:,symbol)-16*(mod((1:16)',2)-1));
    symbol=symbol+1; 
end
[~,MAX_METER]=max(trellis.PathMeter);
DemodRes=zeros(L,1);
DemodRes(end)=MAX_METER;
for i=(L-1):-1:1
    DemodRes(i)=SurvivalPath(DemodRes(i+1),i+1);
end
DemodRes=mod(DemodRes,2);
    