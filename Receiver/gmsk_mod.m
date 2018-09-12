function modulated_data=gmsk_mod(message,SymbolRate,phase_offset,InitialBits)
%% settings
SamplePerSymbol=floor(44100/SymbolRate);
switch nargin
    case 1        
        SamplePerSymbol=8;
        InitialBits=[0;1]; % first element is initial phase offset
        phase_offset=0;
    case 2
        InitialBits=[0;1]; % first element is initial phase offset
        phase_offset=0;
    case 3
        InitialBits=[0;1]; % first element is initial phase offset
end
T=1;
BT=0.3;
Fd=0.25; %pick frequency deviation
%% G(t)
g=@(t) 1/(2*T)*(qfunc(2*pi*BT*(t/T-1/2)/sqrt(log(2)))-qfunc(2*pi*BT*(t/T+1/2)/sqrt(log(2))));
q=@(t) integral(g,-2,t);
%% piece wise Q
t=linspace(0,T,SamplePerSymbol+1);
t=t(1:end-1);
t=[t,T+t,T*2+t]-1.5*T;
q_piece=zeros(length(t),1);
for i=1:length(t)
    q_piece(i)=q(t(i));
end
%% message 
I=[InitialBits;message]; % appending zeros first of data
I=2*I-1;
I=[0;I];
hist_phase=cumsum(I(1:(end-3)));
%% self
self_filt=repmat(q_piece(1:SamplePerSymbol),length(message),1);
self_data=repelem(I(4:end),SamplePerSymbol);
modulated_phase=self_filt.*self_data*4*pi*T*Fd;
%% previuse symbol
pre_filt=repmat(q_piece((SamplePerSymbol+1):(2*SamplePerSymbol)),length(message),1);
pre_data=repelem(I(3:(end-1)),SamplePerSymbol);
modulated_phase=modulated_phase+pre_filt.*pre_data*4*pi*T*Fd;
%% before previuse symbol
bpre_filt=repmat(q_piece((2*SamplePerSymbol+1):end),length(message),1);
bpre_data=repelem(I(2:(end-2)),SamplePerSymbol);
modulated_phase=modulated_phase+bpre_filt.*bpre_data*4*pi*T*Fd;
%% phase_History
modulated_phase=modulated_phase+repelem(hist_phase,SamplePerSymbol,1)*2*pi*T*Fd;
modulated_data=exp(1i*(modulated_phase+phase_offset));




