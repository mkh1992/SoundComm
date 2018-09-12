% Decoder
function DecData=MessageDecoder(CodedData,FEC)
%% convolutional decoding
vitDecoder = comm.ViterbiDecoder(poly2trellis(7, [171 133]), ...
  'InputFormat', 'Hard','OutputDataType','uint32','TracebackDepth',39);
if strcmp(FEC,'1/2')
   MsgBchInter=vitDecoder(CodedData);
% elseif strcmp(FEC,'2/3')
%    vitDecoder.PuncturePatternSource =  'Property';
%    vitDecoder.PuncturePattern = [1;0;1;1]; 
%    MsgBchInter=vitDecoder(CodedData);
elseif strcmp(FEC,'3/4')
   vitDecoder.PuncturePatternSource =  'Property';
   vitDecoder.PuncturePattern = [1;0;1;1;1;0]; 
   MsgBchInter=vitDecoder(CodedData);
elseif strcmp(FEC,'5/6')
   vitDecoder.PuncturePatternSource =  'Property';
   vitDecoder.PuncturePattern = [1;0;1;0;1;1;1;0;1;0]; 
   MsgBchInter=vitDecoder(CodedData);
end
MsgBchInter=MsgBchInter(40:end);
%% deinterleaving
H = comm.ConvolutionalDeinterleaver('NumRegisters',3);
MsgBch=H(MsgBchInter);
MsgBch=MsgBch(13:end);
release(H);
%% BCH decoding
n=31;
k=26;
MsgBch=reshape(MsgBch,29,length(MsgBch)/29)';
[s,~]=size(MsgBch);
MsgBch=[zeros(s,2),MsgBch];
MsgBch=gf(MsgBch);
Msg=bchdec(MsgBch,n,k);
Msg=Msg.x;
Msg=Msg(:,3:end);
[s,~]=size(Msg);
DecData=reshape(Msg',8,3*s)';

