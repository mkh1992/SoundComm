% CODER  
%BCH code(29,24)  ----->   convolutional interleaving  ------>  convolutional FEC code
function CodedData=MessageCoder(HEX,FEC)
message=hexToBinaryVector(HEX,8);
% dec=hex2dec(HEX);
% message=de2bi(dec,8,'left-msb');
[s,~]=size(message);
CombMess=reshape(message',24,ceil(s/3))'; % combine each 3 byte for codding
%% BCH SETTING
CombMess=[zeros(ceil(s/3),2),CombMess]; % appending zeros first of data
CombMess=gf(CombMess);
n=31;
k=26;
BchEnc = bchenc(CombMess,n,k);
MsgBch=BchEnc.x;
MsgBch=MsgBch(:,3:end); % removing Zeros first of code
MsgBch=MsgBch';
MsgBch=MsgBch(:);
%% convolutional interleaving
H = comm.ConvolutionalInterleaver('NumRegisters',3);
MsgBchInter=H([MsgBch;randi([0,1],12,1)]);
%% convolutional Coding FEC
% mother code
convEncoder = comm.ConvolutionalEncoder(poly2trellis(7, [171 133])); % Generates standard 1/2 code
MsgBchInter=[MsgBchInter;randi([0 1],39,1)];
if strcmp(FEC,'1/2')
   CodedData=convEncoder(MsgBchInter);
% elseif strcmp(FEC,'2/3')
%    convEncoder.PuncturePatternSource = 'Property';
%    convEncoder.PuncturePattern = [1;0;1;1]; 
%    CodedData=convEncoder(MsgBchInter);
elseif strcmp(FEC,'3/4')
   convEncoder.PuncturePatternSource = 'Property';
   convEncoder.PuncturePattern = [1;0;1;1;1;0]; 
   CodedData=convEncoder(MsgBchInter);
elseif strcmp(FEC,'5/6')
   convEncoder.PuncturePatternSource = 'Property';
   convEncoder.PuncturePattern = [1;0;1;0;1;1;1;0;1;0]; 
   CodedData=convEncoder(MsgBchInter);
end
