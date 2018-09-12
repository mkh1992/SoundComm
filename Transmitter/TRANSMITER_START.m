%% TEXT to HEX converter
function [HEX,CenterFreq,SymbolRate,FEC]=TRANSMITER_START()
prompt = {'Enter desired message','Link frequency','Symbol Rate:','FEC:'};
title = 'Settings';
dims = [5 35;1 35;1 35;1 35];
definput = {'Write your message here','10e3','500','3/4'};
answer = inputdlg(prompt,title,dims,definput);
if isempty(answer)
   msgbox('you cancelled the proccess','Error','warn');
elseif sum(sum(answer{4,1}==['1/2';'3/4';'5/6']))~=5
    msgbox('desired FEC CodeRate is not supported','Error','error')
elseif ~(str2double(answer{2,1})<18e3 && str2double(answer{2,1})>2e3 )
     msgbox('Center Frequency range must be within 2-18 KHz ','Error','warn')
elseif str2double(answer{3,1})>1e3
     msgbox('Symbol Rate must be Less than 1Ks/s ','Error','warn')
     
else
    CenterFreq=str2double(answer{2,1});
    SymbolRate=str2double(answer{3,1});
    FEC=answer{4,1};
    HEX=answer{1,1}+0;
    HEX=dec2hex(HEX);
end
