function RxData=decoder_soft(Data_fin)

Data_fin=(Data_fin+1)/2;
ttt=poly2trellis(7,[171 133]);
Qcode=quantiz(Data_fin, [0.0094:0.0078:0.9922]);
%Qcode=quantiz(Data_fin, [0.0039:0.0039:0.9945]);
tblen=64; delay=tblen;
RxData=vitdec(Qcode, ttt,tblen,'cont','soft',7 );
