function Data=encoder_soft(Ori_data)


ttt=poly2trellis(7,[171 133]);
Data = convenc(Ori_data,ttt);
