function [FrameStart,FrameEnd,Steps]=PackMan(Data_size,Buf_size,Buf_overlap,index)
% Data_size   : The Array-length of the input data buffer
% Buf_size    : sub-array-size, who can be the FFT-size
% Buf_overlap : Num of samples shift pr. frame. Buf_size/2(pow-calc) or Buf_size/4(FDcomp)

FrameStart=(index-1)*Buf_overlap+1;
FrameEnd=FrameStart+Buf_size-1;
Steps=ceil((Data_size-Buf_size)/Buf_overlap);

