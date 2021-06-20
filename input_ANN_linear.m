function [training_set]=input_ANN_linear(training_input,ff_TDL,sig_length)%ff_TDL=0,sig_length=tlen=2
% function to generate the ANN input
start=1;
finish=ff_TDL+1;
training_set=zeros(ff_TDL+1,sig_length);%1*2
for i=1:sig_length
    training_set(:,i)=training_input(start:finish)';
    start=start+1;
    finish=finish+1;
end