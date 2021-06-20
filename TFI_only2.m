%% GRNN TFI scramblingあり

clear all;
tic;
rand('state',sum(100*clock));
ifft_size=64;
guard_interval=16;
carrier_interval=1;
number_carrier=64;
number_symbol=20;
number_pilot=1;
modulation_index=2;
trans_rate=20*10^6;
Doppler=800;
number_path=15;
path_interval=1;
decay_factor=1;
coding_ratio=1/2;
viterbi_path=2^(7-1);
ebno_min=0;
ebno_max=30;
ebno_step=5;
repeat=40000;
ber_result=zeros(floor((ebno_max-ebno_min)/ebno_step)+1,2);
phase=1;%DFCEの繰り返す回数

threshold=0;
n1=5;
n2=10;
training=[1,n2];%学習させる入力および教師信号番号
viterbi_path=2^(7-1);

 pilot_signal = repmat([1,0], 1, ifft_size/2); 
%%
path_amplitude=zeros(number_path,1);
for tap = 1: number_path
    path_amplitude(tap)=10^(0-(((tap-1)*decay_factor)/20));
end
path_power = path_amplitude.^2;
normalization=sum(path_power);
norm_path_amplitude=path_amplitude./sqrt(normalization);

evm_rms=zeros(floor((ebno_max-ebno_min)/ebno_step)+1,1);
evm_rms_time=zeros(number_symbol,1);
%%


 for ebno=1:((ebno_max-ebno_min)/ebno_step)+1
        eb_no=(ebno-1)*ebno_step+ebno_min;
        
        modified_eb_no=eb_no+10*log10(modulation_index)+10*log10(ifft_size/(ifft_size+guard_interval))+10*log10(number_symbol/(number_symbol+number_pilot));% +10*log10(coding_ratio);
%          modified_eb_no=eb_no+10*log10(modulation_index)+10*log10(1/2)+10*log10(64/80)+10*log10(20/20.5);
      
            %==========%
            %  送信機  %
            %==========%
    for  R=1:repeat
        disp(['Repeat', num2str(R)]);
        
%%  ber特性　-> 5%ほど改善
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %送信信号生成
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bit_generation_vector=permute(floor(rand(coding_ratio*number_symbol*number_carrier*modulation_index,1)*2),[2 1]);
        bit_encode=encoder_soft(bit_generation_vector);
        [temp, int_pattern]=sort(rand(1,size(bit_encode,2)));
        for xxx=1:size(bit_encode,2)
            bit_generation_nr_vector(1,xxx)=bit_encode(1,int_pattern(:,xxx));
        end
        bit_generation_nr = reshape(bit_generation_nr_vector,number_symbol,number_carrier,modulation_index);
      
        data_nr = bit_generation_nr * 2 -1;
        data = (data_nr(:,:,1) + 1i*data_nr(:,:,2))./sqrt(2);
        tx_baseband = [pilot_signal;data];
        time_signal = ifft(tx_baseband').*sqrt(ifft_size);
        
        tx_signal = [time_signal(ifft_size - guard_interval + 1 : ifft_size,:); time_signal];
        serial_signal_tx = reshape(tx_signal,1,size(tx_signal,1)*size(tx_signal,2));
        
        faded_signal = multipath(serial_signal_tx,trans_rate,Doppler,number_path,decay_factor,path_interval);
        
        noise_dis=10^(-modified_eb_no/20);%雑音分散
        noise = sqrt(1/2)*(randn(1,length(faded_signal)) + 1i*randn(1,length(faded_signal))).*noise_dis;
        serial_signal_rx = faded_signal + noise;
        
        
        % 直並列変換
        parallel_signal = reshape(serial_signal_rx,ifft_size + guard_interval,number_symbol + number_pilot);
        %　ガードインタバル除去
        parallel_signal = parallel_signal(guard_interval + 1:ifft_size + guard_interval,:);
        % FFT処理
        fft_signal = (fft(parallel_signal))'./sqrt(ifft_size);
        
        
        H_resp = parallel_signal(:,number_pilot);
        channel_resp = H_resp(1:16,:) + H_resp(33:48,:);
        channel_responce = [channel_resp;zeros(ifft_size*3/4,number_pilot)];
        H = (fft(channel_responce))'./sqrt(ifft_size);
        received_signal = fft_signal(number_pilot + 1: number_pilot + number_symbol,:)./repmat(H,number_symbol,1);


%%
                %復調部
                for k=1:number_symbol
                    for cc=1:number_carrier
                        if modulation_index==1
                            X_a=received_signal(k,cc);
                            Xi_a=sign(real(X));
                            detected_bit_a(k,cc)=(Xi_a+1)/2;
                        elseif modulation_index==2
                            X_a=received_signal(k,cc);
                            Xi_a=sign(real(X_a));
                            Xq_a=sign(imag(X_a));
                            detected_bit_a(k,cc,1)=(Xi_a+1)/2;
                            detected_bit_a(k,cc,2)=(Xq_a+1)/2;
                        end
                    end
                end
                
                %==========%
                %誤り訂正部
                detected_bit_vector_pre_a=reshape(detected_bit_a,1,number_symbol*number_carrier*modulation_index);
                
                %=====%
                %ディインタリーバ
                for xxx=1:size(bit_encode,2)
                    detected_bit_encode_a(1,int_pattern(xxx))=detected_bit_vector_pre_a(1,xxx);
                end
                %=====%
                
                detected_bit_vector_a=decoder_soft(detected_bit_encode_a*2-1);
                total_bit=length(bit_generation_vector);
                %==========%
                %==========%
    
            
            
            %誤り計算部
            error=find(bit_generation_vector(1,1:(total_bit-viterbi_path))-detected_bit_vector_a(1,viterbi_path+1:total_bit));
            number_error=length(error);
            disp(number_error);
            ber=number_error/(total_bit-viterbi_path);%BER計算値
            ber_result(ebno,1)=modified_eb_no;
            ber_result(ebno,2)=ber_result(ebno,2)+ber;
         end
      
        disp(['------------------------']);
    end
    
%     view(net);%ネットワークの仕組みを表示

    
    ber_result(:,2)=ber_result(:,2)/repeat;
    figure(2)
    legend('TFI-GRNN')
    title('BER');
    semilogy(ebno_min:ebno_step:ebno_max,ber_result(:,2)','-mo');
    xlabel('Eb/No[dB]');
    ylabel('BER');
    axis([0,30,10^(-6),10^(0)]);
    grid;
    hold on
tocf(3);
