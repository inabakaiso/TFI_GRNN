%% GRNN TFI

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
Doppler=10;
number_path=15;
path_interval=1;
decay_factor=1;
coding_ratio=1/2;
viterbi_path=2^(7-1);
ebno_min=0;
ebno_max=30;
ebno_step=5;
repeat=2000;
ber_result=zeros(floor((ebno_max-ebno_min)/ebno_step)+1,2);
phase=1;%DFCEの繰り返す回数

threshold=0;
n1=5;
n2=10;
training=[1,n2];%学習させる入力および教師信号番号
viterbi_path=2^(7-1);

evm_rms=zeros(floor((ebno_max-ebno_min)/ebno_step)+1,1);
evm_rms_time=zeros(number_symbol,1);
for spread=5;%(1)/sqrt(2*2);%5 
    for ebno=1:((ebno_max-ebno_min)/ebno_step)+1
        eb_no=(ebno-1)*ebno_step+ebno_min;
        
%         modified_eb_no=eb_no+10*log10(modulation_index)+10*log10(ifft_size/(ifft_size+guard_interval))+10*log10(number_symbol/(number_symbol+number_pilot));% +10*log10(coding_ratio);
         modified_eb_no=eb_no+10*log10(modulation_index)+10*log10(1/2)+10*log10(64/80)+10*log10(20/20.5);
        for R=1:repeat
            disp(['spread:',num2str(spread) 'Eb.No:',num2str(eb_no) '     ' 'Repeat:',num2str(R)]);
            %==========%
            %  送信機  %
            %==========%
            bit_generation_vector=floor(rand(coding_ratio*number_symbol*number_carrier*modulation_index,1)*2)';
            
%             %==========%
%             %誤り訂正部
%             
            bit_encode=encoder_soft(bit_generation_vector);
            %=====%
            %インタリーバ
            [temp,int_pattern]=sort(rand(1,size(bit_encode,2)));
            for xxx=1:size(bit_encode,2)
                bit_generation_nr_vector(1,xxx)=bit_encode(1,int_pattern(xxx));
            end
            %=====%
            bit_generation_nr=reshape(bit_generation_nr_vector,number_symbol,number_carrier,modulation_index);
            %==========%
           
            %データの二値化
            data_nr=bit_generation_nr*2-1;
           
            %パイロット信号発生部
            pilot_signal=repmat([1,0],1,ifft_size/2);
            
            %変調部
            if modulation_index==1
                tx_baseband=[pilot_signal;data_nr];
            elseif modulation_index==2%QPSK
                data=(data_nr(:,:,1)+sqrt(-1).*data_nr(:,:,2))./sqrt(2);
                tx_baseband=[pilot_signal;data];
            end
            
            %%%scarmbling
%             scrambling_code=2*randint(21,64)-1;
            % scrambling_code=2*randi([0,1],21,64)-1;
           
            % tx_baseband=scrambling_code.*tx_baseband;
            
            
            %IFFT処理部
            time_signal=(ifft(tx_baseband')).*sqrt(ifft_size);
                   
                    
            %ガードインターバル挿入  -> サブキャリアの 49-64まで挿入
            %time_signal((ifft_size-guard_interval+1):ifft_size,:)  -> (16,21)
            tx_signal=[time_signal((ifft_size-guard_interval+1):ifft_size,:);time_signal];
      
            %並直列変換
            serial_signal_tx=reshape(tx_signal,1,size(tx_signal,1)*size(tx_signal,2));
            
            %==========%
            %  通信路  %
            %==========%
            %マルチパスフェージング
            faded_signal=multipath(serial_signal_tx,trans_rate,Doppler,number_path,decay_factor,path_interval);
            %ガウス雑音(AWGN)
            noise_dis=10^(-modified_eb_no/20);%雑音分散
            noise=sqrt(1/2)*(randn(1,length(faded_signal))+sqrt(-1)*randn(1,length(faded_signal))).*noise_dis;
            serial_signal_rx=faded_signal+noise;
            
            
            
%             %=======================
%             %     EVM用処理
%             %=======================
%             parallel_signal_faded=reshape(faded_signal,ifft_size+guard_interval,number_symbol+number_pilot);
%             %========ガードインターバル除去=======%
%             parallel_signal_faded=parallel_signal_faded(guard_interval+1:ifft_size+guard_interval,:);
%             %============FFT処理部==========%
%             frequency_signal_faded=fft(parallel_signal_faded)'/sqrt(ifft_size);
%             
%             frequency_signal_data_faded=frequency_signal_faded(number_pilot+1:number_pilot+number_symbol,:);
%             if modulation_index==1
%                 modulated_data_sample=data_nr;
%             else
%                 modulated_data_sample=data;
%             end
%             H_signal=frequency_signal_data_faded./modulated_data_sample;%実際のチャネル時変動の算出
%             
            %==========%
            %  受信機  %
            %==========%
            %直並列変換
            parallel_signal=reshape(serial_signal_rx,ifft_size+guard_interval,number_symbol+number_pilot); %(64,21)
            %ガードインターバル除去
            parallel_signal=parallel_signal(guard_interval+1:ifft_size+guard_interval,:);
         
            %FFT処理部
            frequency_signal=fft(parallel_signal)'./sqrt(ifft_size);%(21,64)

            %%%descrambling
       %    frequency_signal=fft_signal./(scrambling_code);
            
            %パイロット部分抜き出し
            H_resp = parallel_signal(:,number_pilot);
            
            %チャネル推定
            %time windows - [0-16],[33-47]までを抜き出して他は0で埋める
            H_resp_time_average=H_resp(1:16,:)+H_resp(33:48,:); %(1,16)
            % 0で埋め合わせ
            channel_responce = [H_resp_time_average;zeros(ifft_size*3/4,number_pilot)];
            
            %fft処理
%           channel_resp=fft(H_resp_time_average)./sqrt(ifft_size); %1,64
            channel_resp = (fft(channel_responce))'./sqrt(ifft_size);
            %チャネル補償処理 
            received_signal=frequency_signal(number_pilot+1:number_pilot+number_symbol,:)./repmat(channel_resp,number_symbol,1); %(20,64)
 
            %復調部
            for k=1:number_symbol
                for cc=1:number_carrier
                    if modulation_index==1
                        X=received_signal(k,cc);
                        Xi=sign(real(X));
                        detected_bit(k,cc)=(Xi+1)/2;
                    elseif modulation_index==2
                        X=received_signal(k,cc);
                        Xi=sign(real(X));
                        Xq=sign(imag(X));
                        detected_bit(k,cc,1)=(Xi+1)/2;
                        detected_bit(k,cc,2)=(Xq+1)/2;
                    end
                end
            end
            
            %==========%
            %誤り訂正部
            detected_bit_vector_pre=reshape(detected_bit,1,number_symbol*number_carrier*modulation_index);
            
            %=====%
            %ディインタリーバ
            for xxx=1:size(bit_encode,2)
                detected_bit_encode(1,int_pattern(xxx))=detected_bit_vector_pre(1,xxx);
            end
            %=====%
        
            detected_bit_vector=decoder_soft(detected_bit_encode*2-1);
            total_bit=length(bit_generation_vector);
            %==========%
            

            %%%----------------%%%
            %   中間地点誤り率　　%
            %%%----------------%%%
%             length_bit = length(bit_generation_vector);
%             error = length(find(bit_generation_vector(:,1:length_bit-viterbi_path) - detected_bit_vector(:,viterbi_path+1:length_bit)));
%             BER = error/(total_bit - viterbi_path);
%             ber_result(ebno)=ber_result(ebno)+BER;
            


            channel_resp_dec=repmat(channel_resp,number_symbol,1);
            
            for Phase=1:phase
                if Phase==1
                    detected_bit_vector_p=[detected_bit_vector(:,viterbi_path+1:total_bit),detected_bit_vector(:,1:viterbi_path)];
                else
                    detected_bit_vector_p=[detected_bit_vector_a(:,viterbi_path+1:total_bit),detected_bit_vector_a(:,1:viterbi_path)];
                end
                %繰り返し部=========%
                %==========%
                %誤り訂正部
                bit_encode_a=encoder_soft(detected_bit_vector_p);
                %=====%
                %インタリーバ
                for xxx=1:size(bit_encode_a,2)
                    bit_generation_nr_vector_a(1,xxx)=bit_encode_a(1,int_pattern(xxx));
                end
                %=====%
                bit_generation_nr_a=reshape(bit_generation_nr_vector_a,number_symbol,number_carrier,modulation_index);
               
                %==========%
                %データの二値化
                data_nr_a=bit_generation_nr_a*2-1;
                %変調部
                if modulation_index==1
                    tx_baseband_a=[data_nr_a];
                elseif modulation_index==2
                    data_a=(data_nr_a(:,:,1)+sqrt(-1).*data_nr_a(:,:,2))./sqrt(2);
                    tx_baseband_a=[data_a];
                end
                
                %受信サンプル作成
                fft_signal_a=channel_resp_dec.*tx_baseband_a;
                %受信信号からパイロット信号除去
                frequency_signal_a=frequency_signal(number_pilot+1:number_symbol+number_pilot,:);
                %⊿H作成
                for nn=1:number_symbol
                    for f=1:number_carrier
                        del_H(nn,f)=frequency_signal_a(nn,f)/(fft_signal_a(nn,f));
                    end
                end
                %新しいチャネル係数作成
                channel_resp_hat=del_H.*channel_resp_dec;
                %===雑音除去===%
                channel_resp_hat_b=[channel_resp;channel_resp_hat];
                
                for nn=1:number_symbol-1
                    channel_resp_hat_ave(nn,:)=(channel_resp_hat_b(nn,:)+channel_resp_hat_b(nn+1,:)+channel_resp_hat_b(nn+2,:))/3;
                end
                %=============%
                %許容可能なチャネル分散
                for f=1:number_carrier
                    for nn=1:number_symbol-1
                        norm(nn,f)=(abs(channel_resp_hat_ave(nn,f))^2)/(abs(channel_resp_dec(nn,f))^2);
                        if norm(nn,f)>=threshold
                            channel_resp_dec_aa(nn,f)=channel_resp_hat_ave(nn,f);
                        else
                            channel_resp_dec_aa(nn,f)=channel_resp_dec(nn,f);
                        end
                    end
                end
                channel_resp_behind=(channel_resp_dec_aa(number_symbol-1,:)-channel_resp_dec_aa(number_symbol-2,:))+channel_resp_dec_aa(number_symbol-1,:);
                %後端のみ線形推定
                channel_resp_dec_a=[channel_resp_dec_aa;channel_resp_behind];
                
                
                
                %===========ANN学習=========%
                %ANNトレーニング作成
                %training=[1,n1,n2]; %% set at the beginning
                %ANN入力（トレーニング後）
                input=1:1:number_symbol;
                %ターゲット作成
                target=channel_resp_dec_a;
                target(1,:)=channel_resp;%先頭だけはパイロット推定を用いる
                
                % ANN parameters
                sig_length=number_symbol;
                ff_TDL=0;
                % feed forward length
                ff_zeros=zeros(1,ff_TDL);
                
                %% Training Neural Network
                tlen=length(training);
                %trainig length
                
                %イコライザのトレーニング
                %=============================================
                %      イコライザ（ANN-Based Linear）
                %=============================================
                for cc=1:number_carrier
                    if cc==1
                        OOK=real(target(training,1)).';
                    else
                        OOK=[OOK;real(target(training,cc)).'];
                    end
                end
                for cc=1:number_carrier
                    if cc==1
                        OOK=[OOK;imag(target(training,1)).'];
                    else
                        OOK=[OOK;imag(target(training,cc)).'];
                    end
                end
                
                training_input=input_ANN_linear([ff_zeros training],ff_TDL,tlen);
                % training input
                net=newgrnn(training_input, OOK, spread);%デフォルトのspreadは1
                
                %ANNイコライザ
                %=============================================
                %      イコライザ（ANN-Based Linear）
                %=============================================
                %% Simulation of ANN output
                ann_input=input_ANN_linear([ff_zeros input],ff_TDL,sig_length);
                
                ANN_output=sim(net,ann_input);
                channel_resp_dec_new(:,:,1)=reshape(ANN_output(1:number_carrier,:),number_carrier,number_symbol).';
                channel_resp_dec_new(:,:,2)=reshape(ANN_output(number_carrier+1:number_carrier*2,:),number_carrier,number_symbol).';
                
                channel_resp_dec_re=channel_resp_dec_new(:,:,1)+sqrt(-1).*channel_resp_dec_new(:,:,2);
                
                
                %チャネル補償
                received_signal_a=frequency_signal(number_pilot+1:number_symbol+number_pilot,:)./channel_resp_dec_re;%channel_resp_dec_re;
                
                %復調部
                for k=1:number_symbol
                    for cc=1:number_carrier
                        if modulation_index==1
                            X_a=received_signal_a(k,cc);
                            Xi_a=sign(real(X));
                            detected_bit_a(k,cc)=(Xi_a+1)/2;
                        elseif modulation_index==2
                            X_a=received_signal_a(k,cc);
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
                %インタリーバ
                for xxx=1:size(bit_encode_a,2)
                    detected_bit_encode_a(1,int_pattern(xxx))=detected_bit_vector_pre_a(1,xxx);
                end
                %=====%
                
                detected_bit_vector_a=decoder_soft(detected_bit_encode_a*2-1);
                total_bit=length(bit_generation_vector);
                %==========%
                %==========%
            end
            
            
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
    figure(31)
    legend('TFI-GRNN')
    title('BER');
    semilogy(ebno_min:ebno_step:ebno_max,ber_result(:,2)','-bo');
    xlabel('Eb/No[dB]');
    ylabel('BER');
    axis([0,30,10^(-6),10^(0)]);
    grid;
    hold on

end
tocf(3);
