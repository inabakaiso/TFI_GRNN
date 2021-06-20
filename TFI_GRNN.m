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
phase=1;%DFCE�̌J��Ԃ���

threshold=0;
n1=5;
n2=10;
training=[1,n2];%�w�K��������͂���ы��t�M���ԍ�
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
            %  ���M�@  %
            %==========%
            bit_generation_vector=floor(rand(coding_ratio*number_symbol*number_carrier*modulation_index,1)*2)';
            
%             %==========%
%             %��������
%             
            bit_encode=encoder_soft(bit_generation_vector);
            %=====%
            %�C���^���[�o
            [temp,int_pattern]=sort(rand(1,size(bit_encode,2)));
            for xxx=1:size(bit_encode,2)
                bit_generation_nr_vector(1,xxx)=bit_encode(1,int_pattern(xxx));
            end
            %=====%
            bit_generation_nr=reshape(bit_generation_nr_vector,number_symbol,number_carrier,modulation_index);
            %==========%
           
            %�f�[�^�̓�l��
            data_nr=bit_generation_nr*2-1;
           
            %�p�C���b�g�M��������
            pilot_signal=repmat([1,0],1,ifft_size/2);
            
            %�ϒ���
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
            
            
            %IFFT������
            time_signal=(ifft(tx_baseband')).*sqrt(ifft_size);
                   
                    
            %�K�[�h�C���^�[�o���}��  -> �T�u�L�����A�� 49-64�܂ő}��
            %time_signal((ifft_size-guard_interval+1):ifft_size,:)  -> (16,21)
            tx_signal=[time_signal((ifft_size-guard_interval+1):ifft_size,:);time_signal];
      
            %������ϊ�
            serial_signal_tx=reshape(tx_signal,1,size(tx_signal,1)*size(tx_signal,2));
            
            %==========%
            %  �ʐM�H  %
            %==========%
            %�}���`�p�X�t�F�[�W���O
            faded_signal=multipath(serial_signal_tx,trans_rate,Doppler,number_path,decay_factor,path_interval);
            %�K�E�X�G��(AWGN)
            noise_dis=10^(-modified_eb_no/20);%�G�����U
            noise=sqrt(1/2)*(randn(1,length(faded_signal))+sqrt(-1)*randn(1,length(faded_signal))).*noise_dis;
            serial_signal_rx=faded_signal+noise;
            
            
            
%             %=======================
%             %     EVM�p����
%             %=======================
%             parallel_signal_faded=reshape(faded_signal,ifft_size+guard_interval,number_symbol+number_pilot);
%             %========�K�[�h�C���^�[�o������=======%
%             parallel_signal_faded=parallel_signal_faded(guard_interval+1:ifft_size+guard_interval,:);
%             %============FFT������==========%
%             frequency_signal_faded=fft(parallel_signal_faded)'/sqrt(ifft_size);
%             
%             frequency_signal_data_faded=frequency_signal_faded(number_pilot+1:number_pilot+number_symbol,:);
%             if modulation_index==1
%                 modulated_data_sample=data_nr;
%             else
%                 modulated_data_sample=data;
%             end
%             H_signal=frequency_signal_data_faded./modulated_data_sample;%���ۂ̃`���l�����ϓ��̎Z�o
%             
            %==========%
            %  ��M�@  %
            %==========%
            %������ϊ�
            parallel_signal=reshape(serial_signal_rx,ifft_size+guard_interval,number_symbol+number_pilot); %(64,21)
            %�K�[�h�C���^�[�o������
            parallel_signal=parallel_signal(guard_interval+1:ifft_size+guard_interval,:);
         
            %FFT������
            frequency_signal=fft(parallel_signal)'./sqrt(ifft_size);%(21,64)

            %%%descrambling
       %    frequency_signal=fft_signal./(scrambling_code);
            
            %�p�C���b�g���������o��
            H_resp = parallel_signal(:,number_pilot);
            
            %�`���l������
            %time windows - [0-16],[33-47]�܂ł𔲂��o���đ���0�Ŗ��߂�
            H_resp_time_average=H_resp(1:16,:)+H_resp(33:48,:); %(1,16)
            % 0�Ŗ��ߍ��킹
            channel_responce = [H_resp_time_average;zeros(ifft_size*3/4,number_pilot)];
            
            %fft����
%           channel_resp=fft(H_resp_time_average)./sqrt(ifft_size); %1,64
            channel_resp = (fft(channel_responce))'./sqrt(ifft_size);
            %�`���l���⏞���� 
            received_signal=frequency_signal(number_pilot+1:number_pilot+number_symbol,:)./repmat(channel_resp,number_symbol,1); %(20,64)
 
            %������
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
            %��������
            detected_bit_vector_pre=reshape(detected_bit,1,number_symbol*number_carrier*modulation_index);
            
            %=====%
            %�f�B�C���^���[�o
            for xxx=1:size(bit_encode,2)
                detected_bit_encode(1,int_pattern(xxx))=detected_bit_vector_pre(1,xxx);
            end
            %=====%
        
            detected_bit_vector=decoder_soft(detected_bit_encode*2-1);
            total_bit=length(bit_generation_vector);
            %==========%
            

            %%%----------------%%%
            %   ���Ԓn�_��藦�@�@%
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
                %�J��Ԃ���=========%
                %==========%
                %��������
                bit_encode_a=encoder_soft(detected_bit_vector_p);
                %=====%
                %�C���^���[�o
                for xxx=1:size(bit_encode_a,2)
                    bit_generation_nr_vector_a(1,xxx)=bit_encode_a(1,int_pattern(xxx));
                end
                %=====%
                bit_generation_nr_a=reshape(bit_generation_nr_vector_a,number_symbol,number_carrier,modulation_index);
               
                %==========%
                %�f�[�^�̓�l��
                data_nr_a=bit_generation_nr_a*2-1;
                %�ϒ���
                if modulation_index==1
                    tx_baseband_a=[data_nr_a];
                elseif modulation_index==2
                    data_a=(data_nr_a(:,:,1)+sqrt(-1).*data_nr_a(:,:,2))./sqrt(2);
                    tx_baseband_a=[data_a];
                end
                
                %��M�T���v���쐬
                fft_signal_a=channel_resp_dec.*tx_baseband_a;
                %��M�M������p�C���b�g�M������
                frequency_signal_a=frequency_signal(number_pilot+1:number_symbol+number_pilot,:);
                %��H�쐬
                for nn=1:number_symbol
                    for f=1:number_carrier
                        del_H(nn,f)=frequency_signal_a(nn,f)/(fft_signal_a(nn,f));
                    end
                end
                %�V�����`���l���W���쐬
                channel_resp_hat=del_H.*channel_resp_dec;
                %===�G������===%
                channel_resp_hat_b=[channel_resp;channel_resp_hat];
                
                for nn=1:number_symbol-1
                    channel_resp_hat_ave(nn,:)=(channel_resp_hat_b(nn,:)+channel_resp_hat_b(nn+1,:)+channel_resp_hat_b(nn+2,:))/3;
                end
                %=============%
                %���e�\�ȃ`���l�����U
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
                %��[�̂ݐ��`����
                channel_resp_dec_a=[channel_resp_dec_aa;channel_resp_behind];
                
                
                
                %===========ANN�w�K=========%
                %ANN�g���[�j���O�쐬
                %training=[1,n1,n2]; %% set at the beginning
                %ANN���́i�g���[�j���O��j
                input=1:1:number_symbol;
                %�^�[�Q�b�g�쐬
                target=channel_resp_dec_a;
                target(1,:)=channel_resp;%�擪�����̓p�C���b�g�����p����
                
                % ANN parameters
                sig_length=number_symbol;
                ff_TDL=0;
                % feed forward length
                ff_zeros=zeros(1,ff_TDL);
                
                %% Training Neural Network
                tlen=length(training);
                %trainig length
                
                %�C�R���C�U�̃g���[�j���O
                %=============================================
                %      �C�R���C�U�iANN-Based Linear�j
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
                net=newgrnn(training_input, OOK, spread);%�f�t�H���g��spread��1
                
                %ANN�C�R���C�U
                %=============================================
                %      �C�R���C�U�iANN-Based Linear�j
                %=============================================
                %% Simulation of ANN output
                ann_input=input_ANN_linear([ff_zeros input],ff_TDL,sig_length);
                
                ANN_output=sim(net,ann_input);
                channel_resp_dec_new(:,:,1)=reshape(ANN_output(1:number_carrier,:),number_carrier,number_symbol).';
                channel_resp_dec_new(:,:,2)=reshape(ANN_output(number_carrier+1:number_carrier*2,:),number_carrier,number_symbol).';
                
                channel_resp_dec_re=channel_resp_dec_new(:,:,1)+sqrt(-1).*channel_resp_dec_new(:,:,2);
                
                
                %�`���l���⏞
                received_signal_a=frequency_signal(number_pilot+1:number_symbol+number_pilot,:)./channel_resp_dec_re;%channel_resp_dec_re;
                
                %������
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
                %��������
                detected_bit_vector_pre_a=reshape(detected_bit_a,1,number_symbol*number_carrier*modulation_index);
                
                %=====%
                %�C���^���[�o
                for xxx=1:size(bit_encode_a,2)
                    detected_bit_encode_a(1,int_pattern(xxx))=detected_bit_vector_pre_a(1,xxx);
                end
                %=====%
                
                detected_bit_vector_a=decoder_soft(detected_bit_encode_a*2-1);
                total_bit=length(bit_generation_vector);
                %==========%
                %==========%
            end
            
            
            %���v�Z��
            error=find(bit_generation_vector(1,1:(total_bit-viterbi_path))-detected_bit_vector_a(1,viterbi_path+1:total_bit));
            number_error=length(error);
            disp(number_error);
            ber=number_error/(total_bit-viterbi_path);%BER�v�Z�l
            ber_result(ebno,1)=modified_eb_no;
            ber_result(ebno,2)=ber_result(ebno,2)+ber;
        end
        disp(['------------------------']);
    end
    
%     view(net);%�l�b�g���[�N�̎d�g�݂�\��

    
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
