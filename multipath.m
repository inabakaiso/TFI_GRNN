function faded_signal=multipath(input_signal,trans_rate,Doppler,number_path,decay_factor,path_interval)

rand('state',sum(100*clock));%rand関数の初期化
number_symbol=length(input_signal);
Ts=1/trans_rate;%シンボル周期
N_0=16;%素波数
path_amplitude=zeros(number_path,1);
for tap=1:number_path
path_amplitude(tap)=10^(0-(((tap-1)*decay_factor)/20));
randtime=(10^10)*rand(1);% ランダム初期タイム
Omega_m=2*pi*Doppler;
N=4*N_0+2;
xc=zeros(1,number_symbol);
xs=zeros(1,number_symbol);
rand_point=[1:number_symbol]+randtime;
for t=1:N_0
wm=cos(cos(2*pi*t/N)*rand_point*Omega_m*Ts);
xc=xc+cos((pi/N_0)*t)*wm;
xs=xs+sin((pi/N_0)*t)*wm;
end
T=sqrt(2)*cos(rand_point*Omega_m*Ts);
xc=(2*xc+T)*sqrt(1/(2*(N_0+1)));
xs=2*xs*sqrt(1/(2*N_0));
fade(tap,:)=xc+sqrt(-1).*xs;
end
normalization=sum(path_amplitude.^2);
norm_path_amplitude=path_amplitude./sqrt(normalization);
signal1=repmat(input_signal,number_path,1);
signal2=(signal1.*repmat(norm_path_amplitude,1,number_symbol)).*fade;
for a=1:number_path
fade_signal(a,1:number_symbol)=shift(signal2(a,:),-(path_interval*(a-1)));
end
faded_signal=sum(fade_signal);