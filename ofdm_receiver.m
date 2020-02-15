clc, clear, close all
%% TASK 1

fs_dec=12e3;
path='./SeñalesDRM/';  %directory of the signals
signals= dir([path '*.wav']);
file_n=2;
[inp_sig,fs]=audioread([path signals(file_n).name]);
N_fft=2^14;
filter_b_coeffs=[0.0056    0.0036   -0.0022   -0.0101   -0.0119   -0.0015    0.0165    0.0260    0.0118   -0.0241   -0.0553...
    -0.0445    0.0290    0.1471    0.2571    0.3019    0.2571    0.1471    0.0290   -0.0445   -0.0553   -0.0241...
    0.0118    0.0260    0.0165   -0.0015   -0.0119   -0.0101   -0.0022    0.0036    0.0056];
f_axis=-fs/2:fs/N_fft:fs/2-fs/N_fft;
figure, plot(f_axis./1e3,10*log10(abs(fftshift(fft(inp_sig, N_fft)))));
title('Received signal'),ylabel('Power(dB)'), xlabel('Freq (KHz)')
n=0:length(inp_sig)-1;
Ts=1/fs;
inp_sig_bb= inp_sig.*exp(-j*2*pi*Ts*n*12e3).';
figure, plot(f_axis./1e3,10*log10(abs(fftshift(fft(inp_sig_bb, N_fft)))));
%figure, plot(10*log10(abs(pwelch(inp_sig,[ ],[ ],fs)))));
title('Baseband shifted signal'),ylabel('Power(dB)'), xlabel('Freq (KHz)')
q= filter(filter_b_coeffs,1,inp_sig_bb);
figure, plot(f_axis./1e3,10*log10(abs(fftshift(fft(q, N_fft)))));
title('Filtered baseband signal (48KHz)'),ylabel('Power(dB)'), xlabel('Freq (KHz)')

f_axis_dec= -fs_dec/2: fs_dec/N_fft:fs_dec/2 - fs_dec/N_fft;
q_dec=q(1:4:end);
figure, plot(f_axis_dec./1e3,10*log10(abs(fftshift(fft(q_dec, N_fft)))));
title('Filtered baseband signal (12KHz)'),ylabel('Power(dB)'), xlabel('Freq (KHz)')

%% TASK 2

N_modes=4;
Nu=[288, 256 , 176, 112];
Ng=[32, 64,64,88];
Mode= ['Mode A'; 'Mode B'; 'Mode C'; 'Mode D'];
ini_pos=round(length(q_dec)/2); %choosing the start sample, where we are going to estimate the robustness mode
e=5000; %number of samples of the estimation

for i=1:N_modes  %one iteration per mode
    
    for j=ini_pos:ini_pos+e
        
        D(j-ini_pos+1,i)=abs(q_dec(j:Ng(i)-1+j).'*conj(q_dec(j+Nu(i):Ng(i)-1+j+Nu(i))))...
            -1/2*(sum(abs(q_dec(j:Ng(i)-1+j)).^2)+sum(abs(q_dec(j+Nu(i):Ng(i)-1+j+Nu(i))).^2));
        
    end
    %searching for periodic peaks
    [a,locs]=findpeaks(D(:,i), 'MinPeakWidth',25, 'MinPeakHeight', -0.1);
    n_peaks(i)=length(locs);
    [msg id] = lastwarn;
    warning('off',id) %switching off findpeaks warnings
end

figure,
subplot(221), plot(D(:,1)),xlim([0 e]),title('Mode A'),
subplot(222), plot(D(:,2)),xlim([0 e]),title('Mode B'),
subplot(223), plot(D(:,3)),xlim([0 e]),title('Mode C'),
subplot(224), plot(D(:,4)), xlim([0 e]),title('Mode D');


[~,pos]=max(n_peaks);   %Serching for the mode with the most peaks
disp([Mode(pos,:) ' whith Nu ' int2str(Nu(pos)) ' and Ng ' int2str(Ng(pos))])
[a,locs]=findpeaks(D(:,pos), 'MinPeakWidth',25, 'MinPeakHeight', -0.1);
%keeping the position of the maximums

Nu_chosen= Nu(pos);
Ng_chosen=Ng(pos);
acc=Nu_chosen+Ng_chosen; %Number of bits of every block
peak=locs(randi([1 length(locs)],1,1)); %we choose a random peak of the estimation

time_alignment=[ flip(ini_pos+peak:-acc:0) ini_pos+peak+acc:acc:length(q_dec)];
%The vector stores the positions where every OFDM block of the signal starts.


%% TASK 3
ini_pos= ini_pos+peak; %we just do the estimation in the peak that we have chosen before

omega_est=-1/(2*pi)*angle(q_dec(ini_pos:Ng_chosen-1+ini_pos).'*conj(q_dec(ini_pos+Nu_chosen:Ng_chosen-1+ini_pos+Nu_chosen)));
q_dec=q_dec.*exp(-1j*omega_est*2*pi*((0:length(q_dec)-1)/Nu_chosen)).';

for i=1:length(time_alignment)  %Every matrix's row corresponds to one OFDM block
    if time_alignment(i)+acc-1<length(q_dec)
        q_vector_axis(i,:)= q_dec(time_alignment(i)+Ng_chosen:Ng_chosen+time_alignment(i)+Nu_chosen-1);
    end
end


blocks_fft=fft(q_vector_axis, [] , 2); %FFT for every block

mean_power_blocks= sum(abs(blocks_fft),1)./length(blocks_fft); %mean power of all blocks

carrier_number=linspace(-Nu_chosen/2+1,Nu_chosen/2,Nu_chosen);
file_shift=[-5,-2];
mean_power_blocks= circshift(mean_power_blocks,file_shift(file_n));
figure, plot(carrier_number,fftshift(10*log10(mean_power_blocks))), title('Mean power spectrum'), ylabel('mean power (dB)'), xlabel('carrier number')

%% TASK 4

blocks_fft= circshift(blocks_fft,file_shift(file_n),2);
Ns=15;
carrier_index=[14,18,20,24,26,32,36,42,44,49,50,54,56,62,66,68];

for j=1:Ns
    
    frame_sync(j)=sum(sum(blocks_fft(j:Ns:end,carrier_index).*conj(blocks_fft(j:Ns:end,carrier_index))/length(carrier_index)))/Ns;
    
end
figure, stem(0:Ns-1,10*log10(frame_sync)), title('Frame alignment'), ylabel('Mean frame power taking the reference carriers (dB)'), xlabel('Frame number')


[~,loc]= max(frame_sync);

blocks_fft= blocks_fft(loc:end, :);

for j=1:Ns
    
    frame_sync(j)=sum(sum(blocks_fft(j:Ns:end,carrier_index).*conj(blocks_fft(j:Ns:end,carrier_index))/length(carrier_index)))/Ns;
    
end

figure, stem(0:Ns-1,10*log10(frame_sync)), title('Frame alignment after the estimation'), ylabel('Mean frame power taking the reference carriers (dB)'), xlabel('Frame number')




%% TASK 5

p=0:15;
x=2;
y=3;
k0=1;
W=[512,0,512,0,512;0,512,0,512,0;512,0,512,0,512];
Z=[0, 57, 164, 64, 12; 168, 255, 161, 106, 118; 25, 232, 132, 233, 38];
Q1024=12;



for s=0:14
    
    pilot_gain_index(s+1,:)= 1+2*mod(s,3)+6*p;
    
    for k=1:length(pilot_gain_index)
        n=mod(s,y);
        m=floor(s/y);
        p_int=(pilot_gain_index(s+1,k)-k0-n*x)/(x*y);
        pilot_gain_phase1024(s+1,k)= mod(4*Z(n+1,m+1)+ p_int*W(n+1,m+1)+p_int^2*(1+s)*Q1024,1024);
    end
    
end
pilot_gain_phase1024(1,9)=651;

pilot_gain_phase= 2*exp(1j*2*pi*pilot_gain_phase1024/1024);

FAC_index=[0,0,0,0,0,0;
    0,0,0,0,0,0;
    13,25,43,55,67,0;
    15,27,45,57,69,0;
    17, 29, 47, 59, 71,0;
    19, 31, 49, 61, 73,0;
    9, 21, 33, 51, 63, 75;
    11, 23, 35, 53, 65, 77;
    13, 25, 37, 55, 67, 79;
    15, 27, 39, 57, 69, 81;
    17, 29, 41, 59, 71, 83;
    19, 31, 43, 61, 73,0;
    21, 33, 45, 63, 75,0;
    23, 35, 47, 65, 77,0;
    0,0,0,0,0,0];

fac_symbols=[];
for i=1:15
    
    locs=find(FAC_index(i,:)~=0);
    
    if ~isempty(locs)
        ind_aux= FAC_index(i,locs);
        ind2_aux=[];
        
        for j=1:length(ind_aux)
            ind2_aux=[ind2_aux ,find(ind_aux(j)==pilot_gain_index(i+1,:))];
        end
        
        fac_block=blocks_fft(i:15:end, ind_aux);
        cell=blocks_fft(i+1:15:end, ind_aux);
        fac_block=fac_block(1:length(cell.'),:);
        [f,c]=size(fac_block);
        
        if(length(cell)==length(fac_block))
            inv_channel_response= repmat(pilot_gain_phase(i+1,ind2_aux),f,1)./cell;
            fac_symbols_aux= fac_block.*inv_channel_response;
            fac_symbols=[fac_symbols; fac_symbols_aux(:)];
        end
        
    end
    
end


scatterplot(fac_symbols), xlim([-2.5,2.5]), ylim([-2.5,2.5])


%% TASK 6
channel_response=zeros(13,floor(length(blocks_fft.')/15), length(pilot_gain_index));
symbols=[];

carrier_index=[16,48,64];
carrier_phase=[331,651,555];

carrier_index_def=horzcat(pilot_gain_index, repmat(carrier_index, 15,1));

[carrier_index_def, locs]= sort(carrier_index_def,2);
carrier_phase= horzcat(pilot_gain_phase1024, repmat(carrier_phase, 15,1));
for i=1:size(locs,1)
    carrier_phase_def(i,:)= 2*exp(1j*2*pi*carrier_phase(i,locs(i,:))/1024);
end

for s=2:14
    
    cells=blocks_fft(s+1:15:end, carrier_index_def(s+1,:));
    [f,~]=size(cells);
    channel_response= cells./repmat(carrier_phase_def(s+1,:),f,1);
    
    if s~=14
        channel_response2= cells./repmat(carrier_phase_def(s+2,:),f,1);
    end
    
    cell_channel=blocks_fft(s+1:15:end, 1:100);
    
    for k=1:length(channel_response)
        interp_channel_response=interp1(carrier_index_def(s+1,:),channel_response(k,:), 1:100, 'linear');
        if s~=14
            interp_channel_response(carrier_index_def(s+2,:))=channel_response2(k,:);
        end
        symbols=[symbols, 1./interp_channel_response.*cell_channel(k,:)];
    end
    
end


scatterplot(symbols(1:20000)), xlim([-2.5,2.5]), ylim([-2.5,2.5])
