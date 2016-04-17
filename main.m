clc;
clear all;
close all;

Nfft=2048; 
Ng=512; 
Nofdm=2560; 
Nsym=100;
Nps=4; %Pilot Spacing
Np=Nfft/Nps; %Number of pilots per OFDM symbol
Nbps=4;
itr=10;
nn=11; %No. of SNR observations

%Initialising MSE arrays for different SNR observations:
z_linin=zeros(1,nn);
z_splin=zeros(1,nn);
z_lindft=zeros(1,nn);
z_spldft=zeros(1,nn);
zmmse=zeros(1,nn);
zmmse_dft=zeros(1,nn);

for t=1:itr
     z=[];
     z1=[];
     z2=[];
     z3=[];
     zm=[];
     zms=[];
     snr1=[];

    % 16 - QAM - Modulation Scheme 
     M=16;
     hmod = modem.qammod('M',M, 'SymbolOrder','gray');
     Es=1; 
     A=sqrt(3/2/(M-1)*Es); % Signal energy and QAM normalization factor
     for nsym=1:Nsym
        Xp = 2*(randn(1,Np)>0)-1; % Pilot sequence generation
        msgint=randi(1,Nfft-Np,M); % bit generation
        dat_ser = A*modulate(hmod,msgint);
     end

     % serial to parllel conversion
     dat_par=dat_ser.';

     % Pilot Insertion - Comb Type Arrangement
     counter = 0; 
     loc = [];
     for i=1:Nfft
         if mod(i,Nps)==1
             X(i)=Xp(floor(i/Nps)+1); 
             loc=[loc i]; 
             counter = counter+1;
         else
             X(i) = dat_par(i-counter);
         end
     end
     
     % inverse discret Fourier transform (IFFT)
     X_ifft=ifft(X,Nfft);
     
     % Adding Cyclic Prefix - Guard interval
     guard=X_ifft(:,end-511:end); % this is the Cyclic Prefix part to be appended.
     ofdm_par=[guard X_ifft];

     % parallel to serial - Generation of the OFDM Signal 
     ofdm=ofdm_par.';
     
     % Channel code
     dopjakes = doppler.jakes;

     dopgauss1 = doppler.bigaussian;
     dopgauss1.CenterFreqGaussian1 = -0.8;
     dopgauss1.CenterFreqGaussian2 = 0.4;
     dopgauss1.SigmaGaussian1 = 0.05;
     dopgauss1.SigmaGaussian2 = 0.1;
     dopgauss1.GainGaussian1 = sqrt(2*pi*(dopgauss1.SigmaGaussian1)^2);
     dopgauss1.GainGaussian2 = 1/10 * sqrt(2*pi*(dopgauss1.SigmaGaussian2)^2);

     dopgauss2 = doppler.bigaussian;
     dopgauss2.CenterFreqGaussian1 = 0.7;
     dopgauss2.CenterFreqGaussian2 = -0.4;
     dopgauss2.SigmaGaussian1 = 0.1;
     dopgauss2.SigmaGaussian2 = 0.15;
     dopgauss2.GainGaussian1 = sqrt(2*pi*(dopgauss1.SigmaGaussian1)^2);
     dopgauss2.GainGaussian2 = 1/10^1.5 * sqrt(2*pi*(dopgauss1.SigmaGaussian2)^2);
     
     fd=130; %Maximum Doppler Shift
     ts=(7/64)*10^-6; %Sampling Time

     chan = rayleighchan(ts, fd); %Rayleigh channel - Multifading Channel

     % Assign profile-specific properties to channel object.
     chan.PathDelays = [0.0 0.2 0.5 1.6 2.3 5.0] * 1e-6;
     chan.AvgPathGaindB = [-3 0 -2 -6 -8 -10];
     chan.DopplerSpectrum = [dopjakes dopjakes dopjakes dopgauss1 dopgauss2 dopgauss2];
     chan.StoreHistory = 1;
     chan.ResetBeforeFiltering = 0;
     chan.NormalizePathGains = 1;
     
     %Passing OFDM Signal through the created Channel 
     chan_op = filter(chan,ofdm);
     XFG=fft(chan_op);
     
     %Channel Covariance Matrix for MMSE estimation  
     x_test=randn(1,6);
     X_test=fft(x_test,Nfft);
     y_test=filter(chan,x_test);
     Y_test=fft(y_test,Nfft);
     H1=Y_test/X_test;
     h=ifft(H1,6);
     H=fft(h,Nfft);
     ch_length=length(h);
         
     %Reception at the receiver end for various SNR ranging from 5dB - 25dB:
     for SNR =5:2:25
         %AWGN Modelling
         snr1=[snr1 SNR];
         n1=ones(2560,1);
         n1=n1*0.000000000000000001i;%Just to ensure that the function awgn adds 'complex gaussian noise'..
         noise=awgn(n1,SNR);
         variance=var(noise);
         N=fft(noise);
         
         Y_rec=XFG+N;
         y_ser=ifft(Y_rec);
         
         % serial to parallel conversion
         y_par=y_ser.';
         
         % guard interval
         y = y_par(Ng+1:Nofdm);
         
         % FFT
         Y = fft(y);
         
         % channel estimation
         
         %LS Estimator with Linear Interpolator 
         H_est = LS_CE(Y,Xp,loc,Nfft,Nps,'linear');
         err=(H-H_est)*(H-H_est)';
         z=[z err/(Nfft*Nsym)];
         
         %LS Estimator with Linear Interpolator - DFT Based
         h_est = ifft(H_est);
         h_DFT = h_est(1:ch_length);
         H_DFT = fft(h_DFT,Nfft); 
         err=(H-H_DFT)*(H-H_DFT)';
         z3=[z3 err/(Nfft*Nsym)];
         
         %LS Estimator with Spline Cubic Interpolator
         H_est = LS_CE(Y,Xp,loc,Nfft,Nps,'spline');
         err=(H-H_est)*(H-H_est)';
         z1=[z1 err/(Nfft*Nsym)];
         
         %LS Estimator with Spline Cubic Interpolator - DFT Based
         h_est = ifft(H_est);
         h_DFT = h_est(1:ch_length);
         H_DFT = fft(h_DFT,Nfft); 
         err=(H-H_DFT)*(H-H_DFT)';
         z2=[z2 err/(Nfft*Nsym)];  
         
         % MMSE Estimator
         H_est = MMSE_CE(Y,Xp,loc,Nfft,Nps,h,SNR);
         err=(H-H_est)*(H-H_est)';
         zm=[zm err/(Nfft*Nsym)];
         
         % MMSE Estimator - DFT Based
         h_est = ifft(H_est);
         h_DFT = h_est(1:ch_length);
         H_DFT = fft(h_DFT,Nfft); 
         err=(H-H_DFT)*(H-H_DFT)';
         zms=[zms err/(Nfft*Nsym)];
     end
     
     z_linin=z_linin+z;
     z_splin=z_splin+z1;
     z_lindft=z_lindft+z2;
     z_spldft=z_spldft+z2;
     zmmse=zmmse+zm;
     zmmse_dft=zmmse_dft+zms;
end


figure(1)
 semilogy(snr1,(1/itr)*z_linin,'r+:', snr1,(1/itr)*z_splin,'bo:', snr1,(1/itr)*z_lindft,'--xg', snr1,(1/itr)*z_spldft,'--sc');
 legend('LS - Linear Interpolation','LS - Spline Cubic Interpolation','LS - Linear Interpolation(DFT)','LS - Spline Cubic Interpolation(DFT)');
 xlabel('SNR');
 ylabel('MSE');
 grid on
 hold on

figure(2)
 semilogy(snr1,(1/itr)*zmmse,'r+:',snr1,(1/itr)*zmmse_dft,'bo:');
 xlabel('SNR');
 ylabel('MSE');
 legend('MMSE','MMSE - DFT Based');
 hold on
 grid on
 

figure(3)
 semilogy(snr1,(1/itr)*z_linin,'r+:', snr1,(1/itr)*z_splin,'bo:', snr1, (1/itr)*zmmse,'--xg');
 xlabel('SNR');
 ylabel('MSE');
 legend('LS - Linear Interpolation','LS - Spline Cubic Interpolation','MMSE');
 hold on
 grid on 