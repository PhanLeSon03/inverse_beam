clear;
clc;
Nfft = 256;               % FFT length used to implement FFT filterbank
Nfh = Nfft/2+1;           % number of frequency points in [0,Fs/2]
Fs = 48000;
Fsh = Fs/2;
fl = 350;                 % lower and upper cutoff frequencies of filterbank
fu = 20000;
klow = round(fl/Fsh*Nfh);
kup = round(fu/Fsh*Nfh);
f = Fsh/Nfh*(klow:kup)';  % frequencies used to compute T and W
Nf = length(f);
N = 101;                   % number of microphone
dH = 0.01;            % minimum sensor distance in m
gamma = 1;
numCluster = 5;
numPrsnt = 2;


c = 343.2;                 % speed of sound in m/s
f_test = 2000;

x_map = zeros(N,Nf);
n = -(N-1)/2:(N-1)/2;

x_array = n*dH;
theta = linspace(-pi/2, pi/2, 180); % discrete Theta angle 
Ntheta = length(theta);


BP = 0.0307*exp(-1j*3*pi*sin(theta)) + 0.2028*exp(-1j*2*pi*sin(theta)) + ...
    0.1663*exp(-1j*1*pi*sin(theta)) + ...
    0.2004*exp(-1j*0*pi*sin(theta)) + ...
    0.1663*exp(1j*1*pi*sin(theta)) + ...
    0.2028*exp(1j*2*pi*sin(theta)) + ...
    0.0307*exp(1j*3*pi*sin(theta));
BP = BP/max(abs(BP));
BP_ML = zeros(1,length(BP));
BP_ML(69:112) = BP(69:112); 
figure();

plot(theta*180/pi,abs(BP));

hold on 
plot(theta*180/pi,abs(BP_ML),'*r');
ylabel('Gain');
xlabel('Incident angle (\theta)');
set(gcf,'color','w');
legend('BP with side-lobe','BP without side-lobe');
title('reference BP for line array')
set(gca,'FontSize', 12);
axis tight
BP = BP_ML;
%%
FI = zeros(Ntheta,Nf);
figure()

for i=1:Nf
    
  Br = zeros(N,1);  
  Rc = (f(i)*N*dH/c);
  %f(i) = Rc*c/N/dH;
  x_points = Rc*cos(theta);
  %boudary  = find(abs(n)<= Rc  );
  %outside  = find(abs(n)> Rc  );
  boudary = find(abs(n)<= Rc*sind(90-69)  );
  outside  = find(abs(n)> Rc*sind(90-69)  );
  thetaS = asin(n(boudary)/Rc);
  

  % reference beam-pattern
  Br(boudary)= 0.0307*exp(-1j*3*pi*sin(thetaS)) + ...
    0.2028*exp(-1j*2*pi*sin(thetaS)) + ...
    0.1663*exp(-1j*1*pi*sin(thetaS)) + ...
    0.2004*exp(-1j*0*pi*sin(thetaS)) + ...
    0.1663*exp(1j*1*pi*sin(thetaS)) + ...
    0.2028*exp(1j*2*pi*sin(thetaS)) + ...
    0.0307*exp(1j*3*pi*sin(thetaS));
   %Br(outside)= abs(Br(boudary(end)));
   Br(boudary) = Br(boudary)/max(Br(boudary));

    
   % animation plot
   plot(abs(Br)) 
   temp = fftshift(Br(end:-1:1));
   h = fftshift(ifft(temp(end:-1:1))); 
    %   MagMax = max(abs(h));
    %   idxZero = find (abs(h)<0.05*MagMax);
    %   h(idxZero) = 0;
   x_map(:,i) = h;
  
  
  %hold on
  pause(0.01)
  
  %beam plot
  beta = 2*pi*f(i)/c;             % wave number
  D = exp(1j*beta*x_array(ones(1,Ntheta),:).*sin(theta(ones(N,1),:))');
  FI(:,i) = (D*(h));     
end
%%
pos = [0.5 0.5 0.4 0.4];

figure('numbertitle','off','name','beam pattern of full array',...
      'Units','normal','Position',pos);
surf(180/pi*theta,f,abs(FI)');
axis tight
set(gca,'XTick',[0 45 90 135 180]);
view([25,50]);
xlabel('Azimuth \theta in deg');
ylabel('freqeuncy in Hz');
zlabel('Magnitude');
grid on
set(gcf,'color','w');

figure()
k_p = round((16000-fl)/Fsh*Nfh)+1;
p = linspace(-pi, pi, 360);
dBmin = -50;
R_ref = abs([BP BP(end:-1:1)])';
RdB = max(dBmin,20*log10(R_ref));
polarplot(p,RdB','k'); 

hold on

R = abs([FI(end:-1:1,k_p);FI(:,k_p)]);
RdB = max(dBmin,20*log10(R)) ;
polarplot(p,RdB','--r'); 

axis tight;
%title(Plot_Title(iConfig));
set(gca,'FontSize', 12);
legend('expected BP','real BP');
set(gcf,'color','w');
