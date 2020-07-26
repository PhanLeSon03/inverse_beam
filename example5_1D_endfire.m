clear;
clc;
Nfft = 128;               % FFT length used to implement FFT filterbank
Nfh = Nfft/2+1;           % number of frequency points in [0,Fs/2]
Fs = 16000;
Fsh = Fs/2;
fl = 1000;                 % lower and upper cutoff frequencies of filterbank
fu = 6000;
klow = round(fl/Fsh*Nfh);
kup = round(fu/Fsh*Nfh);
f = Fsh/Nfh*(klow:kup)';  % frequencies used to compute T and W
Nf = length(f);
% N = 141;                   % number of microphone
% dH = 0.01;            % minimum sensor distance in m
% numCluster = 9;
% numPrsnt = 1;
% gamma = 0.9;

N = 141;                   % number of microphone
dH = 0.01;            % minimum sensor distance in m
numCluster = 7;
numPrsnt = 1;
gamma = 0.5;

c = 340;                 % speed of sound in m/s
f_test = 2000;

x_map = zeros(N,Nf);
n = -(N-1)/2:(N-1)/2;

x_array = n*dH;
phi = linspace(0, pi, 180); % discrete Theta angle 
Nphi = length(phi);


BP = -0.14 - 0.57*cos(phi(end:-1:1)) + 0.57*cos(phi(end:-1:1)).^2 + 1.15*cos(phi(end:-1:1)).^3; 
BP = BP/max(BP);

figure();
plot(abs(BP));
ylabel('Gain');
xlabel('Incident angle');
title('reference beambattern for line array')
set(gcf,'color','w');
%%

FI = zeros(Nphi,Nf);
figure()

for i=1:Nf
  Br = zeros(N,1);  
  Rc = (f(i)*N*dH/c);
  %f(i) = Rc*c/N/dH;
  %x_points = Rc*cos(phi(end:-1:1));
  boudary  = find(abs(n)<= Rc  );
  phiS = acos(n(boudary)/Rc);
  Br(boudary)= -0.14 - 0.57*cos(phiS) + 0.57*cos(phiS).^2 + 1.15*cos(phiS).^3;
  Br(boudary) = abs(Br(boudary));%/max(Br(boudary));

  
  idxBr = boudary(end);
  %Br(1:boudary(1)-1)=Br(boudary(1)); 
  Br(idxBr+1:end)=1;
  idxBr = boudary(1);
  Br(1:idxBr)=Br(idxBr);
  plot(abs(Br)) 
  temp = fftshift(Br(end:-1:1));
  h = fftshift(ifft(temp(end:-1:1))); 

  x_map(:,i) = h;
  
  
  %hold on
  pause(0.01)
  
  %beam plot
  beta = 2*pi*f(i)/c;             % wave number
  D = exp(1j*beta*x_array(ones(1,Nphi),:).*cos(phi(ones(N,1),:))');
  FI(:,i) = (D*(h));     
end
%%
pos = [0.5 0.5 0.4 0.4];
figure('numbertitle','off','name','beam pattern of full array',...
      'Units','normal','Position',pos);
surf(180/pi*phi,f,abs(FI)');
axis tight
set(gca,'XTick',[0 45 90 135 180]);
view([25,50]);
xlabel('Azimuth \phi in deg');
ylabel('frequency in Hz');
zlabel('Magnitude');
grid on
set(gcf,'color','w');