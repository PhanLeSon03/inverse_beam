clear;
clc;
set(gcf,'color','w');
c = 343.2;
ft = 48000;
dH = 0.01;
NFIR = 256;
fStep = (ft/2)/(NFIR/2+1);

fLow = 200;
fHigh = 20000;

kTUH = floor(fLow/fStep);
kTOH = floor(fHigh/fStep);

f = fStep*(kTUH:kTOH);       % frequencies used to compute W
% Beamformer
N = 200;
n1 = -N/2:1:(N-2)/2;
n2 = n1;
[n1,n2] = meshgrid(n1,n2);
aux=sqrt(n1.^2+n2.^2);
d=dH;
thetaC = pi/8;
thetaC1 = pi/16;
W_2D = zeros(N, N,kTOH-kTUH+1);

dBmax = 40;

B_REF = zeros(100,120);
x = linspace(0,5,100);
t = linspace(0, pi/2, 100);
p = linspace(-pi, pi, 120);
%y=dBmax*sinc(x);
y =sin(pi*x)./(pi*x);
y(1) = 1;
B_REF(:,1) = abs(y);
for i=2:120
    B_REF(:,i) = B_REF(:,1);
end

% plot(B_REF(1:90,1))
% ylabel('gain');
% xlabel('elevation angle');
% axis tight
% set(gcf,'color','w');
%%
% 3D plot
B_REF = max(0,20*log10(B_REF+eps)+dBmax);
Xc = B_REF .* (sin(t')*cos(p));
Yc = B_REF .* (sin(t')*sin(p));
Zc = B_REF .* (cos(t')*ones(1,length(p)));

mesh(Xc,Yc,Zc,B_REF);
title('Expected Beam-pattern: dB');

alph = x(2)./t(2);
   
figure(1)
Sphere = ones(100,120);
B_REF = zeros(100,120);
%calculate the weight for broadband beamformer ============================
for k=kTUH:kTOH

    Rc = k*fStep*N*dH/c;
    boudary  = find(aux> Rc );
    temp = aux/Rc;    
    temp1 = alph*pi*asin(temp);
    temp2 = abs(sin(temp1)./(temp1));
    temp2(boudary)=0;
    Hd = (temp2);
    Hd(N/2+1,N/2+1)=1;

    % 3D plot
%     if k > kTOH -30    
%        
%         Sphere = ones(100,120)*Rc;
%         Line = Hd(101,101:101+floor(Rc));
%         B_REF(:,1) = interp1(1:floor(Rc)+1,Line,1:100);
%         for i=2:120
%             B_REF(:,i) = B_REF(:,1);
%         end
%         Xc = Sphere .* (sin(t')*cos(p));
%         Yc = Sphere .* (sin(t')*sin(p));
%         Zc = Sphere .* (cos(t')*ones(1,length(p)));
% 
%         mesh(Xc,Yc,Zc,B_REF,'FaceLighting','gouraud','LineWidth',0.1);
%         title('Spherical Coordinate at difference Radius');
%         alpha 0.01
%         hold on
%     end    
%    colormap(parula(64))
    mesh(Hd);
    hold on
    axis off
    title(k*fStep)

    h =fftshift(ifft2(rot90(fftshift(rot90(Hd,2)),2)));

    W_2D(:,:,k-kTUH+1) = h;

    pause(0.05);
end

W = ones(N*N,kTOH-kTUH+1);
for iMic=1:N
   for jMic=1:N
      W((iMic-1)*N+jMic,:) = (W_2D(iMic,jMic,:));
   end
end

%verify the beam pattern ==================================================

% Full array
mics_ref = zeros(N*N,2);
for yPos = 1: N
  for xPos = 1:N
      mics_ref(xPos + N*(yPos-1), 1) = (xPos)*dH-(N/2)*dH;   
      mics_ref(xPos + N*(yPos-1), 2) = (yPos)*dH-(N/2)*dH;
  end
end


%%%%%%%%%% plot beam response
Fsh = ft/2;                  % half of sampling frequency
Nh = NFIR/2+1;   
klow = round(fLow/Fsh*Nh);     % low index
kup = round(fHigh/Fsh*Nh);      % high index
fp = [2500 4000 8000 16000];
kf = round(fp/Fsh*Nh);
%fp = Fsh/Nh*kf;              % frequencies rounded to FFT grid
kf = kf-klow+1;              % index used in W matrix corresponding to fp 
pos = [0.045 0.045 0.45 0.45];
% for k = 1:length(fp)
%     figure('numbertitle','off','name','Fixed Array radiation pattern (dB)',...
%                   'Units','normal','Position',pos);
% 
%     %Plotting beam-pattern for 1st looking direction 
%     [R,t,p] = array_pattern_fft(mics_ref,W,fp(k),kf(k)); 
%     R = R/max(R(:));
% 
%     RdB = max(0,10*log10(R+eps)+dBmax);
%     %RdB = max(0,sqrt(R));
%     RB_90 = RdB(end,:);
%     BW  = find(RB_90 > dBmax -3 ); 
%     %polar(p,RdB(end,:)); 
% 
%     % 3D plot
%     Xc = RdB .* (sin(t')*cos(p));
%     Yc = RdB .* (sin(t')*sin(p));
%     Zc = RdB .* (cos(t')*ones(1,length(p)));
% 
%     mesh(Xc,Yc,Zc,RdB);
%        
%     pos(1) = pos(1) + 0.3;
%     
%     ght = title(sprintf('f = %3.2f Hz, unit dB',fp(k)));
%     colormap jet;
% end   

[goc, BS] = Directivity(mics_ref,W,f);
BS_dB = max(0,10*log10(BS+eps)+dBmax);
[Fxx,Fyy] = meshgrid(goc/pi*180,f);

figure('numbertitle','off','name','Fixed Array radiation pattern (dB)',...
                  'Units','normal','Position',[0.1 0.1 0.5 0.5]);
surf(Fxx,Fyy,sqrt(BS));
xlabel('elevator');
ylabel('frequency');
set(gcf,'color','w');
