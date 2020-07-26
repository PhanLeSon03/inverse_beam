set(gcf,'color','w');
c = 343.2;
ft = 32000;
dH = 0.015;
NFIR = 256;
fStep = (ft/2)/(NFIR/2+1);

fLow = 250;
fHigh = 16000;

kTUH = floor(fLow/fStep);
kTOH = floor(fHigh/fStep);

f = fStep*(kTUH:kTOH);       % frequencies used to compute W
% Beamformer
N = 100;
n1 = -(N-1)/2:1:(N-1)/2;
n2 = n1;
[n1,n2] = meshgrid(n1,n2);
aux=sqrt(n1.^2+n2.^2);
d=dH;
thetaC = pi/8;
thetaC1 = pi/16;
W_2D = zeros(N, N,kTOH-kTUH+1);

LP = fir1(4,0.4,'low');
LPW = LP'*LP;
t = pi/180*(-180:2.5:180);%

t = linspace(0, pi/2, 100);
p = linspace(-pi, pi, 100);
B_REF = zeros(100,100);
y = ones(100,1);
y(floor(sin(thetaC1)*100):floor(sin(thetaC)*100)) = 1/10;
y(floor(sin(thetaC)*100)+1:end) = 0;
B_REF(:,1) = abs(y);
for i=2:100
    B_REF(:,i) = B_REF(:,1);
end
dBmax=40;
% 3D plot
B_REF = max(0,20*log10(B_REF+eps)+dBmax);
Xc = B_REF .* (sin(t')*cos(p));
Yc = B_REF .* (sin(t')*sin(p));
Zc = B_REF .* (cos(t')*ones(1,length(p)));
figure()
mesh(Xc,Yc,Zc,B_REF);
    hold on;
    plot3(dBmax*cos(p),dBmax*sin(p),zeros(length(p),1),'b--');
    plot3(dBmax*sin(t),zeros(length(t),1),dBmax*cos(t),'b--');
    plot3(zeros(length(t),1),dBmax*sin(t),dBmax*cos(t),'b--');
title('Reference beam pattern: dB');
set(gcf,'color','w');
%calculate the weight for broadband beamformer ============================ 
 figure()
for k=kTUH:kTOH
    w0=2*pi*ft*k/NFIR;
    k0= w0/c;
    R= k0*dH*sin(thetaC);

    
    Hd = ones(N);%cos(aux.*pi/18).^2; 
    %Hd((aux > R*alpha)|(aux < 0.5*R*alpha)) = 0.0;
    Hd((aux > k*fStep*N*dH*sin(thetaC1)/c)) = 0.1;
    Hd((aux > k*fStep*N*dH*sin(thetaC)/c)) = 0;
    Hd = filter2(LPW,Hd,'same');

     mesh(n1,n2,Hd)
     axis off

    h =fftshift(ifft2(rot90(fftshift(rot90(Hd,2)),2)));
    h=h/sum(sum(h));
    W_2D(:,:,k-kTUH+1) = h;
    set(gcf,'color','w');
 pause(0.05);
end
set(gcf,'color','w');
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
      mics_ref(xPos + N*(yPos-1), 1) = (xPos)*dH-12.0*dH;   
      mics_ref(xPos + N*(yPos-1), 2) = (yPos)*dH-12.0*dH;
  end
end


%%%%%%%%%% plot beam response
dBmax = 40;
Fsh = ft/2;                  % half of sampling frequency
Nh = NFIR/2+1;   
klow = round(fLow/Fsh*Nh);     % low index
kup = round(fHigh/Fsh*Nh);      % high index
fp = [2500 4000 8000 16000];
kf = round(fp/Fsh*Nh);
%fp = Fsh/Nh*kf;              % frequencies rounded to FFT grid
kf = kf-klow+1;              % index used in W matrix corresponding to fp 
pos = [0.045 0.045 0.45 0.45];
for k = 1:length(fp)
    figure('numbertitle','off','name','Fixed Array radiation pattern (dB)',...
                  'Units','normal','Position',pos);

    %Plotting beam-pattern for 1st looking direction 
    [R,t,p] = array_pattern_fft(mics_ref,W,fp(k),kf(k)); 
    R = R/max(R(:));

    RdB = max(0,10*log10(R+eps)+dBmax);
    RB_90 = RdB(end,:);
    BW  = find(RB_90 > dBmax -3 ); 
    %polar(p,RdB(end,:)); 

    % 3D plot
    Xc = RdB .* (sin(t')*cos(p));
    Yc = RdB .* (sin(t')*sin(p));
    Zc = RdB .* (cos(t')*ones(1,length(p)));
    %Zc = RdB .* (ones(length(t),1)*ones(1,length(p)));
    %Xc = dBmax*(sin(t')*cos(p));
    %Yc = dBmax*(sin(t')*sin(p));
    %Zc = RdB;
    mesh(Xc,Yc,Zc,RdB);
    colormap jet;
    
    axis([-dBmax dBmax -dBmax dBmax 0 dBmax]);
    view([25,50]);
    xlabel('x');
    ylabel('y');
    zlabel('z');

    %draw radial grid lines

    hold on;
    plot3(dBmax*cos(p),dBmax*sin(p),zeros(length(p),1),'b--');
    plot3(dBmax*sin(t),zeros(length(t),1),dBmax*cos(t),'b--');
    plot3(zeros(length(t),1),dBmax*sin(t),dBmax*cos(t),'b--');
    hold off;
    pos(1) = pos(1) + 0.3;
    %pos(2) = pos(2) + 0.1;
    ght = title(sprintf('f = %3.2f Hz',fp(k)));
    set(gcf,'color','w');
end   

[goc, BS] = Directivity(mics_ref,W,f);
BS_dB = max(0,10*log10(BS+eps)+dBmax);
[Fxx,Fyy] = meshgrid(goc/pi*180,f);
%%
figure('numbertitle','off','name','Fixed Array radiation pattern (dB)',...
                  'Units','normal','Position',[0.1 0.1 0.5 0.5]);
%surf(Fxx,Fyy,BS_dB);
surf(Fxx,Fyy,sqrt(BS));
xlabel('elevation angle \theta in deg');
ylabel('frequency in Hz');
zlabel('Magnitude');
%zlabel('dB');
set(gcf,'color','w');
%colormap jet;