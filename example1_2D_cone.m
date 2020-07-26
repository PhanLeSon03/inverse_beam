
c = 343.2;
ft = 48000;
dH = 0.015;
NFIR = 256;
fStep = (ft/2)/(NFIR/2+1);

fLow = 1000;
fHigh = 16000;

kTUH = floor(fLow/fStep);
kTOH = floor(fHigh/fStep);

f = fStep*(kTUH:kTOH);       % frequencies used to compute W
% Beamformer
N = 25;
n1 = -(N-1)/2:1:(N-1)/2;
n2 = n1;
[n1,n2] = meshgrid(n1,n2);
aux=sqrt(n1.^2+n2.^2);
d=dH;
thetaC = pi/12;
W_2D = zeros(N, N,kTOH-kTUH+1);



LP = fir1(4,0.4,'low');
LPW = LP'*LP;



t = linspace(0, pi/2, 50);
p = linspace(-pi, pi, 50);
B_REF = zeros(50,50);

y = ones(50,1);
y(floor(sin(thetaC)*50):end) = 0;
B_REF(:,1) = abs(y);
for i=2:50
    B_REF(:,i) = B_REF(:,1);
end
dBmax=40;
% 3D plot
figure()
B_REF = max(0,20*log10(B_REF+eps)+dBmax);
Xc = B_REF .* (sin(t')*cos(p));
Yc = B_REF .* (sin(t')*sin(p));
Zc = B_REF .* (cos(t')*ones(1,length(p)));

mesh(Xc,Yc,Zc,B_REF);
    hold on;
    plot3(dBmax*cos(p),dBmax*sin(p),zeros(length(p),1),'b--');
    plot3(dBmax*sin(t),zeros(length(t),1),dBmax*cos(t),'b--');
    plot3(zeros(length(t),1),dBmax*sin(t),dBmax*cos(t),'b--');
title('Reference beam pattern: dB');
set(gcf,'color','w');


Sphere = ones(50,50);
figure()
%calculate the weight for broadband beamformer ============================
for k=kTUH:kTOH

    R = k*fStep*N*dH/c;
    Hd = ones(N);%cos(aux.*pi/18).^2; 
    Hd(aux > R*sin(thetaC)) = 0.0;
    
    %filter to remove overshoot of invert FFT of cylinder
    Hd = filter2(LPW,Hd,'same');
    Hd(floor(N/2)+1,floor(N/2)+1)=1;

     mesh(n1,n2,Hd)
     axis off


    h =fftshift(ifft2(rot90(fftshift(rot90(Hd,2)),2)));

    W_2D(:,:,k-kTUH+1) = h;
 
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
fp = [3000 4000 8000 16000];
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

    % 3D plot
    Xc = RdB .* (sin(t')*cos(p));
    Yc = RdB .* (sin(t')*sin(p));
    Zc = RdB .* (cos(t')*ones(1,length(p)));

    mesh(Xc,Yc,Zc,RdB);
   
    
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
    colormap jet;
    set(gcf,'color','w');
end   

[goc, BS] = Directivity(mics_ref,W,f);
BS_dB = max(0,10*log10(BS+eps)+dBmax);
[Fxx,Fyy] = meshgrid(goc/pi*180,f);

figure('numbertitle','off','name','Fixed Array radiation pattern (dB)',...
                  'Units','normal','Position',[0.1 0.1 0.5 0.5]);
surf(Fxx,Fyy,BS);
xlabel('elevation angle \theta in deg');
ylabel('frequency in Hz');
zlabel('Magnitude');
%colormap jet;
set(gcf,'color','w');