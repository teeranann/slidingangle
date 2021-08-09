% The calculation of the sliding angle by monitoring the image brightness. 
% Teeranan Nongnual
% teeranan.no@buu.ac.th
% Faculty of Science, Burapha University, THAILAND
% Version 2.1 (August 2021)
% https://github.com/teeranann/slidingangle

% Citations
% T.Nongnual, N.Damnong, S.Srimongkol, S.Phanrangsee,T.Benjalersyarnon, 
% Thai Journal of Mathematics, 2021, 19 (3).

% Grants
% NSTDA (#FDA-CO-2563-12485-TH)

%%
close all
clearvars -except pathnamm
currgenpath=genpath(pwd);
addpath(currgenpath)

%% 1. CALCULATIONS

%% 1.1 Config parameters
edgewidth=100;
cc = lines(25);
baselineright=0;
edgewidthright=350;
fitfrmse=0;

video=1;
saveprntfig=1;
edgelineplot=0;
figinsetrawimg=1;
stackimg=1;
baseslidingplot=0;
%%
edgedet

%% 1.2 Image brightness: Summation of color code in all pixels - Find peaks of image brightness
sumcolor=[];
for n=1:camframetot
    sumcolor(n)=sum(mov{1,n},'all');
end
sumcolordiff=diff(sumcolor);
sumcolorgrad=gradient(sumcolor);
[pkss1,loccolordiffpeaks1]=findpeaks(sumcolordiff,'threshold',0.1*camheight*camwidth,'MinPeakDistance',0.1*camfrate,'MinPeakHeight',0.1*max(sumcolordiff));
[pkss2,loccolordiffpeaks2]=findpeaks(sumcolorgrad,'threshold',0.1*camheight*camwidth,'MinPeakDistance',0.1*camfrate,'MinPeakHeight',0.1*max(sumcolorgrad));
[loccolordiffpeaks,locpkun]=uniquetol([loccolordiffpeaks1 loccolordiffpeaks2],0.03);
pkss=([pkss1 pkss2]);
pkss=pkss(locpkun);
loccolordiffmax=max(loccolordiffpeaks);

%% 1.3 Crop baseline and sliding line, fit SL/BL to get tilt angle
baseline=[EdgeCellSALR{1,1}.x EdgeCellSALR{1,1}.y];
cropbaseline=[];
if baselineright
    fprintf(2,'Baseline cropped on right.!\n')
end
for n=1:length(baseline)
    if baselineright
        if baseline(n,1)>=camwidth-edgewidthright || baseline(n,2)<=edgewidth
            cropbaseline(end+1,1:2)=baseline(n,1:2);
        end
    else
        if baseline(n,1)<=edgewidth || baseline(n,2)<=edgewidth
            cropbaseline(end+1,1:2)=baseline(n,1:2);
        end
    end
end
slidingline=[EdgeCellSALR{loccolordiffmax+1,1}.x EdgeCellSALR{loccolordiffmax+1,1}.y];
cropslidingline=[];
for n=1:length(slidingline)
    if slidingline(n,1)<=edgewidth || slidingline(n,2)<=edgewidth
        cropslidingline(end+1,1:2)=slidingline(n,1:2);
    end
end
[blfit,blfitcoeff]=polyfit(cropbaseline(:,1),cropbaseline(:,2),1);
blfitrsqr = 1 - (blfitcoeff.normr/norm(cropbaseline(:,2) - mean(cropbaseline(:,2))))^2;
[slfit,slfitcoeff]=polyfit(cropslidingline(:,1),cropslidingline(:,2),1);
slfitrsqr = 1 - (slfitcoeff.normr/norm(cropslidingline(:,2) - mean(cropslidingline(:,2))))^2;
slidinganglee= (atan((slfit(1)-blfit(1))/(1+blfit(1)*slfit(1)))) *180 /pi;

% %% Method [1]
% fprintf('sliding time = %1.3f s\n',timest(loccolordiffmax+1));
% fprintf('(1) sliding angle = %1.2f deg, R^2(baseline linear)=%1.4f, R^2(sliding linear)=%1.4f\n',slidinganglee,blfitrsqr,slfitrsqr);

% %% Method [3]
sumcolorpx=sumcolor/camheight/camwidth;
tm1=round((camduration-4.5)*camfrate);
tm2=round((camduration-3)*camfrate);
tm3=round((camduration-1.2)*camfrate);
tm4=camframetot;
xx=timest(tm1:tm4);yy=sumcolorpx(tm1:tm4);
peakflat=islocalmax(smooth(yy), 'FlatSelection', 'first');
peakflatf=find((peakflat==1));
peakflat1=peakflatf(1)+tm1-1;
tm2=peakflat1-5;
tm3=peakflat1+5;
tm1=round((tm2-1.5*camfrate));
[line1,line1coeff]=polyfit(timest(tm1:tm2),sumcolorpx(tm1:tm2),1);
line1rsqr = 1 - (line1coeff.normr/norm(sumcolorpx(tm1:tm2) - mean(sumcolorpx(tm1:tm2))))^2;
[line2,line2coeff]=polyfit(timest(tm3:tm4),sumcolorpx(tm3:tm4),1);
line2rsqr = 1 - (line2coeff.normr/norm(sumcolorpx(tm3:tm4) - mean(sumcolorpx(tm3:tm4))))^2;
xintercept=((line2(2)-line1(2))/(line1(1)-line2(1)));
frameintercept=round(xintercept*camfrate);
rotorstepsz=90/((line2(2)-line1(2))/(line1(1)-line2(1)))/camfrate;

Tangle={};
 for n=1:size(EdgeCellSALR,1)
     Tangle{n,1}=[];
    for p=1:size(EdgeCellSALR{n,1}.x,1)
        if ( EdgeCellSALR{n,1}.x(p)<=edgewidth || EdgeCellSALR{n,1}.y(p)<=edgewidth ) && ~isnan(EdgeCellSALR{n,1}.x(p)) && ~isnan(EdgeCellSALR{n,1}.y(p)) 
           Tangle{n,1}=[Tangle{n,1}; EdgeCellSALR{n,1}.x(p) EdgeCellSALR{n,1}.y(p)];
        end
    end
 end
 for n=1:length(Tangle)
     if isempty(Tangle{n,1})
         [~,Tangle{n,1}(:,2)]=max(abs(diff(double(mov{1,n}(:,1:150))))); %%%%% x=1:150,y=loc(diff)
         Tangle{n,1}(:,1)=1:150;
     end
 end
if fitfrmse
    fitfrmstart=round((slidinganglee-10)/0.1);
    fitfrmend=round((slidinganglee+10)/0.1);
else
    fitfrmstart=1;
    fitfrmend=length(Tangle);
end
allslope_1=[];
allangle=[];
allintercept=zeros(camframetot-1,2);
for n=fitfrmstart:fitfrmend
    warning('off')    
    [fitX,fitXcoef] = fit(Tangle{n,1}(:,1),Tangle{n,1}(:,2),'A*x+B');%,'startpoint',[allslope_1(end,1) allslope_1(end,2)]);
    fitXv = struct(fitX); fitXcoefv = struct(fitXcoef);
    warning('on')
    allslope_1(end+1,1:2)=[fitXv.coeffValues{1,1} fitXv.coeffValues{1,2}]; %%%R2x = fitXcoefv.rsquare;
end
for p=1:length(allslope_1)
    allangle(end+1,1)= abs(atan((blfit(1)-allslope_1(p,1))/(1+blfit(1)*allslope_1(p,1)))) *180 /pi;
end
for p=1:length(allslope_1)-1
    allintercept(p,1)=(allslope_1(p+1,2)-allslope_1(p,2))/(allslope_1(p,1)-allslope_1(p+1,1));
    allintercept(p,2)=allslope_1(p,1)*allintercept(p,2)+allslope_1(p,2);
end

timeststdcurve=timest(fitfrmstart:fitfrmend);

p = polyfit(timeststdcurve(1,:)',allangle(:,1),1);   
poly_line = polyval(p,timeststdcurve(1,:));  

%% 1.4 Method [2]
allanglemin=[];
[~,loci]=min(abs(allangle(:,1)-0.1));
[~,locf]=min(abs((allangle)-85));
if loci>locf
    loci=1;
end
alltanglelin=allangle(loci:locf,1);
timelin=timeststdcurve(1,loci:locf);
sumcolorlin=sumcolor(1,loci:locf);
sumcolordifflin=sumcolordiff(1,loci:locf);
warning('off')
[fitX,fitXcoef] = fit(timelin(1,:)',alltanglelin(:,1),'A*x+B');
fitXv = struct(fitX); fitXcoefv = struct(fitXcoef);
warning('on')
slopelin=fitXv.coeffValues{1,1};
interceptlin=fitXv.coeffValues{1,2};
R2x = fitXcoefv.rsquare;
[p,pp] = polyfit(timelin(1,:)',alltanglelin(:,1),1); 
poly_line = polyval(p,timelin(1,:));  


%% 2. ILLUSTRATIONS
disp('OUTPUTS')

%% 2.1 Config parameters
% Frames to be plotted
frameplot=sort([1 loccolordiffpeaks loccolordiffpeaks+1]); %loccolordiffpeaks
% frameplot=[1 loccolordiffmax loccolordiffmax+1]; %loccolordiffpeaks
% frameplot=[1 loccolordiffmax-9 loccolordiffmax loccolordiffmax+1 round(camframetot*1/2) round(camframetot*3/4) camframetot];

%% 2.2 Summation of color pixels = Brightness
sumcolorpx=sumcolor/camheight/camwidth;
sumcolorpxdiff=sumcolordiff/camheight/camwidth;
if saveprntfig
    if ~exist([pathnamm,fnam(1:end-4)],'dir')
        mkdir([pathnamm,fnam(1:end-4)])
    end
end

%% 2.3 Baseline vs sliding line
if baseslidingplot
    figure;
    scatter(cropbaseline(1:5:end,1),cropbaseline(1:5:end,2),30,'o','LineWidth',1.2,'MarkerEdgeColor','b'); hold all;
    scatter(cropslidingline(1:5:end,1),cropslidingline(1:5:end,2),30,'o','LineWidth',1.2,'MarkerEdgeColor','k')
    plot(cropbaseline(:,1),polyval(blfit,cropbaseline(:,1)),'r','linewidth',1.2);
    plot(cropslidingline(:,1),polyval(slfit,cropslidingline(:,1)),'m','linewidth',1.2);
    axis ij;
    set(gca,'fontsize',14,'linewidth',1.2)
    box on
    axis square
end

%% 2.4 Calculate Tf. xintercept = final time that the brightness is flat and the rotor stops
figure;
scatter(timest(tm1:tm4),sumcolorpx(tm1:tm4),15,'o','LineWidth',0.8,'MarkerEdgeColor','b'); hold on
xinterceptf=round(xintercept*camfrate);
plot(timest(tm1:xinterceptf+10),polyval(line1,timest(tm1:xinterceptf+10)),'r-','linewidth',1.2); hold on
plot(timest(xinterceptf-10:tm4),polyval(line2,timest(xinterceptf-10:tm4)),'LineStyle','-','color',rgb('green'),'linewidth',1.2);
xlabel('{\itt} (s)')
ylabel('{\itB}')
ylim([round(min(sumcolorpx(tm1:tm4))-1) round(max(sumcolorpx(tm1:tm4))+1)])
set(gca,'fontsize',14,'linewidth',1.2)
box on
axis square
xlim([timest(tm1) timest(tm4)])
timsc=timest(tm1:tm4);
[~,fnumxintercept]=min(abs(timest-xintercept));
text((timest(tm4)-timest(tm1))/10*6.5+timest(tm1), round(min(sumcolorpx(tm1:tm4))-1)+0.3, sprintf('{\\itt_F} = %1.2f s', xintercept), 'Fontsize', 12,'color','r');
if saveprntfig
    saveas(gcf,[pathnamm,fnam(1:end-4),'\',fnam(1:end-4),'-tf.fig']);
    openfig([pathnamm,fnam(1:end-4),'\',fnam(1:end-4),'-tf.fig']);
    printcurfig;
    close(gcf);
end

%% 2.5 Baseline with cropped baseline
figure;
for n=1:30:fnumxintercept
    scatter(EdgeCellSALR{n,1}.x,EdgeCellSALR{n,1}.y,15,'.'); hold all;
    axis ij;
    set(gca,'fontsize',14,'linewidth',1.2)
    set(gca,'XTicklabel',[],'YTicklabel',[])
    set(gca,'Xtick',[],'ytick',[])
    box on
end
for n=1:30:fnumxintercept
    scatter(Tangle{n,1}(:,1),Tangle{n,1}(:,2),15,'k.'); hold all;
    set(gca,'fontsize',14,'linewidth',1.2)
    axis ij;
end
xlim([0 camwidth])
ylim([0 camheight])
plot([edgewidth edgewidth camwidth],[camheight edgewidth edgewidth],'b','linewidth',1.2)
set(gca,'XTicklabel',[],'YTicklabel',[])
set(gca,'Xtick',0:100:camwidth,'ytick',fliplr(camheight:-100:0),'tickdir','in')
box on
set(gca,'dataaspectratio',[1 1 1]) % square pixels
if saveprntfig
    saveas(gcf,[pathnamm,fnam(1:end-4),'\',fnam(1:end-4),'-edgeallf.fig']);
    openfig([pathnamm,fnam(1:end-4),'\',fnam(1:end-4),'-edgeallf.fig']);
    printcurfig;
    close(gcf);
end

%% 2.6 EDGE LINE IMAGE of selected images
if edgelineplot
    for nmov=frameplot
        figure;
        scatter(EdgeCellSALR{nmov,1}.x,EdgeCellSALR{nmov,1}.y,15,'.'); hold all;
        axis ij;
        xlim([0 camwidth])
        ylim([0 camheight])
        set(gca,'fontsize',14,'linewidth',1.2)
        set(gca,'XTicklabel',[],'YTicklabel',[])
        set(gca,'Xtick',0:100:camwidth,'ytick',fliplr(camheight:-100:0),'tickdir','in')
        box on
        set(gca,'dataaspectratio',[1 1 1]) % square pixels
        
        txfn=text(round(camwidth/10*7),40,['[#',sprintf('%1.0f/%1.0f]',nmov,camframetot)]);
        set(txfn,'color','k');
        set(txfn,'fontsize',14);
        
        set(gca,'dataaspectratio',[1 1 1]) % square pixels
        if saveprntfig
            saveas(gcf,[pathnamm,fnam(1:end-4),'\',fnam(1:end-4),'-edge-',num2str(nmov),'.fig']);
            openfig([pathnamm,fnam(1:end-4),'\',fnam(1:end-4),'-edge-',num2str(nmov),'.fig']);
            printcurfig;
            close(gcf);
        end
    end
end

%% 2.7 Rotor stepsize deg/frame
figure;
timeilnsize=size(timelin,2);
plot(timelin(1,1:round((timeilnsize-1)/30):timeilnsize),alltanglelin(1:round((timeilnsize-1)/30):timeilnsize,1),'bo','linewidth',1.2); hold all;
plot(timelin(1,:),poly_line,'r','linewidth',1.2);
set(gca,'fontsize',14,'linewidth',1.2,'tickdir','out','ticklength',[0.015,0.015])
set(gcf,'windowstyle','normal')
set(gcf,'resize','off')
axis square
xlabel('{\itt} (s)')
ylabel(['{\it\theta_T} (' ,char(0176), ')'])
xlim([0 timelin(end)])
ylim([0 alltanglelin(end)])
scatter(timest(loccolordiffmax+1),0,50,'|','LineWidth',1.2,'MarkerEdgeColor','r'); hold on
plot([timest(loccolordiffmax+1) timest(loccolordiffmax+1)],[0 timest(loccolordiffmax+1)*p(1)+p(2)],'--','LineWidth',1.2,'Color','r')
scatter(0,timest(loccolordiffmax+1)*p(1)+p(2),50,'_','LineWidth',1.2,'MarkerEdgeColor',[0, 0.5, 0]); hold on
plot([0 timest(loccolordiffmax+1)],[timest(loccolordiffmax+1)*p(1)+p(2) timest(loccolordiffmax+1)*p(1)+p(2)],'--','LineWidth',1.2,'Color',[0, 0.5, 0])
set(gca,'XTick',0:5:xintercept)
text(xintercept*3/5, 15, sprintf('{\\itt_S} = %1.2f s', timest(loccolordiffmax+1)), 'Fontsize', 12,'color','r');
text(xintercept*3/5, 10, [sprintf('{\\it\\theta_S} = %1.1f%',timest(loccolordiffmax+1)*p(1)+p(2)),char(0176)], 'Fontsize', 12,'color',[0, 0.5, 0],'VerticalAlignment', 'top');
text(1, 75, sprintf('{\\it\\theta_T} = %1.2f{\\itt} + (%1.2f), {\\itR}^2 = %1.4f\n{\\itt_{lag}} = %1.2f s', p(1),p(2),R2x,-p(2)/p(1)), 'Fontsize', 12,'color','k');
if saveprntfig
    saveas(gcf,[pathnamm,fnam(1:end-4),'\',fnam(1:end-4),'-tiltangle.fig']);
    openfig([pathnamm,fnam(1:end-4),'\',fnam(1:end-4),'-tiltangle.fig']);
    printcurfig;
    close(gcf);
end

%% 2.8 Brightness plot
figure;
yyaxis left
Ax1 =plot(timest(1,:),sumcolorpx);
set(Ax1,'color','b','linewidth',1.2);
ylabel('{\itB}') % left y-axis
set(gca,'fontsize',14,'linewidth',1.2,'tickdir','out','ticklength',[0.015,0.015],'YColor','b');
axis square
hstart=0;hend=255;
set(gca,'YLim',[hstart hend]);
set(gca,'YTick',[0 64 128 192 255]);
set(gca,'TickLabelInterpreter', 'tex');
% set(gca,'YTickLabel',{'_{(black)}0', 64, 128,192,'_{(white)}255'})
text(-1.5, 245, 'white', 'Fontsize', 10,'color','b','horizontalalignment','right');
text(-1.5, 10, 'black', 'Fontsize', 10,'color','b','horizontalalignment','right');

yyaxis right
Ax2 =plot(timest(1,2:end),sumcolorpxdiff); hold on
set(Ax2,'color',[0, 0.5, 0],'linewidth',1.2);
ylabel(['{\itB}',char(39)]) % right y-axis

xlabel('{\itt} (s)')
set(gca,'fontsize',14,'linewidth',1.2,'tickdir','out','ticklength',[0.015,0.015],'YColor',[0, 0.5, 0])
axis square

hhstart=-5;
hhend=10+10*round(max(sumcolorpxdiff)/10);
hhend=max([hhend 15]);
set(gca,'YLim',[hhstart hhend]);
set(gca,'YTick',hhstart:5:hhend)
scatter(timest(loccolordiffpeaks+1),sumcolorpxdiff(loccolordiffpeaks)+1,40,'v','LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor',rgb('orange')); hold on
scatter(timest(loccolordiffmax+1),sumcolorpxdiff(loccolordiffmax)+1,40,'v','LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor','r'); hold on
xlim([0 xintercept])
set(gca,'XTick',0:5:xintercept)
fprintf('Final time of rotation:\n tf = %1.3f s, frame#(t=tf) = %1.0f\n', xintercept,fnumxintercept)
text(timest(loccolordiffmax+1)+2,sumcolorpxdiff(loccolordiffmax)+1.8, sprintf('{\\itt_S} = %1.2f s', timest(loccolordiffmax+1)), 'Fontsize', 12,'color','r','verticalalignment','top',...
    'horizontalalignment','left');
if saveprntfig
    saveas(gcf,[pathnamm,fnam(1:end-4),'\',fnam(1:end-4),'-brightness.fig']);
    openfig([pathnamm,fnam(1:end-4),'\',fnam(1:end-4),'-brightness.fig']);
    printcurfig;
    close(gcf);
end

%% 2.9 RAW IMAGE of selected images
for nmov=frameplot
    figure
    imagesc(mov{nmov})
    set(gca,'XTicklabel',[],'YTicklabel',[])
    set(gca,'Xtick',0:100:camwidth,'ytick',fliplr(camheight:-100:0),'tickdir','in')
    box on
    set(gca,'dataaspectratio',[1 1 1]) % square pixels
    axis ij
    axis equal
    colormap gray
    hold on
    xlim([0 camwidth])
    ylim([0 camheight])
    txfn=text(round(camwidth/10*7.5),40,['[#',sprintf('%1.0f/%1.0f]',nmov,camframetot)]);
    set(txfn,'color','k');
    set(txfn,'fontsize',14);
    txfn=text(round(camwidth/10*7.5),80,sprintf('[{\\itt} = %2.2f s]',timest(nmov)));
    set(txfn,'color','k');
    set(txfn,'fontsize',14);
    txfn=text(round(camwidth/10*7.5),120,['{\it\theta_T} = ',sprintf('%1.2f%',allangle(nmov)),char(0176)]);
    set(txfn,'color','k');
    set(txfn,'fontsize',14);
    if nmov== loccolordiffmax+1
        txfn=text(round(camwidth/10*7.5),160,['{\it\theta_S} = ',sprintf('%1.1f%',timest(loccolordiffmax+1)*p(1)+p(2)),char(0176)]);
        set(txfn,'color','r');
        set(txfn,'fontsize',14);
    end
    box on
    set(gca,'linewidth',1.2)
    if figinsetrawimg
        axes('position',[.13 .69 camwidth/camheight*0.2 0.2],'box','on')
        scatter(EdgeCellSALR{nmov,1}.x,EdgeCellSALR{nmov,1}.y,15,'.'); hold all;
        axis ij;
        xlim([0 camwidth])
        ylim([0 camheight])
        set(gca,'fontsize',14,'linewidth',1.2)
        set(gca,'XTicklabel',[],'YTicklabel',[])
        set(gca,'Xtick',0:100:camwidth,'ytick',fliplr(camheight:-100:0),'tickdir','in','TickLength',[0.025 0.05])
        box on
        set(gca,'dataaspectratio',[1 1 1]) % square pixels
    end
    if saveprntfig
        saveas(gcf,[pathnamm,fnam(1:end-4),'\',fnam(1:end-4),'-im-',num2str(nmov),'.fig']);
        openfig([pathnamm,fnam(1:end-4),'\',fnam(1:end-4),'-im-',num2str(nmov),'.fig']);
        printcurfig;
        close(gcf);
    end
end
hold off

%% 2.10 Summary results
fprintf('Tilt angle vs time:\n TA = m t + c\n m = %1.4f deg/s, c = %1.4f deg, stepsize = %1.4f deg/frame, R^2 = %1.4f\n', p(1),p(2),p(1)/camfrate,R2x)
fprintf('Rotor lag time:\n tlag = %1.2f s (x-intercept)\n', -p(2)/p(1))
fprintf('Sliding time:\n ts = %1.3f s (frame# = %1.0f)\n', timest(loccolordiffmax+1), loccolordiffmax+1)
fprintf('Sliding angle:\n (SA = TA(t=ts)) = %1.1f deg\n', timest(loccolordiffmax+1)*p(1)+p(2))

%% 2.11 Video output
if video
    writerObj = VideoWriter([pathnamm,fnam(1:end-4),'\',fnam(1:end-4),'-slidingangle.avi']); % Name it.
    writerObj.FrameRate = 25; % How many frames per second.
    open(writerObj);
end
if video
    figure
    for cframe=1:min([camframetot, round((90-p(2))/(p(1)/camfrate))]) %frame=1:showevfr:numframes
        imagesc(mov{cframe});
        set(gca,'XTicklabel',[],'YTicklabel',[])
        set(gca,'Xtick',0:100:camwidth,'ytick',fliplr(camheight:-100:0),'tickdir','in')
        box on
        set(gca,'dataaspectratio',[1 1 1]) % square pixels
        
        axis ij
        axis equal
        hold on
        colormap gray
        xlim([1 camwidth]);
        ylim([1 camheight]);

        txfn=text(round(camwidth/10*7.5),40,['[#',sprintf('%1.0f/%1.0f]',cframe,camframetot)]);
        set(txfn,'color','k');
        set(txfn,'fontsize',14);
        txfn=text(round(camwidth/10*7.5),80,sprintf('[{\\itt} = %2.2f s]',timest(cframe)));
        set(txfn,'color','k');
        set(txfn,'fontsize',14);
        
        txfn=text(round(camwidth/10*7.5),120,['{\it\theta_T} = ',sprintf('%1.2f%',allangle(cframe)),char(0176)]);
        set(txfn,'color','k');
        set(txfn,'fontsize',14);
        txfn=text(round(camwidth/10*7.5),160,['{\it\theta_S} = ',sprintf('%1.1f%',timest(loccolordiffmax+1)*p(1)+p(2)),char(0176)]);

        set(txfn,'color','r');
        set(txfn,'fontsize',14);
        set(gca,'linewidth',1.2)
        if figinsetrawimg
            axes('position',[.13 .67 camwidth/camheight*0.2 0.2],'box','on')
            scatter(EdgeCellSALR{cframe,1}.x,EdgeCellSALR{cframe,1}.y,15,'.'); hold all;
            axis ij;
            xlim([0 camwidth])
            ylim([0 camheight])
            set(gca,'fontsize',14,'linewidth',1.2)
            set(gca,'XTicklabel',[],'YTicklabel',[])
            set(gca,'Xtick',0:100:camwidth,'ytick',fliplr(camheight:-100:0),'tickdir','in','TickLength',[0.025 0.05])
            box on
            set(gca,'dataaspectratio',[1 1 1]) % square pixels
        end
        drawnow

        fr = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, fr);
        clf   %del figure.
    end
    close(writerObj);
    close(gcf)
end

%% 2.12 Stack image
imgplotorder=sort(unique([1:round(loccolordiffmax/3):loccolordiffmax loccolordiffmax-1 loccolordiffmax loccolordiffmax+1 loccolordiffmax:(camframetot-loccolordiffmax)/6:camframetot]));
imgplotorder=round(imgplotorder(1:9));
rzxaxismov=0.85;%rzxaxismov=4/5;
if stackimg
    figure('position',[0 0 round(camwidth*rzxaxismov) camheight])
    tlo = tiledlayout(3,3,'TileSpacing','none','Padding','none');
    for i=1:9
        nexttile, imagesc(mov{imgplotorder(i)}(:,camwidth-round(camwidth*rzxaxismov)+1:end)); hold on
        box on
        set(gca,'dataaspectratio',[1 1 1]) % square pixels
        axis ij
        axis equal
        colormap gray
        xlim([0 round(camwidth*rzxaxismov)])
        ylim([0 camheight])
        set(gca,'linewidth',1.2)
        if i>3 && i<7
            set(gca,'xcolor','r','ycolor','r','linewidth',1.5)
            patch([0 round(camwidth*rzxaxismov) round(camwidth*rzxaxismov) 0],[0 0 camheight camheight],'r','facealpha',0.05,'edgealpha',0.05)
        else
            set(gca,'xcolor','k','ycolor','k')
        end
    end
    set(tlo.Children,'Xtick',[],'ytick',[],'xticklabel',[],'yticklabel',[]);
    if saveprntfig
        saveas(gcf,[pathnamm,fnam(1:end-4),'\',fnam(1:end-4),'-im-9box.fig']);
        openfig([pathnamm,fnam(1:end-4),'\',fnam(1:end-4),'-im-9box.fig']);
        printcurfig;
        close(gcf);
    end
end
rmpath(currgenpath)
disp(' ')