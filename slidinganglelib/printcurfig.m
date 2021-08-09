% currentfig1=gcf;
% fnam=currentfig1.FileName;
% 
% saveas(gcf,'Temp.fig')

currentfig=gcf;
currentfig.InvertHardcopy='on';
currentfig.Renderer='painters';
set(gcf,'Renderer','painters');
set(gcf, 'DefaultAxesFontName', 'Arial');
% set(gcf,'Renderer','opengl');

% A = print2array(currentfig);
% B = sum(A,3);
% loc1=find(sum(B,1)~=0);

% print(currentfig, [currentfig.FileName(1:end-4) '.jpg'], '-djpeg')
% a=sum(imread('x.jpg'),3);
% loc1=find((sum(a,1)./size(a,1))~=3*255, 1 );
% loc2=find((sum(a,1)./size(a,1))~=3*255, 1, 'last' );
% loc3=find((sum(a,2)./size(a,2))~=3*255, 1 );
% loc4=find((sum(a,2)./size(a,2))~=3*255, 1, 'last' );
% 
% wd=loc2-loc1;
% hg=loc4-loc3;
% wdadj=500/hg*wd;
% ax.Position=[0 0 (loc2-loc1)/size(a,2) (loc4-loc3)/size(a,1)];

% currentfig.Position(1)=0;
% currentfig.Position(2)=0;
% set(gcf,'PaperUnits','points')
% set(gcf,'PaperSize',[560 420])
% set(gcf,'PaperPositionMode','auto')
% sza=currentfig.PaperPosition;
% set(gcf,'PaperUnits','points')
% set(gcf,'PaperSize',[sza(4) sza(3)])
% set(gcf,'Resize','off')


% set(gcf,'Position',[0 0 wdadj 500])
% ax1=currentfig.CurrentAxes;
% ax1.Units='points';
% ax1.Position(1)=ax1.Position(1)+15;


% set(gcf,'Position',[0 0 wdadj 500])
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];

% print(resize,'bestfit'R)
% print(currentfig, [currentfig.FileName(1:end-4) '.jpg'], '-djpeg')
% print(currentfig, [currentfig.FileName(1:end-4) '.bmp'], '-dbmp', '-r300')
% print(currentfig, [currentfig.FileName(1:end-4) '.eps'], '-depsc2', '-r300')
% print(currentfig, [currentfig.FileName(1:end-4) '-2.emf'], '-dmeta', '-r300')
% print(currentfig, [currentfig.FileName(1:end-4) '.png'], '-dpng')
exportgraphics(currentfig, [currentfig.FileName(1:end-4) '.pdf'],'ContentType','vector','Resolution',300)
exportgraphics(currentfig, [currentfig.FileName(1:end-4) '.emf'],'ContentType','vector','Resolution',300)
exportgraphics(currentfig, [currentfig.FileName(1:end-4) '.jpg'],'Resolution',300)
exportgraphics(currentfig, [currentfig.FileName(1:end-4) '.png'],'Resolution',300)
% delete Temp.fig
% close all
% 
% openfig(fnam);


% export_fig FIGEXP -dmeta -r300 -transparent
% movefile('FIGEXP.png',[currentfig.FileName(1:end-4) '-2.emf'],'f');