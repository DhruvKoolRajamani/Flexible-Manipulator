function drawmanip(y,m,l1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    l2 = l1/2;
    l3 = l2/2;
    l1_cg = 2*l1/3;
    l2_cg = 2*l2/3;
    l3_cg = 2*l3/3;

    th1 = y(1);
    th2 = y(2);
    th3 = y(3);
    
    th12 = th1 + th2;
    th123 = th12 + th3;
    
    L_mw = 2*sqrt(m);
    
    x = 0;
    y = 0;
    
    px1 = l1*sin(th1);
    py1 = -l1*cos(th1);
    
    px2 = l1_cg*sin(th1) + l2*sin(th12);
    py2 = -l1_cg*cos(th1) - l2*cos(th12);
    
    px3 = l1_cg*sin(th1) + l2_cg*sin(th12) + l3*sin(th123);
    py3 = -l1_cg*cos(th1) - l2_cg*cos(th12) - l3*cos(th123);
    
    plot([-10 10],[0 0],'b','LineWidth',2)
    hold on
    
    %plot([x px1 px2 px3],[y py1 py2 py3],'r','LineWidth',L_mw)
    plot([x px1],[y py1],'r','LineWidth',L_mw);
    plot([px1 px2],[py1 py2],'b','LineWidth',L_mw);
    plot([px2 px3],[py2 py3],'g','LineWidth',L_mw);
    
%     xlim([-5 5]);
     ylim([-8 8]);
%    set(gca,'Color','k','XColor','w','YColor','w')
%    set(gcf,'Position',[10 900 800 400])
%    set(gcf,'Color','k')
    
    drawnow
    hold off
end

