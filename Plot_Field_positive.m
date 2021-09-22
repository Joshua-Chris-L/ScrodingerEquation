function DWs = Dws( )
clear all ; close all ;

nx = 81 ; ny  = 81 ;

dx = 2.0e-9 ;

nn = 1e9 ;

Nmax = 10 ;

RN = 1/(4*pi) ;

RR = 1/(pi*40*40) ; 

vidfile = VideoWriter('testmovie.mp4','MPEG-4');
vidfile.FrameRate = 3 ;
open(vidfile);

for j = 1 : Nmax 
    TT(j) = (j-1)*500.0*5.0e-14*nn ;  
    chi = num2str(j) ;
    MST = importdata(strcat(strcat('mag_Field_', chi),'.txt')) ; 
     
    %MST = importdata('mag_26.txt') ;
    
    ind = MST.data(:,1:2); 
    sss = MST.data(:,3:8);
    qqq = MST.data(:,9:10);

    II = reshape(ind(:,1),nx,ny)';    
    JJ = reshape(ind(:,2),nx,ny)';
    
    magx_a = reshape(sss(:,1),nx,ny)';  
    magy_a = reshape(sss(:,2),nx,ny)';
    magz_a = reshape(sss(:,3),nx,ny)'; 

    magx_b = reshape(sss(:,4),nx,ny)';  
    magy_b = reshape(sss(:,5),nx,ny)';
    magz_b = reshape(sss(:,6),nx,ny)'; 
    
    qqq_a = reshape(qqq(:,1),nx,ny)';
    qqq_b = reshape(qqq(:,2),nx,ny)';
    
    QQa(j) = 0.0 ;
    QQb(j) = 0.0 ;
    Rxa(j) = 0.0 ;
    Rya(j) = 0.0 ;
    Rxb(j) = 0.0 ;
    Ryb(j) = 0.0 ;
    
    mmagx = 0.5*(magx_a + magx_b) ;
    nmagx = 0.5*(magx_a - magx_b) ;
    
    mmagy = 0.5*(magy_a + magy_b) ;
    nmagy = 0.5*(magy_a - magy_b) ;
    
    mmagz = 0.5*(magz_a + magz_b) ;
    nmagz = 0.5*(magz_a - magz_b) ;
    
    
    for ii = 1 : nx 
        for jj = 1 : ny
            xx = II(ii,jj) ;
            yy = JJ(ii,jj)  ;
            rr = sqrt(xx*xx +yy*yy) ;
            QQa(j) = QQa(j) + RN*qqq_a(ii,jj)*dx*dx ; 
            QQb(j) = QQb(j) + RN*qqq_b(ii,jj)*dx*dx ; 
            
            Rxa(j) = Rxa(j) + RN*xx*dx*qqq_a(ii,jj)*dx*dx ; 
            Rxb(j) = Rxb(j) + RN*xx*dx*qqq_b(ii,jj)*dx*dx ;
            
            Rya(j) = Rya(j) + RN*yy*dx*qqq_a(ii,jj)*dx*dx ; 
            Ryb(j) = Ryb(j) + RN*yy*dx*qqq_b(ii,jj)*dx*dx ;
            
        end 
    end 
    
    
    figure
    pcolor(dx*II,dx*JJ,nmagz); shading interp; axis off ; colorbar ; hold on;
    %quiver(dx*II,dx*JJ,mmagx, mmagy, 'k'); shading interp; axis off ;
    set(gca,'Fontsize',14,'fontweight','b');
    set(gcf,'Renderer','opengl'); title('m_z')
         
         
         

 

    F(j) = getframe(gcf);
    writeVideo(vidfile, F(j));
    
    
    
end 

close(vidfile)


ti = 2 ; tf = 4 ;

vx_a = nn*(Rxa(tf)/QQa(tf) - Rxa(ti)/QQa(ti))/(TT(tf) - TT(ti))
vx_b = nn*(Rxb(tf)/QQb(tf) - Rxb(ti)/QQb(ti))/(TT(tf) - TT(ti))


vy_a = nn*(Rya(tf)/QQa(tf) - Rya(ti)/QQa(ti))/(TT(tf) - TT(ti))
vy_b = nn*(Ryb(tf)/QQb(tf) - Ryb(ti)/QQb(ti))/(TT(tf) - TT(ti))





figure 
subplot(321) ; plot(TT, QQa, '-or','linewidth', 2.5,'MarkerSize',4); hold on ; 
set(gca,'Fontsize',32,'linewidth', 2.0, 'fontweight','n');   
xlabel('Time (ns)','Fontsize',32,'fontweight','n'); 
ylabel('Q_a','Fontsize',32,'fontweight','n');

subplot(322) ; plot(TT, QQb, '-ob','linewidth', 2.5,'MarkerSize',4); hold on ; 
set(gca,'Fontsize',32,'linewidth', 2.0, 'fontweight','n');   
xlabel('Time (ns)','Fontsize',32,'fontweight','n'); 
ylabel('Q_b','Fontsize',32,'fontweight','n');

subplot(323) ; plot(TT, nn*Rxa./QQa, '-or','linewidth', 2.5,'MarkerSize',4); hold on ; 
set(gca,'Fontsize',32,'linewidth', 2.0, 'fontweight','n');   
xlabel('Time (ns)','Fontsize',32,'fontweight','n'); 
ylabel('R_a^x(nm)','Fontsize',32,'fontweight','n');

subplot(324) ; plot(TT, nn*Rxb./QQb, '-ob','linewidth', 2.5,'MarkerSize',4); hold on ; 
set(gca,'Fontsize',32,'linewidth', 2.0, 'fontweight','n');   
xlabel('Time (ns)','Fontsize',32,'fontweight','n'); 
ylabel('R_b^x (nm)','Fontsize',32,'fontweight','n');

subplot(325) ; plot(TT, nn*Rya./QQa, '-or','linewidth', 2.5,'MarkerSize',4); hold on ; 
set(gca,'Fontsize',32,'linewidth', 2.0, 'fontweight','n');   
xlabel('Time (ns)','Fontsize',32,'fontweight','n'); 
ylabel('R_a^y (nm)','Fontsize',32,'fontweight','n');

subplot(326) ; plot(TT, nn*Ryb./QQb, '-ob','linewidth', 2.5,'MarkerSize',4); 
set(gca,'Fontsize',32,'linewidth', 2.0, 'fontweight','n');   
xlabel('Time (ns)','Fontsize',32,'fontweight','n'); 
ylabel('R_b^y (nm)','Fontsize',32,'fontweight','n');






end

