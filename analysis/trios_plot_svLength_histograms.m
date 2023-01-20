function trios_plot_svLength_histograms()
    SV_LENGTHS=[500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000];
    X_TICK_LABELS={'P1','P2','C','T','0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'};
    
    pbsv=load('scratch_cunial-read-length-distribution_2842314_pbsv_svLengths.histogram');
    sniffles1=load('scratch_cunial-read-length-distribution_2842314_sniffles1_svLengths.histogram');
    sniffles2=load('scratch_cunial-read-length-distribution_2842314_sniffles2_svLengths.histogram');
    
    % Total number of SVs
    [nRows,nColumns]=size(pbsv);
    totals_pbsv=zeros(nRows,1);
    totals_pbsv(1)=sum(pbsv(3,:));
    totals_pbsv(2)=sum(pbsv(4,:));
    totals_pbsv(3)=sum(pbsv(2,:));
    totals_pbsv(4)=sum(pbsv(1,:));
    for i=[5:nRows]
        totals_pbsv(i)=sum(pbsv(i,:));
    endfor
    [nRows,nColumns]=size(sniffles1);
    totals_sniffles1=zeros(nRows,1);
    totals_sniffles1(1)=sum(sniffles1(3,:));
    totals_sniffles1(2)=sum(sniffles1(4,:));
    totals_sniffles1(3)=sum(sniffles1(2,:));
    totals_sniffles1(4)=sum(sniffles1(1,:));
    for i=[5:nRows]
        totals_sniffles1(i)=sum(sniffles1(i,:));
    endfor
    [nRows,nColumns]=size(sniffles2);
    totals_sniffles2=zeros(nRows,1);
    totals_sniffles2(1)=sum(sniffles2(3,:));
    totals_sniffles2(2)=sum(sniffles2(4,:));
    totals_sniffles2(3)=sum(sniffles2(2,:));
    totals_sniffles2(4)=sum(sniffles2(1,:));
    for i=[5:nRows]
        totals_sniffles2(i)=sum(sniffles2(i,:));
    endfor
    figure(3);
    subplot(1,3,1); bar(totals_pbsv);
    xticklabels(X_TICK_LABELS);
    axis square; grid on; ylabel('Total number of SVs');
    title('pbsv','fontsize', 18);
    subplot(1,3,2); bar(totals_sniffles1);
    xticklabels(X_TICK_LABELS);
    axis square; grid on; ylabel('Total number of SVs');
    title('sniffles1','fontsize', 18);
    subplot(1,3,3); bar(totals_sniffles2);
    xticklabels(X_TICK_LABELS);
    axis square; grid on; ylabel('Total number of SVs');
    title('sniffles2','fontsize', 18);
    
    
    % Length histograms
    figure(1); 
    subplot(1,3,1); hold on;
    plot(SV_LENGTHS,pbsv(2,:),'.--b');
    plot(SV_LENGTHS,pbsv(3,:),'.--k');
    plot(SV_LENGTHS,pbsv(4,:),'.--k');
    plot(SV_LENGTHS,pbsv(1,:),'.-r');
    axis square; grid on; xlabel('SV length'); ylabel('Number of SVs');
    legend('child','parent 1','parent 2','baseline');
    title('pbsv: SVs in the baseline and in C,M,F','fontsize', 18);
    set(gca, 'fontsize', 18);
    
    subplot(1,3,2); hold on;
    plot(SV_LENGTHS,sniffles1(2,:),'.--b');
    plot(SV_LENGTHS,sniffles1(3,:),'.--k');
    plot(SV_LENGTHS,sniffles1(4,:),'.--k');
    plot(SV_LENGTHS,sniffles1(1,:),'.-r');
    axis square; grid on; xlabel('SV length'); ylabel('Number of SVs');
    legend('child','parent 1','parent 2','baseline');
    title('sniffles1: SVs in the baseline and in C,M,F','fontsize', 18);
    set(gca, 'fontsize', 18);
    
    subplot(1,3,3); hold on;
    plot(SV_LENGTHS,sniffles2(2,:),'.--b');
    plot(SV_LENGTHS,sniffles2(3,:),'.--k');
    plot(SV_LENGTHS,sniffles2(4,:),'.--k');
    plot(SV_LENGTHS,sniffles2(1,:),'.-r');
    axis square; grid on; xlabel('SV length'); ylabel('Number of SVs');
    legend('child','parent 1','parent 2','baseline');
    title('sniffles2: SVs in the baseline and in C,M,F','fontsize', 18);
    set(gca, 'fontsize', 18);

    
    figure(2); 
    subplot(1,3,1); hold on;
    plot(SV_LENGTHS,pbsv(1,:),'.-r');
    [nRows,nColumns]=size(pbsv);
    for i=[5:nRows]
        plot(SV_LENGTHS,pbsv(i,:),'.-k');
    endfor
    axis square; grid on; xlabel('SV length'); ylabel('Number of SVs');
    legend('baseline');
    title('pbsv','fontsize', 18);
    set(gca, 'fontsize', 18);
    
    subplot(1,3,2); hold on;
    plot(SV_LENGTHS,sniffles1(1,:),'.-r');
    [nRows,nColumns]=size(sniffles1);
    for i=[5:nRows]
        plot(SV_LENGTHS,sniffles1(i,:),'.-k');
    endfor
    axis square; grid on; xlabel('SV length'); ylabel('Number of SVs');
    legend('baseline');
    title('sniffles1','fontsize', 18);
    set(gca, 'fontsize', 18);
    
    subplot(1,3,3); hold on;
    plot(SV_LENGTHS,sniffles2(1,:),'.-r');
    [nRows,nColumns]=size(sniffles2);
    for i=[5:nRows]
        plot(SV_LENGTHS,sniffles2(i,:),'.-k');
    endfor
    axis square; grid on; xlabel('SV length'); ylabel('Number of SVs');
    legend('baseline');
    title('sniffles2','fontsize', 18);
    set(gca, 'fontsize', 18);
endfunction
