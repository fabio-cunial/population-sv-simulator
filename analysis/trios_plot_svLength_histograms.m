function trios_plot_svLength_histograms()
    SV_LENGTHS=[500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000];
    
    pbsv=load('scratch_cunial-read-length-distribution_2842314_pbsv_svLengths.histogram');
    sniffles1=load('scratch_cunial-read-length-distribution_2842314_sniffles1_svLengths.histogram');
    sniffles2=load('scratch_cunial-read-length-distribution_2842314_sniffles2_svLengths.histogram');
    
    figure(1); 
    subplot(1,3,1); hold on;
    semilogy(SV_LENGTHS,pbsv(2,:),'.--k');
    semilogy(SV_LENGTHS,pbsv(3,:),'.--k');
    semilogy(SV_LENGTHS,pbsv(4,:),'.--k');
    semilogy(SV_LENGTHS,pbsv(1,:),'.-r');
    axis square; grid on; xlabel('SV length'); ylabel('Number of SVs');
    legend('individual 1','individual 2','individual 3');
    title('pbsv: SVs in the baseline and in C,M,F','fontsize', 18);
    set(gca, 'fontsize', 18);
    
    subplot(1,3,2); hold on;
    semilogy(SV_LENGTHS,sniffles1(2,:),'.--k');
    semilogy(SV_LENGTHS,sniffles1(3,:),'.--k');
    semilogy(SV_LENGTHS,sniffles1(4,:),'.--k');
    semilogy(SV_LENGTHS,sniffles1(1,:),'.-r');
    axis square; grid on; xlabel('SV length'); ylabel('Number of SVs');
    legend('individual 1','individual 2','individual 3');
    title('sniffles1: SVs in the baseline and in C,M,F','fontsize', 18);
    set(gca, 'fontsize', 18);
    
    subplot(1,3,3); hold on;
    semilogy(SV_LENGTHS,sniffles2(2,:),'.--k');
    semilogy(SV_LENGTHS,sniffles2(3,:),'.--k');
    semilogy(SV_LENGTHS,sniffles2(4,:),'.--k');
    semilogy(SV_LENGTHS,sniffles2(1,:),'.-r');
    axis square; grid on; xlabel('SV length'); ylabel('Number of SVs');
    legend('individual 1','individual 2','individual 3');
    title('sniffles2: SVs in the baseline and in C,M,F','fontsize', 18);
    set(gca, 'fontsize', 18);

    
    figure(2); 
    subplot(1,3,1); hold on;
    semilogy(SV_LENGTHS,pbsv(1,:),'.-r');
    [nRows,nColumns]=size(pbsv);
    for i=[5:nRows]
        semilogy(SV_LENGTHS,pbsv(i,:),'.-k');
    endfor
    axis square; grid on; xlabel('SV length'); ylabel('Number of SVs');
    legend('baseline');
    title('pbsv','fontsize', 18);
    set(gca, 'fontsize', 18);
    
    subplot(1,3,2); hold on;
    semilogy(SV_LENGTHS,sniffles1(1,:),'.-r');
    [nRows,nColumns]=size(sniffles1);
    for i=[5:nRows]
        semilogy(SV_LENGTHS,sniffles1(i,:),'.-k');
    endfor
    axis square; grid on; xlabel('SV length'); ylabel('Number of SVs');
    legend('baseline');
    title('sniffles1','fontsize', 18);
    set(gca, 'fontsize', 18);
    
    subplot(1,3,3); hold on;
    semilogy(SV_LENGTHS,sniffles2(1,:),'.-r');
    [nRows,nColumns]=size(sniffles2);
    for i=[5:nRows]
        semilogy(SV_LENGTHS,sniffles2(i,:),'.-k');
    endfor
    axis square; grid on; xlabel('SV length'); ylabel('Number of SVs');
    legend('baseline');
    title('sniffles2','fontsize', 18);
    set(gca, 'fontsize', 18);
endfunction
