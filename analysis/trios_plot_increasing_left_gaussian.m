function plot_increasing_left_gaussian()
    # --------------------------------------------------------------------------
    # Sampling several length coverages
    figure(1);

    # Top figure
    subplot(2,1,1);
    bins=load('scratch_cunial-read-length-distribution_2842314_bins_bin_.histogram');
    plot(bins(:,1),bins(:,2))
    axis square; grid on; xlabel('Read length'); ylabel('Number of reads'); 
    %title('Observed');
    title('All reads from child 2842314');
    set(gca, 'fontsize', 18);

    # Bottom figure
    subplot(2,1,2)
    LEGEND={'0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0'};
    id='0.0'; A0=load(sprintf('scratch_cunial-read-length-distribution_2842314_reads_c%s_coverage_%s.fastq.histogram',id,id));
    id='0.1'; A1=load(sprintf('scratch_cunial-read-length-distribution_2842314_reads_c%s_coverage_%s.fastq.histogram',id,id));
    id='0.2'; A2=load(sprintf('scratch_cunial-read-length-distribution_2842314_reads_c%s_coverage_%s.fastq.histogram',id,id));
    id='0.3'; A3=load(sprintf('scratch_cunial-read-length-distribution_2842314_reads_c%s_coverage_%s.fastq.histogram',id,id));
    id='0.4'; A4=load(sprintf('scratch_cunial-read-length-distribution_2842314_reads_c%s_coverage_%s.fastq.histogram',id,id));
    id='0.5'; A5=load(sprintf('scratch_cunial-read-length-distribution_2842314_reads_c%s_coverage_%s.fastq.histogram',id,id));
    id='0.6'; A6=load(sprintf('scratch_cunial-read-length-distribution_2842314_reads_c%s_coverage_%s.fastq.histogram',id,id));
    id='0.7'; A7=load(sprintf('scratch_cunial-read-length-distribution_2842314_reads_c%s_coverage_%s.fastq.histogram',id,id));
    id='0.8'; A8=load(sprintf('scratch_cunial-read-length-distribution_2842314_reads_c%s_coverage_%s.fastq.histogram',id,id));
    id='0.9'; A9=load(sprintf('scratch_cunial-read-length-distribution_2842314_reads_c%s_coverage_%s.fastq.histogram',id,id));
    id='1.0'; A10=load(sprintf('scratch_cunial-read-length-distribution_2842314_reads_c%s_coverage_%s.fastq.histogram',id,id));
    plot(A0(:,1), [ A0(:,2), A1(:,2), A2(:,2), A3(:,2), A4(:,2), A5(:,2), A6(:,2), A7(:,2), A8(:,2), A9(:,2), A10(:,2) ]);
    axis square; grid on; xlabel('Read length'); ylabel('Number of reads'); 
    title('Sampled reads from child 2842314');
    legend(LEGEND,'location','eastoutside');
    set(gca, 'fontsize', 18);


    # --------------------------------------------------------------------------
    # Creating ground truth
    figure(2);

    subplot(3,2,1);
    plot(bins(:,1),bins(:,2));
    axis square; grid on; xlabel('Read length'); ylabel('Number of reads'); 
    title('All reads from child 2842314');
    set(gca, 'fontsize', 14);

    subplot(3,2,2);
    maxCoverageChild=load('scratch_cunial-read-length-distribution_2842314_long_coverage_2842314_reads.fastq.histogram');
    plot(maxCoverageChild(:,1),maxCoverageChild(:,2));
    axis square; grid on; xlabel('Read length'); ylabel('Number of reads'); 
    title('Long reads from child 2842314');
    set(gca, 'fontsize', 14);

    subplot(3,2,4);
    maxCoverageParent1=load('scratch_cunial-read-length-distribution_2842314_long_coverage_1146737_reads.fastq.histogram');
    plot(maxCoverageParent1(:,1),maxCoverageParent1(:,2));
    axis square; grid on; xlabel('Read length'); ylabel('Number of reads'); 
    title('Long reads from parent 1146737');
    set(gca, 'fontsize', 14);

    subplot(3,2,6);
    maxCoverageParent2=load('scratch_cunial-read-length-distribution_2842314_long_coverage_1701906_reads.fastq.histogram');
    plot(maxCoverageParent2(:,1),maxCoverageParent2(:,2));
    axis square; grid on; xlabel('Read length'); ylabel('Number of reads'); 
    title('Long reads from parent 1701906');
    set(gca, 'fontsize', 14);
endfunction
