function analyze_trios(MATRIX_DIR)
    CALLERS={'pbsv', 'sniffles1', 'sniffles2'};
    MEASURES={'precision', 'recall', 'f1'};
	TITLES={'Precision', 'Recall', 'F1'};
	QUESTION2_LEFT_COVERAGES=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 3.0, 4.0];
	LEGEND={'0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0', '1.1', '1.2', '1.3', '1.4', '1.5', '2.0', '3.0', '4.0'};
	SV_LENGTHS=[0, 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000];
    FONTSIZE=12;
	COVERAGE_LINES={'-k','-r','-g','-b','-y','-m','-c', '--k','--r','--g','--b','--y','--m','--c', '-.k','-.r','-.g','-.b','-.y','-.m','-.c', ':k',':r',':g',':b',':y',':m',':c'};
	
	totalRecall=zeros(length(CALLERS),length(QUESTION2_LEFT_COVERAGES));
    for clr = [1:length(CALLERS)]
		for matrix = [1:3]
			for ms = [1:length(MEASURES)]
				B=zeros(length(QUESTION2_LEFT_COVERAGES),length(SV_LENGTHS));
	            try
	                fid=fopen(sprintf('%s/long_coverage_%s_matrix%d_%s.txt',MATRIX_DIR,CALLERS{clr},matrix,MEASURES{ms}));
	            catch
	                % NOP
	            end_try_catch
                if (fid>=0)
                    while (true)
                        str=fgetl(fid);
                        if (str==-1)
                            break
                        endif
                        A=str2num(str);
                        coverage=A(1); svLength=A(2); value=A(3);
						B(find(QUESTION2_LEFT_COVERAGES==coverage),find(SV_LENGTHS==svLength))=value;
                    endwhile
                    fclose(fid);
                endif
				
				% SV length analysis
				figure(clr*3+matrix);
                subplot(1,length(MEASURES),ms); hold on;
                for coverage = [1:length(QUESTION2_LEFT_COVERAGES)]
                    plot(SV_LENGTHS,B(coverage,:), COVERAGE_LINES{coverage});
                endfor
                xlabel('SV length'); axis square; grid on; set(gca, 'fontsize', FONTSIZE);
				if (matrix == 1)
					string='>=';
				elseif (matrix == 2)
					string='=';
				else
					string='<=';
				endif
                title({TITLES{ms},CALLERS{clr},string}, 'fontsize', FONTSIZE);
				legend(LEGEND,'location','eastoutside');
				
				% Global analysis
				if (ms == 2 && matrix == 1)  % Only recall for now
					totalRecall(clr,:)=B(:,1)';
				endif
			endfor
		endfor
    endfor
	
	% Global analysis
	figure(100);
	for clr = [1:length(CALLERS)]
		subplot(1,length(CALLERS),clr);
		plot(QUESTION2_LEFT_COVERAGES,totalRecall(clr,:),'.-');
		axis([0,1,min(totalRecall(clr,:)),max(totalRecall(clr,:))]);
		axis square; grid on;
		xlabel('Short-read coverage of each haplotype'); ylabel('Recall');
		title(sprintf('%s, all SVs',CALLERS{clr}));
		set(gca, 'fontsize', FONTSIZE);
		maximum=totalRecall(clr,10);
		minimum=totalRecall(clr,1);
		printf('%s delta: %f',CALLERS{clr}, 100*(maximum-minimum)/minimum );
	endfor
endfunction