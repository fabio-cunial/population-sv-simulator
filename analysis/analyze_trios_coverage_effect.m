function analyze_trios_coverage_effect(MATRIX_DIR)
    CALLERS={'pbsv', 'sniffles2'};
	RIGHT_COVERAGES=[2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5];  % Of each haplotype
    MIN_COVERAGE=min(RIGHT_COVERAGES).*2; MAX_COVERAGE=max(RIGHT_COVERAGES).*2;
	LEGEND={'5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15'};
	SV_LENGTHS=[0, 100, 200, 400, 600, 800, 1000, 2000, 4000, 6000, 8000, 10000];
    MIN_SV_LENGTH=min(SV_LENGTHS); MAX_SV_LENGTH=max(SV_LENGTHS);
    FONTSIZE=32;
	COVERAGE_LINES={'-k','-r','-g','-b','-y','-m','-c', '--k','--r','--g','--b','--y','--m','--c', '-.k','-.r','-.g','-.b','-.y','-.m','-.c', ':k',':r',':g',':b',':y',':m',':c'};
    CALLERS_LINES={'.-b', '.-r', '.--b', '.--r'};
	
    # 1 - Analysis without genotypes
    MEASURES={'precision', 'recall', 'f1'};
    TITLES={'Precision', 'Recall', 'F1'};
	totalRecall=zeros(length(CALLERS),length(RIGHT_COVERAGES));
    for clr = [1:length(CALLERS)]
		for matrix = [1:3]
			for ms = [1:length(MEASURES)]
				B=zeros(length(RIGHT_COVERAGES),length(SV_LENGTHS));
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
						B(find(RIGHT_COVERAGES==coverage),find(SV_LENGTHS==svLength))=value;
                    endwhile
                    fclose(fid);
                endif
				
				% SV length analysis
				figure(clr*3+matrix);
                subplot(1,length(MEASURES),ms); hold on;
                for coverage = [1:length(RIGHT_COVERAGES)]
                    plot(SV_LENGTHS,B(coverage,:), COVERAGE_LINES{coverage});
                endfor
                xlabel('SV length'); axis([MIN_SV_LENGTH,MAX_SV_LENGTH,0,1]); axis square; grid on; set(gca, 'fontsize', FONTSIZE);
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
	figure(100); hold on;
	for clr = [1:length(CALLERS)]
		plot(RIGHT_COVERAGES.*2,totalRecall(clr,:),CALLERS_LINES{clr});
		axis([MIN_COVERAGE-1,MAX_COVERAGE+1,0,1]); axis square; grid on;
		xlabel('Coverage', 'fontsize', FONTSIZE); ylabel('Recall', 'fontsize', FONTSIZE);
		%title(sprintf('%s, all SVs',CALLERS{clr}));
		set(gca, 'fontsize', FONTSIZE);
		maximum=totalRecall(clr,length(RIGHT_COVERAGES));
		minimum=totalRecall(clr,1);
		printf('%s delta: %f',CALLERS{clr}, 100*(maximum-minimum)/minimum );
	endfor
    
    # 2 - Genotype concordance
	totalGtConcordance1=zeros(length(CALLERS),length(RIGHT_COVERAGES));
    totalGtConcordance2=zeros(length(CALLERS),length(RIGHT_COVERAGES));
    for clr = [1:length(CALLERS)]
		for matrix = [1:3]
            # Recall matrix
            B=zeros(length(RIGHT_COVERAGES),length(SV_LENGTHS));
            try
                fid=fopen(sprintf('%s/long_coverage_%s_matrix%d_%s.txt',MATRIX_DIR,CALLERS{clr},matrix,'recall'));
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
					B(find(RIGHT_COVERAGES==coverage),find(SV_LENGTHS==svLength))=value;
                endwhile
                fclose(fid);
            endif
            mRecall=B;
            # TP matrix
            B=zeros(length(RIGHT_COVERAGES),length(SV_LENGTHS));
            try
                fid=fopen(sprintf('%s/long_coverage_%s_matrix%d_%s.txt',MATRIX_DIR,CALLERS{clr},matrix,'tp'));
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
					B(find(RIGHT_COVERAGES==coverage),find(SV_LENGTHS==svLength))=value;
                endwhile
                fclose(fid);
            endif
            mTP=B;
            # GT matrix
            B=zeros(length(RIGHT_COVERAGES),length(SV_LENGTHS));
            try
                fid=fopen(sprintf('%s/long_coverage_%s_matrix%d_%s.txt',MATRIX_DIR,CALLERS{clr},matrix,'gt'));
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
					B(find(RIGHT_COVERAGES==coverage),find(SV_LENGTHS==svLength))=value;
                endwhile
                fclose(fid);
            endif
            mGT=B;
            % Number of true calls for each SV length
            nTrueCallsPerLength=mTP(length(RIGHT_COVERAGES),:) ./ mRecall(length(RIGHT_COVERAGES),:);
				
			% SV length analysis
			figure(200+clr*3+matrix);
            subplot(1,2,1); hold on;
            for coverage = [1:length(RIGHT_COVERAGES)]
                plot(SV_LENGTHS, mGT(coverage,:)./nTrueCallsPerLength, COVERAGE_LINES{coverage});
            endfor
            xlabel('SV length'); axis([MIN_SV_LENGTH,MAX_SV_LENGTH,0,1]); axis square; grid on; set(gca, 'fontsize', FONTSIZE);
			if (matrix == 1)
				string='>=';
			elseif (matrix == 2)
				string='=';
			else
				string='<=';
			endif
            title({'Recall with exact genotype',CALLERS{clr},string}, 'fontsize', FONTSIZE);
			legend(LEGEND,'location','eastoutside');
            subplot(1,2,2); hold on;
            for coverage = [1:length(RIGHT_COVERAGES)]
                plot(SV_LENGTHS, mGT(coverage,:)./mTP(coverage,:), COVERAGE_LINES{coverage});
            endfor
            xlabel('SV length'); axis([MIN_SV_LENGTH,MAX_SV_LENGTH,0,1]); axis square; grid on; set(gca, 'fontsize', FONTSIZE);
			if (matrix == 1)
				string='>=';
			elseif (matrix == 2)
				string='=';
			else
				string='<=';
			endif
            title({'Genotype concordance',CALLERS{clr},string}, 'fontsize', FONTSIZE);
			legend(LEGEND,'location','eastoutside');
			
            % Global genotype concordance analysis
            if (matrix == 1)
                totalGtConcordance1(clr,:)=mGT(:,1)./nTrueCallsPerLength(1);
                totalGtConcordance2(clr,:)=mGT(:,1)./mTP(:,1);
            endif
		endfor
    endfor
	
	% Global analysis
	figure(100); hold on;
	for clr = [1:length(CALLERS)]
		plot(RIGHT_COVERAGES.*2,totalGtConcordance1(clr,:),CALLERS_LINES{length(CALLERS)+clr});
		axis([MIN_COVERAGE-1,MAX_COVERAGE+1,0,1]); axis square; grid on;
        set(gca, 'xtick', MIN_COVERAGE-1:1:MAX_COVERAGE+1);
        set(gca, 'ytick', 0:0.05:1);
        set(gca, 'fontsize', FONTSIZE);
		%xlabel('Coverage'); ylabel('Recall');
		%title(sprintf('%s, recall with exact genotype, all SVs',CALLERS{clr}));
		%set(gca, 'fontsize', FONTSIZE);
		maximum=totalGtConcordance1(clr,length(RIGHT_COVERAGES));
		minimum=totalGtConcordance1(clr,1);
		printf('%s delta: %f',CALLERS{clr}, 100*(maximum-minimum)/minimum );
	endfor
    legend({'pbsv','sniffles2','pbsv, exact genotypes','sniffles2, exact genotypes',},'location','southeast','fontsize', FONTSIZE);
	figure(300+2);
	for clr = [1:length(CALLERS)]
		subplot(1,length(CALLERS),clr);
		plot(RIGHT_COVERAGES.*2,totalGtConcordance2(clr,:),'.-');
		axis([MIN_COVERAGE,MAX_COVERAGE,0,1]); axis square; grid on;
		xlabel('Coverage'); ylabel('GT concordance');
		title(sprintf('%s, genotype concordance, all SVs',CALLERS{clr}));
		set(gca, 'fontsize', FONTSIZE);
		maximum=totalGtConcordance2(clr,length(RIGHT_COVERAGES));
		minimum=totalGtConcordance2(clr,1);
		printf('%s delta: %f',CALLERS{clr}, 100*(maximum-minimum)/minimum );
	endfor
endfunction