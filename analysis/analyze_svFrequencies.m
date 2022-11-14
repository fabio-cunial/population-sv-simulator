
MATRIX_DIR='/Users/fcunial/Downloads/performanceMatrices_svFrequencies';
CALLERS={'pbsv', 'sniffles1', 'sniffles2'};
SVFREQUENCIES={'0.005', '0.01', '0.05', '0.1', '0.25', '0.5', '1.0'};
READ_LENGTHS=[10000, 12500, 15000, 17500, 20000, 22500];
COVERAGES=[4, 8, 12, 16, 20];
MEASURES={'fn', 'recall'};
N_INDIVIDUALS=300;
DELTA=1000;


% 1. Per-individual plots
COVERAGE_LINES={'.b','.c','.y','.m','.r'};
LEGEND={'4', '8', '12', '16', '20'};
for clr = [1:length(CALLERS)]
    for svf = [1:length(SVFREQUENCIES)]
        figure((clr-1)*(1+length(SVFREQUENCIES))+1+svf);
        for ms = [1:length(MEASURES)]
            x=zeros(length(COVERAGES),N_INDIVIDUALS*length(READ_LENGTHS));
            y=zeros(length(COVERAGES),N_INDIVIDUALS*length(READ_LENGTHS));
            lastX=zeros(length(COVERAGES),1);
            fid=-1;
            try
                fid=fopen(sprintf('%s/%s_svf%s_matrix_%s.txt',MATRIX_DIR,CALLERS{clr},SVFREQUENCIES{svf},MEASURES{ms}));
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
                    readLength=A(1); coverage=A(2);
                    for i=[3:length(A)]
                        j=find(COVERAGES==coverage);
                        lastX(j)=lastX(j)+1;
                        x(j,lastX(j))=readLength;
                        y(j,lastX(j))=A(i);
                    endfor
                endwhile
                fclose(fid);
            endif
            subplot(1,length(MEASURES),ms); hold on;
            for coverage = [1:length(COVERAGES)]
                WOBBLE=(rand(1,lastX(coverage))-0.5)*DELTA;
                plot(x(coverage,[1:lastX(coverage)])+WOBBLE, y(coverage,1:lastX(coverage)), COVERAGE_LINES{coverage});
            endfor
            xlabel('avg read length'); axis square; grid on; %legend(LEGEND);
            title(sprintf('%s <=%s %s',CALLERS{clr},SVFREQUENCIES{svf},MEASURES{ms}));
        endfor
    endfor
endfor
lastFigure=100;


% 2. Merged plots
COVERAGE_LINES={'-.b','-.c','-.y','-.m','-.r'};
LEGEND={'4 merge', '8 merge', '12 merge', '16 merge', '20 merge'};
for clr = [1:length(CALLERS)]
    for svf = [1:length(SVFREQUENCIES)]
        figure(lastFigure+(clr-1)*(1+length(SVFREQUENCIES))+1+svf);
        for ms = [1:length(MEASURES)]
            x=zeros(length(COVERAGES),length(READ_LENGTHS));
            y=zeros(length(COVERAGES),length(READ_LENGTHS));
            lastX=zeros(length(COVERAGES),1);
            fid=-1;
            try
                fid=fopen(sprintf('%s/%s_svf%s_matrix_merge_%s.txt',MATRIX_DIR,CALLERS{clr},SVFREQUENCIES{svf},MEASURES{ms}));
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
                    readLength=A(1); coverage=A(2);
                    j=find(COVERAGES==coverage);
                    lastX(j)=lastX(j)+1;
                    x(j,lastX(j))=readLength;
                    y(j,lastX(j))=A(3);
                endwhile
                fclose(fid);
            endif
            subplot(1,length(MEASURES),ms); hold on;
            for coverage = [1:length(COVERAGES)]
                WOBBLE=(rand(1,lastX(coverage))-0.5)*DELTA;
                plot(x(coverage,[1:lastX(coverage)])+WOBBLE, y(coverage,1:lastX(coverage)), COVERAGE_LINES{coverage});
            endfor
            xlabel('avg read length'); axis square; grid on; %legend(LEGEND);
            title(sprintf('%s MERGE AND JOINT <=%s %s',CALLERS{clr},SVFREQUENCIES{svf},MEASURES{ms}));
        endfor
    endfor
endfor


% 3. Joint plots
COVERAGE_LINES={'-*b','-*c','-*y','-*m','-*r'};
LEGEND={'4 merge', '8 merge', '12 merge', '16 merge', '20 merge', '4 joint', '8 joint', '12 joint', '16 joint', '20 joint'};
for clr = [1:length(CALLERS)]
    for svf = [1:length(SVFREQUENCIES)]
        figure(lastFigure+(clr-1)*(1+length(SVFREQUENCIES))+1+svf);
        for ms = [1:length(MEASURES)]
            x=zeros(length(COVERAGES),length(READ_LENGTHS));
            y=zeros(length(COVERAGES),length(READ_LENGTHS));
            lastX=zeros(length(COVERAGES),1);
            fid=-1;
            try
                fid=fopen(sprintf('%s/%s_svf%s_matrix_joint_%s.txt',MATRIX_DIR,CALLERS{clr},SVFREQUENCIES{svf},MEASURES{ms}));
            catch
                % NOP
            end_try_catch
            fid
            sprintf('%s/%s_svf%s_matrix_joint_%s.txt',MATRIX_DIR,CALLERS{clr},SVFREQUENCIES{svf},MEASURES{ms})
            if (fid>=0)
                while (true)
                    str=fgetl(fid);
                    if (str==-1)
                        break
                    endif
                    A=str2num(str);
                    readLength=A(1); coverage=A(2);
                    j=find(COVERAGES==coverage);
                    lastX(j)=lastX(j)+1;
                    x(j,lastX(j))=readLength;
                    y(j,lastX(j))=A(3);
                endwhile
                fclose(fid);
            endif
            subplot(1,length(MEASURES),ms); hold on;
            for coverage = [1:length(COVERAGES)]
                WOBBLE=(rand(1,lastX(coverage))-0.5)*DELTA;
                plot(x(coverage,[1:lastX(coverage)])+WOBBLE, y(coverage,1:lastX(coverage)), COVERAGE_LINES{coverage});
            endfor
            xlabel('avg read length'); axis square; grid on; %legend(LEGEND);
            title(sprintf('%s MERGE AND JOINT <=%s %s',CALLERS{clr},SVFREQUENCIES{svf},MEASURES{ms}));
        endfor
    endfor
endfor
