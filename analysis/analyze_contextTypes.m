
MATRIX_DIR='/Users/fcunial/Downloads/performanceMatrices_contextTypes';
CALLERS={'pbsv', 'sniffles1', 'sniffles2'};
CONTEXT_TYPES=[0,0; 1,1; 4,4;];
READ_LENGTHS=[10000, 12500, 15000, 17500, 20000, 22500];
COVERAGES=[4, 8, 12, 16, 20];
MEASURES={'tp', 'fp', 'fn', 'precision', 'recall', 'f1'};
N_INDIVIDUALS=300;
DELTA=1000;
FONTSIZE=12;


% 1. Per-individual plots
COVERAGE_LINES={'.b','.c','.y','.m','.r'};
LEGEND={'4', '8', '12', '16', '20'};
[nRows,nColumns]=size(CONTEXT_TYPES);
for clr = [1:length(CALLERS)]
    for row = [1:nRows]
        figure((clr-1)*(1+nRows)+1+row);
        ct1=CONTEXT_TYPES(row,1); ct2=CONTEXT_TYPES(row,2);
        for ms = [1:length(MEASURES)]
            x=zeros(length(COVERAGES),N_INDIVIDUALS*length(READ_LENGTHS));
            y=zeros(length(COVERAGES),N_INDIVIDUALS*length(READ_LENGTHS));
            lastX=zeros(length(COVERAGES),1);
            fid=-1;
            try
                fid=fopen(sprintf('%s/%s_rs%d_re%d_matrix_%s.txt',MATRIX_DIR,CALLERS{clr},ct1,ct2,MEASURES{ms}));
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
            subplot(2,length(MEASURES)/2,ms); hold on;
            for coverage = [1:length(COVERAGES)]
                WOBBLE=(rand(1,lastX(coverage))-0.5)*DELTA;
                plot(x(coverage,[1:lastX(coverage)])+WOBBLE, y(coverage,1:lastX(coverage)), COVERAGE_LINES{coverage});
            endfor
            xlabel('avg read length'); axis square; grid on; set(gca, 'fontsize', FONTSIZE);
            title({MEASURES{ms},CALLERS{clr},sprintf('start=%d end=%d',ct1,ct2)}, 'fontsize', FONTSIZE);
        endfor
        for i = [1:3]
            subplot(2,length(MEASURES)/2,i); xlim([READ_LENGTHS(1)-DELTA, READ_LENGTHS(length(READ_LENGTHS))+DELTA]);
        endfor
        for i = [4:length(MEASURES)]
            subplot(2,length(MEASURES)/2,i); axis([READ_LENGTHS(1)-DELTA, READ_LENGTHS(length(READ_LENGTHS))+DELTA, 0, 1]);
        endfor
        subplot(2,length(MEASURES)/2,1); legend(LEGEND,'location','south','orient','horizontal');
    endfor
endfor
lastFigure=100;


% 2. Merged plots
COVERAGE_LINES={'-.b','-.c','-.y','-.m','-.r'};
LEGEND={'4 merge', '8 merge', '12 merge', '16 merge', '20 merge'};
for clr = [1:length(CALLERS)]
    for row = [1:nRows]
        figure(lastFigure+(clr-1)*(1+nRows)+1+row);
        ct1=CONTEXT_TYPES(row,1); ct2=CONTEXT_TYPES(row,2);
        for ms = [1:length(MEASURES)]
            x=zeros(length(COVERAGES),length(READ_LENGTHS));
            y=zeros(length(COVERAGES),length(READ_LENGTHS));
            lastX=zeros(length(COVERAGES),1);
            fid=-1;
            try
                fid=fopen(sprintf('%s/%s_rs%d_re%d_matrix_merge_%s.txt',MATRIX_DIR,CALLERS{clr},ct1,ct2,MEASURES{ms}));
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
            subplot(2,length(MEASURES)/2,ms); hold on;
            for coverage = [1:length(COVERAGES)]
                plot(x(coverage,[1:lastX(coverage)]), y(coverage,1:lastX(coverage)), COVERAGE_LINES{coverage});
            endfor
            xlabel('avg read length'); axis square; grid on; set(gca, 'fontsize', FONTSIZE);
            title({MEASURES{ms},CALLERS{clr},sprintf('start=%d end=%d MERGE AND JOINT',ct1,ct2)}, 'fontsize', FONTSIZE);
        endfor
        for i = [1:3]
            subplot(2,length(MEASURES)/2,i); xlim([READ_LENGTHS(1)-DELTA, READ_LENGTHS(length(READ_LENGTHS))+DELTA]);
        endfor
        for i = [4:length(MEASURES)]
            subplot(2,length(MEASURES)/2,i); axis([READ_LENGTHS(1)-DELTA, READ_LENGTHS(length(READ_LENGTHS))+DELTA, 0, 1]);
        endfor
        subplot(2,length(MEASURES)/2,1); legend(LEGEND,'location','south','orient','horizontal');
    endfor
endfor


% 3. Joint plots
COVERAGE_LINES={'-*b','-*c','-*y','-*m','-*r'};
LEGEND={'4 merge', '8 merge', '12 merge', '16 merge', '20 merge', '4 joint', '8 joint', '12 joint', '16 joint', '20 joint'};
for clr = [1:length(CALLERS)]
    for row = [1:nRows]
        figure(lastFigure+(clr-1)*(1+nRows)+1+row);
        ct1=CONTEXT_TYPES(row,1); ct2=CONTEXT_TYPES(row,2);
        for ms = [1:length(MEASURES)]
            x=zeros(length(COVERAGES),length(READ_LENGTHS));
            y=zeros(length(COVERAGES),length(READ_LENGTHS));
            lastX=zeros(length(COVERAGES),1);
            fid=-1;
            try
                fid=fopen(sprintf('%s/%s_rs%d_re%d_matrix_joint_%s.txt',MATRIX_DIR,CALLERS{clr},ct1,ct2,MEASURES{ms}));
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
            subplot(2,length(MEASURES)/2,ms); hold on;
            for coverage = [1:length(COVERAGES)]
                plot(x(coverage,[1:lastX(coverage)]), y(coverage,1:lastX(coverage)), COVERAGE_LINES{coverage});
            endfor
            xlabel('avg read length'); axis square; grid on; set(gca, 'fontsize', FONTSIZE);
            title({MEASURES{ms},CALLERS{clr},sprintf('start=%d end=%d MERGE AND JOINT',ct1,ct2)}, 'fontsize', FONTSIZE);
        endfor
        for i = [1:3]
            subplot(2,length(MEASURES)/2,i); xlim([READ_LENGTHS(1)-DELTA, READ_LENGTHS(length(READ_LENGTHS))+DELTA]);
        endfor
        for i = [4:length(MEASURES)]
            subplot(2,length(MEASURES)/2,i); axis([READ_LENGTHS(1)-DELTA, READ_LENGTHS(length(READ_LENGTHS))+DELTA, 0, 1]);
        endfor
        subplot(2,length(MEASURES)/2,1); legend(LEGEND,'location','south','orient','horizontal');
    endfor
endfor
