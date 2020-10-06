%%% NCBI nt(nucleotide) position locator %%%

% Developer: Yanshi Hu, M.S., Ph.D. Candidate
% Affiliation: Department of Bioinformatics, State Key Laboratory of Plant Physiology and Biochemistry, College of Life Sciences, Zhejiang University


% ----------------USAGE-------------------%
% INPUT FILE FORMAT:  xls, xlsx
% FORMAT INSIDE THE FILE SHOULD BE: gene symbol, mutated aa sequence, score
% (NOTE: score is calculated as follows: 
%                 positive number refers to the rank of the mutated aa from left to right, 
%                 with non-positive number from right to left. Zero is given when the most right one is the mutated aa and -1 is assigned if the mutation is the second one from right to left.
%  )

% essential variables: 
%             q: original info table
%    codontable: condon table for 20 aa in human (aa: amino acid)


% -----------------------START-------------------------%

% Selecting & Reading file into Matlab
[filename,path] = uigetfile({'*.xls;*.xlsx','Query Files (*.xls,*.xlsx)'},'Select Query File');
if (filename==0), return, end
filename = fullfile(path,filename);
[~,~,q] = xlsread(filename);
load('codon_table.mat')

pause(3);
disp('Welcome to use NCBI nucleotide locator!');
pause(5);
disp('--------------Program starts--------------');
pause(1);
% Finding exact NCBI Gene page for a specific gene symbol
a = {};
b = {}; % gene symbol in the requested table
c = {}; % source code for Entrez Gene page of a specific gene
ntn = ''; % nt after mutation
ntp = ''; % nt position
ntseq = ''; % nt sequence
ntorigin = ''; % original nt before mutation
finalist = {}; % final info table --------- %% FORMAT: gene symbol, mutant sequence, aa substitution info, mRNA variable(s), chromosome mutation info %%---------
lsnum = 1; % final info table index
chroinfo = {}; % chromosome info

for i = 1:length(q(:, 1))
    i
    alist = {}; % 'T190I' pattern storage
    aalist = {}; % 'S>F' pattern storage

    % Processing for the requested table %   
    a = regexp(q{i, 1}, '((?<=[£¨\(]+).*(?=[\)£©]+))', 'match');
    if isempty(a) % Currently, incompatible with this kind of gene format!
        finalist(lsnum, 1:5) = [q(i, 1), {upper(q{i, 2})}, {''}, {''}, {''}]; 
        lsnum = lsnum + 1;
        continue
    end
    if isempty(intersect(a{1}, '>')) % 'T190I'pattern
        alist = regexp(a{1}, '\d+|\w', 'match'); % storing 'T', '190' & 'I' separately
    else % 'S>F'pattern
        aalist = regexp(a{1}, '\w', 'match'); % storing 'S' & 'F'
    end

    % Extracting mRNA NCBI Reference Sequence IDs for a specific gene in the requested table %
    b = regexp(q{i, 1}, '^[^\s(£¨]+', 'match'); % extracting gene symbol
    
    if strcmp(b{1}(1), '¦Â')
        qqq = b{1}(2:length(b{1}));
        b{1} = ['%CE%B2', qqq];
    end        
    c = urlread(['https://www.ncbi.nlm.nih.gov/gene/?term=(', b{1}, '%5Bgene%5D)+AND+(Homo+sapiens%5Borgn%5D)']); % Finding NCBI Gene page
    d = regexp(c, '<h3 class="genomerefseqs gene_std_toggler" id="refseqset-independent">.*<h3 class="genomerefseqs gene_std_toggler" id="annotatedrefseq-title">', 'match');
    if isempty(d)
        e = regexp(c, '(?<=>ID:\s+)\d*', 'match');
        if ~isempty(e)
            c = urlread(['https://www.ncbi.nlm.nih.gov/gene/', e{1}]);
            d = regexp(c, '<h3 class="genomerefseqs gene_std_toggler" id="refseqset-independent">.*<h3 class="genomerefseqs gene_std_toggler" id="annotatedrefseq-title">', 'match');
        else % When the above strategy isn't suitable.
            e = urlread(['https://www.ncbi.nlm.nih.gov/gene?term=(', b{1},')+AND+(Homo+sapiens%5Bporgn:__txid9606%5D)']); % Alternative search strategy
            d = regexp(e, '<h3 class="genomerefseqs gene_std_toggler" id="refseqset-independent">.*<h3 class="genomerefseqs gene_std_toggler" id="annotatedrefseq-title">', 'match');
            if isempty(d)
                d = regexp(e, '(?<=>ID:\s+)\d*', 'match');
                if ~isempty(d)
                    c = urlread(['https://www.ncbi.nlm.nih.gov/gene/', d{1}]);
                    d = regexp(c, '<h3 class="genomerefseqs gene_std_toggler" id="refseqset-independent">.*<h3 class="genomerefseqs gene_std_toggler" id="annotatedrefseq-title">', 'match');
                else
                    finalist(lsnum, 1:5) = [q(i, 1), {upper(q{i, 2})}, {''}, {''}, {''}]; % No gene of Homo sapiens found in NCBI Entrez Gene database!
                    lsnum = lsnum + 1;
                    continue
                end
            end
        end
    end
        
    e = regexp(d{1}, '(?<=>)NM_[\d\.]*', 'match')'; % e: cell matrix for all mRNA NCBI Reference Sequence IDs of a specific gene

    % 'T190I' pattern OR 'S>F' pattern selection    
    if isempty(aalist) % 'T190I' pattern processing
        for j = 1:length(e)
            f = urlread(['https://www.ncbi.nlm.nih.gov/nuccore/', e{j, 1}]);
            g = regexp(f, '(?<=<meta name="ncbi_uidlist" content=")\d*', 'match'); % gi No.
            gbinfo = urlread(['https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=genbank&id=', g{1}, '&maxplex=1']); % genebank info
            cdss = regexp(gbinfo, '(?<=CDS\s+)\d*(?=..)', 'match'); % start nt position for CDS
            chroinfo = regexp(gbinfo, '(?<=chromosome=")\S+(?=")', 'match'); % chromosome info
            % Due to no need for aa sequence mapping
    %         f = regexp(gbinfo, '(?<=/translation=")[^"]*(?=")', 'match'); 
    %         aaseq = regexprep(f{1}, '\s+', ''); % aa sequence
            % Due to no need for aa sequence mapping
            ntfile = urlread(['https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=',g{1}]);
            ntfile1 = regexprep(ntfile, '>NM_.*mRNA\n', '');
            ntseq = regexprep(ntfile1, '\n', ''); % nt sequence

            % Starting requested sequences comparison
            ntstart = str2double(cdss{1}) + (str2double(alist{1, 2}) -1) * 3; % Locating the starting nt for the mutated aa 
            if ntstart > length(ntseq) || ntstart + 2 > length(ntseq)
                continue
            end
            norcodon = ntseq(ntstart:ntstart + 2); % normal codon for this aa
            ntn = {}; % cell array for mutant codon storage
            ntorigin = {}; % cell array for original codon storage
            ntp = []; % matrix array for mutated nt storage
            nnt = 0; % initial index for the codon usage selection
            for ii = 1 : 2 : 40 % finding the mutated aa & its codon in codontable
                if isequal(codontable(ii, 1), alist(1, 3))
                    cc = codontable(ii + 1, :);
                    cc(cellfun(@isempty, cc)) = []; % using cellfun function to remove all empty elements in cells
                    for iii = 1:length(cc)
                        if norcodon(1) ~= cc{1, iii}(1) && norcodon(2) == cc{1, iii}(2) && norcodon(3) == cc{1, iii}(3)
                            nnt = nnt + 1;
                            ntn{nnt} = cc{1, iii}(1);
                            ntorigin{nnt} = norcodon(1);
                            ntp(nnt) = ntstart;
                        else if norcodon(1) == cc{1, iii}(1) && norcodon(2) ~= cc{1, iii}(2) && norcodon(3) == cc{1, iii}(3)
                                nnt = nnt + 1;
                                ntn{nnt} = cc{1, iii}(2);
                                ntorigin{nnt} = norcodon(2);
                                ntp(nnt) = ntstart + 1;
                            else if norcodon(1) == cc{1, iii}(1) && norcodon(2) == cc{1, iii}(2) && norcodon(3) ~= cc{1, iii}(3)
                                    nnt = nnt + 1;
                                    ntn{nnt} = cc{1, iii}(3);
                                    ntorigin{nnt} = norcodon(3);
                                    ntp(nnt) = ntstart + 2; 
                                end
                            end
                        end
                    end
                end
            end 
            if j == 1
                finalist(lsnum, 1:2) = [q(i, 1), {upper(q{i, 2})}];
                finalist(lsnum, 3) = a;
            else
                finalist(lsnum, [1 2 3]) = {'', '', ''};
            end
            finalist(lsnum, 4) = e(j);
            for nntt = 1 : nnt
                finalist(lsnum, 4 + nntt) = {['chr', chroinfo{1}, '_', num2str(ntp(nntt)), '_', ntorigin{nntt}, '-', ntn{nntt}]};
            end
            lsnum = lsnum + 1;      
        end
        
    else % 'S>F' pattern processing
        for j = 1:length(e)
            f = urlread(['https://www.ncbi.nlm.nih.gov/nuccore/', e{j, 1}]);
            g = regexp(f, '(?<=<meta name="ncbi_uidlist" content=")\d*', 'match'); % gi No.
            gbinfo = urlread(['https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=genbank&id=', g{1}, '&maxplex=1']); % genebank info
            cdss = regexp(gbinfo, '(?<=CDS\s+)\d*(?=..)', 'match'); % start nt position for CDS
            chroinfo = regexp(gbinfo, '(?<=chromosome=")\S+(?=")', 'match'); % chromosome info
            f = regexp(gbinfo, '(?<=/translation=")[^"]*(?=")', 'match'); 
            aaseq = regexprep(f{1}, '\s+', ''); % aa sequence
            ntfile = urlread(['https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=',g{1}]);
            ntfile1 = regexprep(ntfile, '>NM_.*mRNA\n', '');
            ntseq = regexprep(ntfile1, '\n', ''); % nt sequence

            % Starting requested sequences comparison
               % mutated aa rank determination
            qq = q{i, 2};
            if q{i, 3} > 0
                qq(q{i, 3}) = aalist{1}; % original aa query sequence before mutation
                initaa = strfind(aaseq, upper(qq)); % OR initaa = findstr(qq, aaseq); OR initaa = findstr(aaseq, qq);  %% initial aa location for query aa sequence in the whole aa sequence
                if isempty(initaa)
                    if j ~= length(e)
                        continue
                    else
                        finalist(lsnum, 1:5) = [q(i, 1), {upper(q{i, 2})}, {''}, {''}, {''}]; % No matching sequence!
                        lsnum = lsnum + 1;
                        break
                    end
                end
                aalen = initaa + q{i, 3} - 1; % mutated aa location
                ntstart = str2double(cdss{1}) + (initaa + q{i, 3} - 2) * 3; % Locating the starting nt for the mutated aa 
            else
                qq(length(q{i, 2}) + q{i, 3}) = aalist{1}; % original aa query sequence before mutation
                initaa = strfind(aaseq, upper(qq)); % initial aa location for query aa sequence in the whole aa sequence
                if isempty(initaa)
                    if j ~= length(e)
                        continue
                    else
                        finalist(lsnum, 1:5) = [q(i, 1), {upper(q{i, 2})}, {''}, {''}, {''}]; % No matching sequence!
                        lsnum = lsnum + 1;
                        break
                    end
                end
                aalen = initaa + length(q{i, 2}) + q{i, 3} - 1; % mutated aa location
                ntstart = str2double(cdss{1}) + (initaa + length(q{i, 2}) + q{i, 3} - 2) * 3; % Locating the starting nt for the mutated aa                
            end
               % locating the original query sequence in the whole aa sequence
            
            norcodon = ntseq(ntstart:ntstart + 2); % normal codon for this aa
            ntn = {}; % cell array for mutant codon storage
            ntorigin = {}; % cell array for original codon storage
            ntp = []; % matrix array for mutated nt storage
            nnt = 0; % initial index for the codon usage selection
            for ii = 1 : 2 : 40 % finding the mutated aa & its codon in codontable
                if isequal(codontable(ii, 1), aalist(1, 2))
                    cc = codontable(ii + 1, :);
                    cc(cellfun(@isempty, cc)) = []; % using cellfun function to remove all empty elements in cells
                    for iii = 1:length(cc)
                        if norcodon(1) ~= cc{1, iii}(1) && norcodon(2) == cc{1, iii}(2) && norcodon(3) == cc{1, iii}(3)
                            nnt = nnt + 1;
                            ntn{nnt} = cc{1, iii}(1);
                            ntorigin{nnt} = norcodon(1);
                            ntp(nnt) = ntstart;
                        else if norcodon(1) == cc{1, iii}(1) && norcodon(2) ~= cc{1, iii}(2) && norcodon(3) == cc{1, iii}(3)
                                nnt = nnt + 1;
                                ntn{nnt} = cc{1, iii}(2);
                                ntorigin{nnt} = norcodon(2);
                                ntp(nnt) = ntstart + 1;
                            else if norcodon(1) == cc{1, iii}(1) && norcodon(2) == cc{1, iii}(2) && norcodon(3) ~= cc{1, iii}(3)
                                    nnt = nnt + 1;
                                    ntn{nnt} = cc{1, iii}(3);
                                    ntorigin{nnt} = norcodon(3);
                                    ntp(nnt) = ntstart + 2; 
                                end
                            end
                        end
                    end
                end        
            end
            if j == 1
                finalist(lsnum, 1:2) = [q(i, 1), {upper(q{i, 2})}];
            else
                finalist(lsnum, [1 2]) = {'', ''};
            end
            finalist(lsnum, 3) = {[aalist{1}, num2str(aalen), aalist{2}]};
            finalist(lsnum, 4) = e(j);
            for nntt = 1 : nnt
                finalist(lsnum, 4 + nntt) = {['chr', chroinfo{1}, '_', num2str(ntp(nntt)), '_', ntorigin{nntt}, '-', ntn{nntt}]};
            end
            lsnum = lsnum + 1;      
        end
    end
end

% Save file & final data is written into excel
[filename, pathname] = uiputfile(...
 {'*.xls';'*.xlsx';'*.*'},...
 'Save as');
xlswrite([pathname, filename], finalist);

disp('--------------Program ends--------------');
    

% --------------------------------------END----------------------------------------------%
