
% Darlan Conterno Minussi
% Inclass assignment 9

% The accession number for human NOTCH1 mRNA is AF308602
% 1. Read the information from this entry into matlab

notch = getgenbank('AF308602');

% 2. Write code that runs a blast query on the first 500 base pairs of this
% gene against the refseq_rna database

[RID, RT0E] = blastncbi(notch.Sequence(1:500), 'blastn', 'Database', 'refseq_rna');

notchblast = getblast(RID, 'WaitTime', RT0E);

% 3. Find the three highest scoring hits from other species and identify
% the length of the alignment and fraction of matches/mismatches.

notchblast_table = struct2table(notchblast.Hits);

for i = 1:size(notchblast_table,1)
    notchblast_table.score(i) = notchblast.Hits(i).HSPs.Score;
end

% removing Homo Sapiens
notchtable_wouth = notchblast_table(cellfun(@isempty, regexp(notchblast_table.Name, 'Homo sapiens')), :);

% Exporting 6 values because the first 3 one are for the same species just
% different transcript variants, therefore, making sure 3 different species
% are being shown in the output.
% 
findAlignment(notchtable_wouth, notchblast, 6);


% 4. Run the same query against the database est_human. Comment on the
% sequences that you find. 

[RID, RT0E] = blastncbi(notch.Sequence(1:500), 'blastn', 'Database', 'est_human');
notchblast_est = getblast(RID, 'WaitTime', RT0E);

% They all show small human sequences that are from cDNAs that came from
% assays where the expression of mRNA was being measured in different
% tissues. The alignment is an indirect indication of expression of that
% particular mRNA in the tissue

% Definition from: https://www.ncbi.nlm.nih.gov/genbank/dbest/
% Expressed Sequence Tags (ESTs) are short (usually <1000 bp), single-pass sequence reads from mRNA (cDNA). 
% Typically they are produced in large batches. 
% They represent a snapshot of genes expressed in a given tissue and/or at a given developmental stage. 
% They are tags (some coding, others not) of expression for a given cDNA library.