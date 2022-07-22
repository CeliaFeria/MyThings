
require "./CTTCTT_information.rb"

get_gff3_gene(get_exon(get_embl(get_AGI_Locus('./ArabidopsisSubNetwork_GeneList.txt'))))

get_gff3_chrom(get_exon(get_embl(get_AGI_Locus('./ArabidopsisSubNetwork_GeneList.txt'))))

report_no_CTT(get_exon(get_embl(get_AGI_Locus('./ArabidopsisSubNetwork_GeneList.txt'))))


