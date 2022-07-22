
#This code search for CTT repeats, creates a gff3 (chromosome or genes) and create a report file with the genes without CTT repeats


require 'bio'
require 'rest-client' 


#Obtain the information of an URL

def fetch(url, headers = {accept: "*/*"}, user = "", pass="")
  response = RestClient::Request.execute({
    method: :get,
    url: url.to_s,
    user: user,
    password: pass,
    headers: headers})
  return response
  
  rescue RestClient::ExceptionWithResponse => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  rescue RestClient::Exception => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  rescue Exception => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
end 

#Create an array with gene names

  # @param path[String] the path of the file where are the names of the genes
  # @return [locus_code] array with all the genes in the file
def get_AGI_Locus(path)
  locus_code= []
  locus = File.open(path, mode: 'r')
  locus.readlines.each do |line| 
    code = line.strip.split("\n")
  locus_code |= code
  end 
  return locus_code
end  



#Create a hash with gene names and embl information

  # @param all_locus[String] array with genes
  # @return [embl] hash with the name of the gene as key and Bio::EMBL object as value
def get_embl(all_locus)
  embl = Hash.new
  all_locus.each do |locus|
      response = fetch("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{locus}")
      next if response.nil?
      embl[locus]= Bio::EMBL.new(response.body)
  end
  return embl
end



#Add new features with positions on genes of CTT repeats  

  # @param embl[Hash] hash with the name of the gene as key and Bio::EMBL object as value
def get_exon(embl)
  
  embl.each do |key, info|
    info = info.to_biosequence #we need to convert the database entry to a Bio::Sequence in order to add the new features
     all_starts_f = [] #array to keep the forward locations without repetitions
     all_starts_r = [] #array to keep the reverse locations without repetitions
    
    info.features.each do |feature|
      next if feature.feature != "exon"
      feature.locations.each do |location|
      exon_seq=info.subseq(location.from,location.to)
      exon_start = location.from #keep the start position of the exon in order to have then the location of the repeats in the gene
        
      next if exon_seq.nil?
        
        if location.strand == 1
          #When the strand is positive we only have to search for the cttctt repeats
          #In order to get the the overlapping matches I used lookahead assertion (?=)
          #In order to get all the positions I used the code from this website, that generates an array with all start positions
          #https://stackoverflow.com/questions/5241653/ruby-regex-match-and-get-positions-of/5241843#5241843
          start_f = exon_seq.enum_for(:scan, /(?=(cttctt))/).map { Regexp.last_match.begin(0)} 
          next if start_f[0].nil?
          
          start_f.each do |x|
          #@Start_f has the positions starting in 0, but the feature position start in 1, we could sum +1 in order to get th same position.
          #  Instead of sum 1 at the start_f positions I leave it like that and then when exon_start is added (start in 1) the problem is solved
          next if all_starts_f.include?(x + exon_start) #skip repeated locations
          add_feature = Bio::Feature.new(feature = 'CTTCTT', position = "#{x+exon_start}..#{x + exon_start+5}")
          add_feature.append(Bio::Feature::Qualifier.new('repeat_motif', 'CTTCTT'))
          add_feature.append(Bio::Feature::Qualifier.new('strand', '+'))
          info.features << add_feature 
            all_starts_f.push(x + exon_start)
          end
          
        end
        
        if location.strand == -1
          #@When the strand is negative we can search in the complement strand the cttctt repeats or we can search 
          #  for the complement of cttctt  ( aagaag) in the positive strand 
          #I used the second option, because I find it easier
          #In order to get the the overlapping matches I used lookahead assertion (?=)
          #In order to get all the positions I used the code from this website, that generates an array with all start positions
          #https://stackoverflow.com/questions/5241653/ruby-regex-match-and-get-positions-of/5241843#5241843
          start_r = exon_seq.enum_for(:scan, /(?=(aagaag))/).map { Regexp.last_match.begin(0)} #is not necessary to add 1 because when we sum the location of the exon is like sum 1 already
          next if start_r[0].nil?
          
          start_r.each do |x|
          #@start_r has the positions starting in 0, but the feature position start in 1, we could sum +1 in order to get th same position.
          #  Instead of sum 1 at the start_r positions I leave it like that and then when exon_start is added (start in 1) the problem is solved
          next if all_starts_r.include?(x + exon_start) #skip repeated locations
          add_feature = Bio::Feature.new(feature = 'CTTCTT',location = "complement(#{x + exon_start}..#{x + exon_start +5})") #complement in order to indicate that they are in the negative strand
          add_feature.append(Bio::Feature::Qualifier.new('repeat_motif', 'CTTCTT'))
          add_feature.append(Bio::Feature::Qualifier.new('strand', '-'))
          info.features << add_feature 
            all_starts_r.push(x + exon_start)
          end
        end
     end
    end
  end

end







#Create a gff3 file with the positions of CTTCTT repeats on the genes

  # @param new_embl[Hash] hash with the name of the gene as key and Bio::Sequence object with new features added
  # @return [gff3_exon] a GFF3 file 
def get_gff3_gene(new_embl)
  
#seqid - ID of the gene
#type - type of feature. Must be a term or accession from the SOFA sequence ontology
#start - Start position of the feature, with sequence numbering starting at 1.
#end - End position of the feature, with sequence numbering starting at 1.
#score - A floating point value.
#strand - defined as + (forward) or - (reverse).
#phase - . 
#attributes - A semicolon-separated list of tag-value pairs, providing additional information.
  gff3_exon = File.open('gff3_exon', mode: 'w')
  gff3_exon.write("##gff-version 3\n") #the first line of the gff3 file
  type = "SO:0000657" #Sequence Ontology repeat-region
  source = "Ensembl"
  score = "."
  phase = "."
  
  new_embl.each do |key, info|
    info = info.to_biosequence
    count = 0
    info.features.each do |feature|
      next if feature.feature != "CTTCTT"
      count += 1
      seqid = "#{key}" #this is the seqID used in Ensembl, so I used the same name
      feature.locations.each do |location|
      start= location.from
      end_= location.to  
      strand = feature.assoc["strand"]
      attributes = "ID=#{feature.assoc.keys[0] + "_#{count}"};Name=#{feature.feature};"
      
      gff3_exon.write("#{seqid}\t#{source}\t#{type}\t#{start}\t#{end_}\t#{score}\t#{strand}\t#{phase}\t#{attributes}\n")
        
      end
    end 
  end
  gff3_exon.close
end

#Create a gff3 file with the positions of CTTCTT repeats on the chromosomes

  # @param new_embl[Hash] hash with the name of the gene as key and Bio::Sequence object with new features added
  # @return [gff3_chrom] a GFF3 file 
def get_gff3_chrom(new_embl)
#seqid - ID of the chromosome
#type - type of feature. Must be a term or accession from the SOFA sequence ontology
#start - Start position of the feature, with sequence numbering starting at 1.
#end - End position of the feature, with sequence numbering starting at 1.
#score - A floating point value.
#strand - defined as + (forward) or - (reverse).
#phase - . 
#attributes - A semicolon-separated list of tag-value pairs, providing additional information.
  gff3_chrom = File.open('gff3_chrom', mode: 'w')
  gff3_chrom.write("##gff-version 3\n") #the first line of the gff3 file
  type = "SO:0000657"
  source = "Ensembl"
  score = "."
  phase = "."
  
  new_embl.each do |key, info|
    info = info.to_biosequence
    count = 0
    info.features.each do |feature|
      next if feature.feature != "CTTCTT"
      count += 1
      seqid = info.primary_accession.split(":")[2] #the number of the chromosome
      feature.locations.each do |location|
      chr_start = info.primary_accession.split(":")[3] #start position of the chromosome
      
      #@We add to the positions the start position of the chromosome.
      #  Both positions start in 1 so we need to rest 1 in order to obtain the correct position
      start= location.from.to_i + chr_start.to_i - 1
      end_= location.to.to_i  + chr_start.to_i - 1
      strand = feature.assoc["strand"]
      attributes = "ID=#{feature.assoc.keys[0] + "_#{count}"};Name=#{feature.feature};"
      
      gff3_chrom.write("#{seqid}\t#{source}\t#{type}\t#{start}\t#{end_}\t#{score}\t#{strand}\t#{phase}\t#{attributes}\n")
        
      end
    end 
  end
  gff3_chrom.close
end

#Create a report file with the genes without CTT repeats

  # @param new_embl[Hash] hash with the name of the gene as key and Bio::Sequence object with new features added
  # @return [report] a report file
def report_no_CTT(new_embl)
  
  report = File.open('report_no_CTT', mode: 'w')
      report.write("*******************************************************************************\n")
      report.write("*******************************REPORT FILE*************************************\n")
      report.write("*******************************************************************************\n")
      report.write("*************These are the gene on the list without CTT repeats**************\n")
   new_embl.each do |key, info|
    info = info.to_biosequence
    count = 0
    info.features.each do |feature|
      next if feature.feature != "CTTCTT"
      count += 1
    end
       if count == 0 #when the count is 0, we add the key of the gene, because that means that it have 0 CTT repeats
         report.write("#{key}\n")
       end
    end
end


