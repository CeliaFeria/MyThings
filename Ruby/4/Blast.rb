
#This code search for orthologs between two databases

require 'bio'



#Read the type of information in order to define de database type

  # @param data[String] the path of the file where is the information that is going to be used to create the database
  # @return [data_type] the type of database
def database_type(data)
  #We have different type of data, protein data for S.pombe and nucleotids for Arabidopsis
  database = Bio::FlatFile.auto(data)
  entry = database.next_entry.to_biosequence.guess #read the first entry type, that can be Bio::Sequence::NA, or Bio::Sequence::AA
  
  if entry==Bio::Sequence::NA
    data_type='nucl'
  else 
    data_type='prot'
  end
  
  return data_type
end 



#Create a database with the information present in a file

  # @param data[String] the path of the file where is the information that is going to be used to create the database
  # @param data_type[String] the type of database that is going to be created
def create_database(data, data_type)
  makedatabase = "makeblastdb -in #{data} -dbtype #{data_type} -out #{data}"
  system(makedatabase)
end

#Read the type of information in database and the query file in order to return the blast type

  # @param database[String] the path of the file where is the information that is going to be used as database
  # @param query[String] the path of the file where is the information that is going to be used as query
  # @return [type] the type of blast that is going to be performed
def blast_type(database, query)
    #we only have this two cases
    if database_type(database) =='prot' && database_type(query) =='nucl'
    type='blastx'  #translate the nucleotide query in order to search in a protein database
    elsif database_type(database) =='nucl' && database_type(query) =='prot'
      type='tblastn'#translate the database
    end
 return type
end




#Perform the blast

  # @param blast_type[String] the type of blast that will be perform
  # @param database1[String] the path of the database
  # @param file_database2[String] the path of the file where is the information that is going to be used as query
  # @return [best_hits] a hash with the best hit for each query when the e-value is less than 1e-6
def blast(blast_type1, database1,file_database2)
  best_hits = Hash.new
  factory1 = Bio::Blast.local(blast_type1, database1)
  
  query = Bio::FlatFile.auto(file_database2)

  query.each_entry do |entry|

  report = factory1.query(entry)
  next if report.hits.empty?
  next if report.hits[0].evalue > 1e-6 #https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0101850
    prot1 = entry.definition.split("|").first.strip #select only the name of the gene 
    prot2 = report.hits[0].definition.split("|").first.strip  
    
    best_hits[prot1] = prot2
end 
  
  return best_hits
end 


#Perform the blast

  # @param blast1[Hash] a hash with the best hits obtained in the first blast
  # @param blast2[Hash] a hash with the best hits obtained in the second blast
  # @return [best_reciprocal] a hash with the orthologs 
def orthologs(blast1, blast2)
  best_reciprocal = Hash.new
  blast1.each do |key,value|
    blast2.each do |key2,value2|
      if key == value2 && value == key2 
        best_reciprocal[key] = value
      end
    end
  end
  return best_reciprocal
end




#Create a report file with the orthologs

  # @param orthologs[Hash] hash with the name of the orthologs 
  # @return [report] a report file
def report(orthologs)
  report = File.open('report_orthologs', mode: 'w')
      report.write("*******************************************************************************\n")
      report.write("*******************************REPORT FILE*************************************\n")
      report.write("*******************************************************************************\n")
      report.write("*************These are the gene on the list with reciporcla hits***************\n\n")
      
  
      orthologs.each do |key,values|
      report.write("\t\t\t#{key}\t\t#{values}\n")
      end
  report.close
end


