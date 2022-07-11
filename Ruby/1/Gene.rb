
#require "./Cross.rb"

class Gene
  attr_accessor :GeneID 
  attr_accessor :GeneName
  attr_accessor :MutantPhenotype
  attr_accessor :LinkedGenes
  @@genes = []
  
  def initialize (params={})
    @GeneID = params.fetch(:GeneID)
    @GeneName= params.fetch(:GeneName)
    @MutantPhenotype= params.fetch(:MutantPhenotype)
    @LinkedGenes = Cross.Linked_genes.find{|key, value| key == @GeneName} #add information about the linked genes
    
    @@genes << self
  end
  
    
  def Gene.openfile(path)
  gene = File.open(path, mode: 'r')
  gene.readlines[1..-1].each do |line|
  geneID, geneName, phenotype = line.strip.split("\t")
  if geneID =~ /A[Tt]\d[Gg]\d\d\d\d\d/
    Gene.new({:GeneID => geneID,
    :GeneName => geneName,
    :MutantPhenotype => phenotype})
    
  else 
    puts("The gene ID is incorrect")
    break
  end 
  end
  gene.close
  end
  
   def Gene.all_genes 
    return @@genes
  end

end


