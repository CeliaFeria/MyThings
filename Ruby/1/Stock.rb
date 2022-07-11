
class Stock_info  

  attr_accessor :SeedStock  
  attr_accessor :MutantGeneID
  attr_accessor :LastPlant
  attr_accessor :Room
  attr_accessor :StorageGrams
  attr_accessor :Gene_associated
  
   @stocks = Hash.new
   @@all_stocks = []
  
  def initialize(params={})
    @SeedStock = params.fetch(:SeedStock,  "ABC")
    @MutantGeneID = params.fetch(:MutantGeneID,  "ABC")
    @LastPlant = params.fetch(:LastPlant,  "ABC")
    @Room = params.fetch(:Room,  "CAM")
    @StorageGrams = params.fetch(:StorageGrams,  123)
    @StorageGrams = @StorageGrams.to_f
    @Gene_associated = Gene.all_genes.find{|x| x.GeneID == @MutantGeneID}
    

    @@all_stocks << self
  end
  
  def Stock_info.openfile (path)
    stock_file = File.open(path, mode: 'r')
    
    @title = stock_file.readline
    
    stock_file.readlines[0..-1].each do |line|
    seed, mutant, plant, place, storage = line.strip.split("\t")
    stock = Stock_info.new({:SeedStock => seed,
      :MutantGeneID => mutant,
      :LastPLant => plant,
      :Room => place,
      :StorageGrams => storage})
      
    @stocks[seed] = stock
    end
    stock_file.close
    end

    def Stock_info.all_stocks
      return @@all_stocks
    end

    def Stock_info.grams_remaining #update information
    @stocks.each do |key,values|
    @stocks[key].StorageGrams -= 7
    @stocks[key].LastPlant = Time.now.strftime('%d/%m/%Y')
      if @stocks[key].StorageGrams <= 0.0
        @stocks[key].StorageGrams = 0.0
        puts "WARNING: we have run out of Seed Stock #{key}"
      end
    end
end
  
  def Stock_info.update #create new file with update information
    f_new = File.open('new', mode: 'w')
    f_new.write(@title)
    @stocks.each do |key,values|
    f_new.write("#{@stocks[key].SeedStock}\t#{@stocks[key].MutantGeneID}\t#{@stocks[key].LastPlant}\t#{@stocks[key].Room}\t#{@stocks[key].StorageGrams}\n")
    end
    f_new.close
   end

end




