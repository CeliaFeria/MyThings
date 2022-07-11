
 class Cross
  attr_accessor :Parent1
  attr_accessor :Parent2
  attr_accessor :F2_Wild
  attr_accessor :F2_P1
  attr_accessor :F2_P2
  attr_accessor :F2_P1P2
  @@linked = Hash.new
  
  def initialize (params={})
    @Parent1 = params.fetch(:Parent1, "ABC")
    @Parent2= params.fetch(:Parent2, "ABC")
    @F2_Wild= params.fetch(:F2_Wild, "123")
    @F2_P1 = params.fetch(:F2_P1, "123")
    @F2_P2= params.fetch(:F2_P2,"123")
    @F2_P1P2 = params.fetch(:F2_P1P2,"123")
  end
  
  def Cross.openfile(path)
  cross = File.open(path, mode: 'r')
  cross.readlines[1..-1].each do |line|
  parent1, parent2, f2_w, f2_p1, f2_p2, f2_p1p2 = line.strip.split("\t")
  Cross.new({:Parent1 => parent1,
    :Parent2 => parent2,
    :F2_Wild => f2_w,
    :F2_P1 => f2_p1,
    :F2_P2 => f2_p2,
    :F2_P1P2 => f2_p1p2})
    end
  cross.close

end

    
  #unlinked genes have the same probability to pass to F2, so in order to test if they are linked we suposse that they have a 9,3,3,1 proportion
  def Cross.chi_square(path)
  cross = File.open(path, mode: 'r')
  cross.readlines[1..-1].each do |line|
  @proportion= []
  @result = []
  parent1, parent2, f2_w, f2_p1, f2_p2, f2_p1p2 = line.strip.split("\t")
    f_obs = [f2_w,f2_p1,f2_p2,f2_p1p2].map{|x| x.to_f}
    @suma = f_obs.sum
    f_obs.each do |gene|
    f_exp = [9,3,3,1]
    f_exp.each do |expected|
    if f_obs.index(gene) == f_exp.index(expected)
      a = gene
      b = (expected*@suma)/16
     prp = ((a-b)**2)/b 
      #puts prp
     @proportion.push(prp)

    end 

    end

    end
    @result= @proportion.sum
    #puts result
       if @result > 7.8 #3 degree of freedom
        Stock_info.all_stocks.each do |x| 
          if x.SeedStock == parent1
          @gene1 = x.Gene_associated.GeneName
          end 
          end
        Stock_info.all_stocks.each do |x|
          if x.SeedStock == parent2
          @gene2 = x.Gene_associated.GeneName
          end
          end 
        
        puts "Recording: #{@gene1} is genetically linked to #{@gene2} with a chisquare #{@result} "
        @@linked[@gene1] = @gene2
        @@linked[@gene2] = @gene1
         
        puts "\n\nFinal Report:\n\n#{@gene1} is linked to #{@gene2}\n#{@gene2} is linked to #{@gene1}"
       end
  end
  end 
    
  def Cross.Linked_genes
    return @@linked
  end 
  
end


