
require "./Web.rb"

class Network
  
attr_accessor :Interactions
attr_accessor :GO_info
attr_accessor :KEGG_info 

#@@Interaction_Network = []
@interacctions = Hash.new
  
  

  
def initialize (params={})
    @Interactions = params.fetch(:Interactions, "ABC")
    @GO_info = params.fetch(:GO_info, "ABC")
    @KEGG_info = params.fetch(:KEGG_info, "ABC")
  
end
  
  
def Network.information(info)
  info.each do |key, value|
  all = Network.new({:Interactions => info[key], 
    :GO_info => get_GO(value), 
    :KEGG_info =>  get_KEGG_info(value)})
      
    @interacctions[key] = all
  end
end 
  
  def Network.all_info
    return @interacctions
  end
  
  def Network.report
  report = File.open('report', mode: 'w')
      report.write("*******************************************************************************\n")
      report.write("*******************************REPORT FILE*************************************\n")
      report.write("*******************************************************************************\n")
      
      @interacctions.each do |key, values|
        report.write("\n\nInteractions\n")
        report.write("#{key}\t")
        @interacctions[key].Interactions. each do |int|
        report.write("#{int}\t")
        end
        
        report.write("\n\nGO information\n")
        @interacctions[key].GO_info.each do |go|
        report.write("#{go}\t")
        end
        
        report.write("\n\nKEGG information\n")
        @interacctions[key].KEGG_info.each do |key,value|
        report.write("ID:#{key}\tName:#{value}\n")
        end
      end 
    report.close
  end
  
  
end 


