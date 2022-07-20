
require 'rest-client' 

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



def get_AGI_Locus(path)
  locus_code= []
  locus = File.open(path, mode: 'r')
  locus.readlines[1..60].each do |line| #I tried the code only with 60 genes
    code = line.strip.split("\n")
  locus_code |= code
  end 
  return locus_code
end  


def get_interacctions(locus)
    res = fetch("http://bar.utoronto.ca:9090/psicquic/webservices/current/search/interactor/#{locus}/?firstResult=0&maxResults=30&format=tab25");  #restritions: 9090-> species, maxresult -> 30
    body = res.body.split("\n") #each interaction is separated by \n
    all_interacctions = []
    all_int_gen = Hash.new
    body.each do |elem|
      elem = elem.split("\t")
      elem[-1].to_s =~ /(\d.\d+)/ #select score
      score = $1
      if score.to_f > 0.5  #significative medium-high interaction above 0.4
        elem[2] =~ /(A[Tt]\d[Gg]\d\d\d\d\d)/ 
        gen1 = $1
        elem[3] =~ /(A[Tt]\d[Gg]\d\d\d\d\d)/
        gen2 = $1
        next if gen1.nil?||gen2.nil?
        #puts gen1,gen2

        if gen1.upcase != locus.upcase #I use upcase because all locus have a T instead of the t in the gen interactors
          all_interacctions.push(gen1) #add to all_interactions the genes that interact with "locus gene"
        else
          all_interacctions.push(gen2)
        end 
      end
    end
    if all_interacctions[0].nil? == false
    all_int_gen[locus] = all_interacctions #associate the interactions with the "locus gene"
    #puts all_int_gen
    end
  return all_int_gen
end





def get_complex_interactions(all_locus)
interactions = Hash.new
int_2 = []
all_locus.each do |locus|
  int_1 = get_interacctions(locus)
  #print("***********INTERACCIONES 1*************************")
  #puts int_1
  #print("************ INTERACCIONES 2 ********************")
  

  int_1.each do |key, values|
     include_int = []
     new_int = []
    values.each do |value|
      int_2.push(get_interacctions(value)) #interactions of the intermediate genes
    end
     #puts int_2
  
    int_2.each do |x|
    x.each do |key2, values2|
      #info = []
      all_locus = all_locus.map{|x| x.upcase} 
        if all_locus.include?(key2.upcase) #if key2 is include in all_locus that means that there is a relation between two genes on the list
          if key2.upcase != key.upcase  
          #print("included another gen on the list\n")
          interactions[key] = key2
          end
        end
    
      values2.each do |value|
        if all_locus.include?(value.upcase)
          if value.upcase != key2.upcase && value.upcase != key.upcase #eliminate repited interactions
          #print("included\n")
         # next if include_int.include?(value)
          #include_int.push(value)
            next if new_int.include?(value)
            new_int.push(key2)
            new_int.push(value)
          end
      end
    
    end
            next if new_int[0].nil?
            interactions[key]= new_int
            #puts "\ngetting interacctions"
            #puts interacctions
            #puts "---------"
        end
  end
end
  end
 return interactions
end


    



def get_Network(interactions)
  aux_has = Hash.new
  aux_has = interactions #auxiliar hash because without it I get an error, i cant modify the hash and iterate
  network = Hash.new #new hash with the network information
  interactions_keys = interactions.keys.map{|x| x.upcase}
  array = []
  aux_has.each do |key, values2|
    count = 0
      values2.each do |value2|
      if interactions_keys.include?(value2.upcase) #include the information of the "locus gene" that interact with other "locus gene"
        next if array.include?(value2)
        array.push(value2)
        bigger_int = interactions.keys.select{|x| x.upcase == value2.upcase} 
         next if interactions[bigger_int[0]].nil? 
         interactions[bigger_int[0]].each do |x| 
           if interactions_keys.include?(x.upcase) #include only the interactions between genes on the list
           next if array.include?(x)
           array.push(x)
           end
         end
         network[key] = array #include all the interactions associated to a "locus gene"
         next if network[bigger_int[0]].nil?
         network[bigger_int[0]] = NilClass
         network.reject! { |key, value| value == NilClass } #eliminate the key and values of the "locus gene" included in other list if it is on the network hash
        count = count +1
      end 
      end 
    if count == 0 
      network[key] = values2 # if there is no conexion between "locus gene" only had the information without addition 
    end
  end
  #print("interacciones final")
   #puts interactions
  return network
end

def get_KEGG_info (locus)
  kegg_info = Hash.new
  locus.each do |locus|
  res = fetch("http://togows.org/entry/genes/ath:#{locus}/pathways.json");  #ath = Arabidopsis thaliana
  body = JSON.parse(res.body)
  body[0].each do |key, values|
    next if values.nil?
    next if kegg_info.keys.include?(key)
    kegg_info[key] = values
  end
  end
    return kegg_info
end 


def get_GO (locus)
  go_info = []
  locus.each do |locus|
  res = fetch("http://togows.org/entry/uniprot/#{locus}/dr.json"); 
  #body = res.body
  body = JSON.parse(res.body)
  #puts body
  for elem in body[0]["GO"].each
    if elem[1] =~ /^P:/
      go = elem[1]
      next if go_info.include?(go)
    go_info.push(go)
    end
  end
  end
  return go_info
end 
     


