=begin
Explanation of my code:

I retreived the positions of exons in the gene as well as their global position in genome in one script, putting them in separate files.

I queried for "CTTCTT" sequence on both plus and minus strands and had it in mind when retreiving the global positions of the repeat.

In case of the minus strand, first, the sequence was reversed-complemented and CTTCTT was queried. The inicial position of the repeat would correspond
to the end of the chromosome minus the position where the repeat "ends" in the gene. The end of the repeat would be the previous plus the lenght of
the repeat.

I tried to make the code in such a way that only changing the variable "repeat" with the sequence to query, it would work universaly

I turned off all the puts statements that were guiding me when writting the code and just put the counter
=end


require 'bio'
require 'net/http' 

#Teachers code for function fetching stuff
def fetch(uri_str)  # this "fetch" routine does some basic error-handling.  

  address = URI(uri_str)  # create a "URI" object (Uniform Resource Identifier: https://en.wikipedia.org/wiki/Uniform_Resource_Identifier)
  response = Net::HTTP.get_response(address)  # use the Net::HTTP object "get_response" method
                                               # to call that address

  case response   # the "case" block allows you to test various conditions... it is like an "if", but cleaner!
    when Net::HTTPSuccess then  # when response is of type Net::HTTPSuccess
      # successful retrieval of web page
      return response  # return that response object
    else
      raise Exception, "Something went wrong... the call to #{uri_str} failed; type #{response.class}"
      # note - if you want to learn more about Exceptions, and error-handling
      # read this page:  http://rubylearning.com/satishtalim/ruby_exceptions.html  
      # you can capture the Exception and do something useful with it!
      response = false
      return response  # now we are returning False
    end 
end





#Variables

no_repeat = [] #Genes without cttctt repeat in exons
repeat ="cttctt" #Repeat sequence to query
has_repeats = File.new("has_repeats.gff", "w")

list = File.open('gene_list.txt', 'r')
list=list.read.split()

#Here will be the chromosomes gff
chromosomes_file = File.new("chromosomes_file.gff", "w")

i=1

list.each do |gene| #loop of genes. Retreive info
  address = URI("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene}")
  response = fetch(address)
  record = response.body
  entry = Bio::EMBL.new(record) #Create EMBL object
  bioseq = entry.to_biosequence
  
  true_position_p = [] #position of repeat in the exon plus strand
  true_position_m = [] #position of repeat in the exon minus strand
  
  chromosomes= []
  
  #In the database, chromosome coordinates are in the "ac" field and they have the following structure: chromosome:TAIR10:n:ci:ce:1
  #Where n is number of chromosome, ci is the begining of the gene and ce is the end. So I split this information in separate arrays
  chromosomes=(entry.ac[0]).split(":")
  #print "#{entry.ac[0]}\n"
  name_chr = chromosomes[1] #chromosome name
  number_chr = chromosomes[2] #number
  begin_chr = chromosomes[3]#start
  end_chr = chromosomes[4]  #end
  
  #Loop the exons, retreive position
  entry.features.each do |feature|
    if feature.feature == "exon"
      hash = feature.assoc.to_s
      position = feature.position.to_s
      
  
      
      #Select exons
      if /exon_id/.match(hash)
          unless /[A-Z]/.match(position)
            next unless feature.feature =="exon"
            #Improve the annotation. The information of exons is stored similarly as in chromosome
            a = position.tr('complement()','')
            a = a.split("..")
            begining = a[0].to_i #Initial position of the exon
            finish =a[1].to_i #Final position of the exon
            
            #If strand is minus, perform reverse_complement of bioruby
            if /complement/.match(position)  
              bioseq_rev = bioseq.reverse_complement
              exon = bioseq_rev.subseq(begining,finish)
              #exon = bioseq_rev [begining..finish]
              
              strand = "minus"
            else 
              exon = bioseq.subseq(begining,finish)
              strand = "plus"
            end
            
            #Look for the repeats.
            #repeat_position_exon = (0 ... exon.length).find_all { |i| exon[i,repeat.length].match repeat }
            repeat_position_exon = (0 ... exon.length).find_all { |i| exon[i,repeat.length]==repeat }

            repeat_position_gene = []
            repeat_position_exon.each do |initial_position|
              aux = []
              initial_position = initial_position + begining + 1 #because in informatics the 0 is not the same as in informatics
              final_position = initial_position + repeat.length - 1 
              aux.push initial_position
              aux.push final_position 
              repeat_position_gene.push aux
            
            end
            if strand == "plus"
                repeat_position_gene.each do |coor|
                true_position_p.push coor
              end
            else
              repeat_position_gene.each do |coor|
              true_position_m.push coor
              end
            end
            true_position_p = true_position_p.uniq  
            true_position_m = true_position_m.uniq 
          end 
      end
    end
  end
  
  if true_position_p.any? or true_position_m.any?
     
    #Create like in Lesson7
    true_position_p.each do |begin_end|
     
      f1 = Bio::Feature.new('myrepeat',begin_end)
      f1.append(Bio::Feature::Qualifier.new('repeat_motif', 'CTTCTT'))
      f1.append(Bio::Feature::Qualifier.new('strand', '+'))
      bioseq.features << f1 
    end
    true_position_m.each do |begin_end|
      f2 = Bio::Feature.new('myrepeat',begin_end)
      f2.append(Bio::Feature::Qualifier.new('repeat_motif', 'CTTCTT'))
      f2.append(Bio::Feature::Qualifier.new('strand', '-'))
      bioseq.features << f2 
    end

#Write GFF
    entry.features.each do |feature|
      hash = feature.assoc
      position = feature.position
      #next unless feature.feature=="exon"
      if /repeat_motif/.match(hash.to_s)
        #puts hash
        begining = position[0].to_i #Initial position
        finish =position[1].to_i #Final position 
        if /"-"/.match(hash.to_s) # -strand
          #print "#{gene}\t.\t#{repeat}\t#{begining}\t#{finish}\t.\t-\t.\n"
          has_repeats.puts"#{gene}\t.\t#{repeat}\t#{begining}\t#{finish}\t.\t-\t.\n"
        else # +strand
          #print "#{gene}\t.\t#{repeat}\t#{begining}\t#{finish}\t.\t+\t.\n"
          has_repeats.puts "#{gene}\t.\t#{repeat}\t#{begining}\t#{finish}\t.\t+\t.\n"
        end
        
        #Same but for chromosomes
        #Here I added some extra . because GFF format for loading in ENSEMBL has to be very specific, there should be 9 columns in GFF to visualize it. And points show "empty" fieds
        #Those plus and minus 1 is because there was an outset, probably because nucleotide 1 is 1, not 0 like in programming.
        if /"-"/.match(hash.to_s)
          chr_minus_ini = end_chr.to_i - finish +1
          chr_minus_end =  chr_minus_ini + repeat.length + (-1)
          #print "#{name_chr}:#{number_chr}\t#{gene}\trepeat_region\t #{chr_minus_ini}\t #{chr_minus_end} \t-\t.\t.\n\n"
          chromosomes_file.puts "#{number_chr}\t.\trepeat_region\t#{chr_minus_ini}\t#{chr_minus_end}\t.\t-\t.\t.\n"
        else
          chr_plus_ini = begining.to_i + begin_chr.to_i - 1
          chr_plus_end = chr_plus_ini + repeat.length 
          #print "#{name_chr}:#{number_chr}\t#{gene}\trepeat_region\t #{chr_plus_ini}\t #{chr_plus_end} \t+\t.\t.\n\n"
          chromosomes_file.puts "#{number_chr}\t.\trepeat_region\t#{chr_plus_ini}\t#{chr_plus_end}\t.\t+\t.\t.\n"
        end
      end
    end
  end
    
  
#For those genes without our query repeat in the exons, I put them in another file
  if true_position_p.empty? and true_position_m.empty?
    no_repeat.push gene
  end
i+=1
puts "#{i}\n"
end

#Output files with no target repeat in exons
has_no_repeats = File.new("no_repeat", "w")
has_no_repeats.print no_repeat.to_s