
require "./Gene.rb"
require "./Stock.rb"
require "./Cross.rb"

Gene.openfile('./gene_information.tsv')

Stock_info.openfile('./seed_stock_data.tsv')
Stock_info.grams_remaining
Stock_info.update

Cross.openfile('./cross_data.tsv')
Cross.chi_square('./cross_data.tsv')


