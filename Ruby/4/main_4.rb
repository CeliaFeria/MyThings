
require "./Blast.rb"

database1= create_database('./Arabidopsis', database_type('./Arabidopsis'))
database2 = create_database('./Spombe', database_type('./Spombe'))

file2 = './Spombe'
file1 = './Arabidopsis'
blast1 = blast(blast_type('./Arabidopsis','./Spombe'), './Arabidopsis', file2)
blast2 = blast(blast_type(file2,file1), './Spombe', file1)

ort = orthologs(blast1,blast2)

report(ort)


