library(ape)


#########################################
# Obtencion de parametros del gen IGH
#########################################
gen_IGH <- read.dna("gene1.txt",format ="fasta")
gen_IGH_1 = as.matrix(gen_IGH[1])

le_gen_IGH_1 = length(gen_IGH_1)#Longitud del gen

inigen1 = 105586437#Position en el cromosoma
fingen1 = 105865678#Position en el cromosoma

inicut1 = 105834401#Position en el cromosoma
finicut1 = 105865678#Position en el cromosoma

ini_seq = inicut1-inigen1+1##Posicion en el geen
fin_seq = finicut1-inigen1+1##Posicion en el geen



################# Parametros Gen IGH ############
ini_seq_rev = le_gen_IGH_1-ini_seq+1
fin_seq_rev = le_gen_IGH_1-fin_seq+1
le_gen_IGH_1 = length(gen_IGH_1)#Longitud del gen
#################################################
############### Estimaciones Gen IGH ############
#Numero de segmentos
num_seg_gen_IGH_1 = ini_seq_rev-fin_seq_rev+2

#Longitud maxima y minima de los segmentos
#Tipo de corte 1------Nucleotido de corte
min_seg_gen_IGH_1 = le_gen_IGH_1-ini_seq_rev
max_seg_gen_IGH_1 = le_gen_IGH_1-fin_seq_rev +1
min_seg_gen_IGH_1
max_seg_gen_IGH_1







######################################################
#####Union del Gen con parte del cromosma que incluye la region de corte##
######################################################
gen_CEBPD <- read.dna("gene2.txt",format ="fasta")#Gen
gen_CEBPD

crom_CEBPD = read.dna("Cromosoma.txt", format="fasta")#Cromosoma

#Apartir del nucleotido 9 inicia la parte que le faltaba al gen (region de corte)
seg.cromo = crom_CEBPD[9:(9+(61*10^3)-1)]
gen_CEBPD.seg.cromo = gen_CEBPD

gen_CEBPD.seg.cromo[(length(gen_CEBPD)+1):(length(gen_CEBPD)+length(seg.cromo))]=seg.cromo
gen_CEBPD.seg.cromo#Gen CEBPD junto con la region de corte

write.dna(as.matrix(gen_CEBPD.seg.cromo), file="jojo.txt", format = "fasta")
fio <- read.dna("jojo.txt",format ="fasta")
#######




### Parámetros del Gen CEBPD ######
Inicut_gen_CEBPD=length(gen_CEBPD)+(18*10^3)
Finicut_gen_CEBPD=length(gen_CEBPD.seg.cromo)
######
############### Estimaciones Gen CEBPD ############
#Numero de segmentos
num_seg_gen_CEBPD = Finicut_gen_CEBPD-Inicut_gen_CEBPD+2

#Longitud maxima y minima de los segmentos
#Tipo de corte 1------Nucleotido de corte
min_seg_gen_CEBPD = Inicut_gen_CEBPD-1
max_seg_gen_CEBPD = Finicut_gen_CEBPD

###############################
###Resultados Finales ##########
min_seg_gen_IGH_1+min_seg_gen_CEBPD#Longitud minima de los nuevos genes
max_seg_gen_IGH_1+max_seg_gen_CEBPD#Longitud máxima de los nuevos genes

num_seg_gen_CEBPD*num_seg_gen_IGH_1#Numero total de segmentos
################################





#### Save a sample of all retults in fasta format #####
set.seed(35)
#fin_seq_rev:(ini_seq_rev+1)
samp_cut_IGH1 = sample(fin_seq_rev:(ini_seq_rev+1),32)
#(Inicut_gen_CEBPD-1):Finicut_gen_CEBPD
samp_cut_CEBPD = sample((Inicut_gen_CEBPD-1):Finicut_gen_CEBPD,32)


for(k in samp_cut_IGH1){
  for(i in samp_cut_CEBPD){
    gen_IGH_1_cut_k = gen_IGH_1[k:le_gen_IGH_1]
    
    gen_CEBPD_cut_i = gen_CEBPD.seg.cromo[1:i]
    
    gen_CEBPD_IGH_1 = gen_CEBPD_cut_i
    gen_CEBPD_IGH_1[(length(gen_CEBPD_cut_i)+1):(length(gen_CEBPD_cut_i)+length(gen_IGH_1_cut_k))]=gen_IGH_1_cut_k
    write.dna(as.matrix(gen_CEBPD_IGH_1), file=paste("Gen_UnCEBPD-",i,"_IGH1-",k,".txt",sep=""),format="fasta")
  }
}

