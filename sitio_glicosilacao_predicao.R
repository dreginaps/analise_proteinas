# ------------------------------------------------------------------------------
# Código: Predição de glicosilação para arquivo .fasta
# Criado por: Dávila Regina Pacheco Silva
# Data: 24/07/2025
# Objetivos:
# - Criar função para detectar motivo N[^P][ST] (N-Glicosilação)
# - Ler arquivo Fasta
# - Aplicar a função às sequências contidas no arquivo Fasta
# - Gerar tabela contendo dados
# - Exportar tabela para arquivo .csv
#
# Pré-requisitos:
# - Os dados devem estar em um arquivo em formato .fasta
# - A pasta com o arquivo .fasta DEVE estar na mesma pasta onde este R Script está
#
# Outras possibilidades de análises:
# https://www.nature.com/articles/srep34595
# https://services.healthtech.dtu.dk/services/NetNGlyc-1.0/
# 
# ------------------------------------------------------------------------------



# Utilizar pacote do CRAN seqinr
# Instalar e carregar pacote:

install.packages("seqinr")
library(seqinr)

# Alocar o arquivo na pasta onde está salvo: 

current_path = getActiveDocumentContext()$path
setwd(dirname(current_path))
print(getwd())

# Criar função para detectar motivo N[^P][ST]:

detect_nglyco <- function(seq) {
  seq_str <- paste(seq, collapse = "")
  matches <- gregexpr("N[^P][ST]", seq_str)[[1]]
  if (matches[1] == -1) return(integer(0))
  return(matches)
}

# ler arquivo .fasta:
# caso o nome do arquivo esteja diferente, ajustar no código

fasta_file <- "seq_proteins_verif_glic.fasta"
seqs <- read.fasta(fasta_file, seqtype = "AA")

# aplicar a função para as sequências do arquivo .fasta

results <-  lapply(names(seqs), function(name) {
  pos <- detect_nglyco(seqs[[name]])
  data.frame(
    SequenceID = name,
    N_Glyco_Positions = paste(pos, collapse = ", ")
  )
})

# gerar tabela:

final_df <- do.call(rbind, results)

# visualizar tabela:

print(final_df)


# exportar para csv:

write.csv(final_df, "nglyco_prediction_results.csv", row.names = FALSE)
