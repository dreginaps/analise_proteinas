# Calcular massa média de sequências de peptídeos
# ------------------------------------------------------------------------------
# Código: Calcular massa média de sequências de peptídeos em arquivo .csv
# Criado por: Dávila Regina Pacheco Silva
# Data: 20/07/2025
# Objetivos:
# - Ler arquivo Fasta
# - Corrigir nome de coluna do arquivo, se necessário
# - Retirar stop códon das sequências (para arquivos provenientes de RNAseq)
# - Calcular as massas médias das sequências
# - Exportar resultado para arquivo .csv
# - Agrupar massas em faixas específicas para gerar de gráficos
# - Contar sequências por agrupamentos
# - Gerar dataframe com os dados
# - Gerar histograma
# - Gerar gráfico de densidade
# - Gerar gráfico de barras
# - Exportar gráficos para .svg e .jpeg
#
# Pré-requisitos:
# - Os dados devem estar em um arquivo em formato .csv
# - A identificação das sequências deve estar em coluna nomeada como "ID"
# - As sequências devem estar em coluna nomeada como "Sequencia"
# - A pasta com o arquivo .csv DEVE estar na mesma pasta onde este R Script está
# 
# ------------------------------------------------------------------------------


# Pacotes necessários

library(rstudioapi)
library(Peptides)
library(dplyr)
library(ggplot2)


# Alocar o arquivo na pasta onde está salvo 

current_path = getActiveDocumentContext()$path
setwd(dirname(current_path))
print(getwd())

# Importar o arquivo .csv, considerando ";" como separadores e "," como símbolo de casa decimal

sp_massas <- read.csv("sp_massas.csv",, sep=";", dec=",")
head(sp_massas)

# Renomenado coluna com caractere especial (descomente e ajuste, caso necessário)

#sp_massas <- sp_massas %>% 
#  rename(ID = X..ID)
#head(sp_massas)


# Retirando sinal de stop códon das sequências (*):

sp_massas <- sp_massas %>%
  mutate(sequencia_limpa = sub("\\*$", "", Sequencia))

head(sp_massas$sequencia_limpa)


# Calculando as massas médias das sequências utilizando mw() do pacote Peptides:

sp_massas <- sp_massas %>%
  mutate(
    massa_Da = sapply(sequencia_limpa, function(seq)
      round(
        mw(
          seq = seq,
          monoisotopic = FALSE,
          avgScale = "expasy",
          label = "none",
          aaShift = NULL
        ),
        0  # arredondamento para 0 casas decimais
      )
    )
  )


head(sp_massas)


# Salvando para o arquivo "resultado_com_massa":

write.csv(sp_massas, "resultado_com_massa.csv")



# Agrupando massas em faixas específicas

# Importar o dado contendo o dado de massa
#read.csv("nomedoarquivo.csv")


# Entre 7 kDa e 11 kDa (não inclui 11k):

btwn_sete_onze <- sp_massas %>%
  filter(Massa >= 7000 & Massa < 11000)

# Entre 11 kDa e 15 kDa (não inclui 15k):

btwn_onze_quinze <- sp_massas %>%
  filter(Massa >= 11000 & Massa < 15000)

# Entre 15 kDa e 20 kDa (não inclui 20k):

btwn_quinze_vinte <- sp_massas %>%
  filter(Massa >= 15000 & Massa < 20000)

# Entre 20 kDa e 50 kDa (não inclui 50k):

btwn_vinte_cinq <- sp_massas %>%
  filter(Massa >= 20000 & Massa < 50000)

# Entre 50 kDa e 100 kDa (não inclui 100k):

btwn_cinq_cem <- sp_massas %>%
  filter(Massa >= 50000 & Massa < 100000)

# Acima de 100 kDa:

abv_cem <-  sp_massas %>%
  filter(Massa >= 100000)

# Contando linhas de cada agrupamento (quantas observações de transcritos em cada faixa):
# As somas de transcritos por faixa são atribuídas a variáveis

trnsc_sete_a_onze <- nrow(btwn_sete_onze)
trnsc_onze_a_quinze <- nrow(btwn_onze_quinze)
trsnc_quinze_a_vinte <- nrow(btwn_quinze_vinte)
trnsc_vinte_a_cinq <- nrow(btwn_vinte_cinq)
trnsc_cinq_a_cem <- nrow(btwn_cinq_cem)
trnsc_acima_cem <- nrow(abv_cem)

# Criando um novo dataframe com os dados:

tb_faixas <- tibble(faixa = c("7k - 11k", "11k - 15k", "15k - 20k", "20k - 50k", "50k - 100k", "100k +"),
                    transcritos = c(trnsc_sete_a_onze, trnsc_onze_a_quinze, trsnc_quinze_a_vinte, trnsc_vinte_a_cinq, trnsc_cinq_a_cem, trnsc_acima_cem))

# Transformando as faixas em fatores (números), para ficarem na ordem correta:

tb_faixas <- tb_faixas %>%
  mutate(
    faixa = factor(
      faixa,
      levels = c("7k - 11k", "11k - 15k", "15k - 20k",
                 "20k - 50k", "50k - 100k", "100k +")
    )
  )

tb_faixas

# Histograma com os dados não agrupados em faixas de massas:

histogram_graph <- ggplot(sp_massas, aes(x = Massa)) +
  geom_histogram(binwidth = 1000, fill = "salmon", color = "white") +
  labs(x = "Massa (Da)", y = "Frequência") +
  theme_minimal()


density_graph <- ggplot(sp_massas, aes(x = Massa)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(x = "Massa (Da)", y = "Densidade estimada") +
  theme_minimal()

bar_graph <- ggplot(tb_faixas, aes(x = faixa, y = transcritos, fill = faixa)) +
  geom_bar(stat = "identity") +
  labs(x = "Faixa de massa", y = "Nº de transcritos") +
  scale_fill_brewer(palette = "Dark2")
theme_minimal()

# Salvando gráficos em .svg:

ggsave("faixa_proteinas_histograma.svg", histogram_graph)
ggsave("faixa_proteinas_densidade.svg", density_graph)
ggsave("faixa_proteinas_barras.svg", bar_graph)

# Salvando gráficos em .jpeg:

ggsave("Fig1_faixa_proteinas_histograma.jpeg", histogram_graph)
ggsave("Fig2_faixa_proteinas_densidade.jpeg", density_graph)
ggsave("Fig3_faixa_proteinas_barras.jpeg", bar_graph)