
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CALCULATE TOTAL EXPENDITURE ON CONSUMPTION OF FAMILIES RELATED TO AGRICULTURAL ACTIVITIES BY STATE
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: EXTRACT AND AGGREGATE THE STATISTIC OF INTEREST DIRECTLY ON THE RAW2CLEAN FILE TO AVOID STORING UNNECESSARY FILES
# 2: SCRIPT ADAPTED FROM "Leitura dos Microdados - R.R" and "Tabela de Despesa Geral.R" STORED ON DOCUMENTATION FOLDER AND DOWNLOADED FROM IBGE





# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "pof_raw2clean script", log = T)

# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(readxl) # to read excel files
library(tidyverse)  # manipulate tables, works with sf
library(sjlabelled) # label columns, preferred than Hmisc::label because has function to clear labels when necessary
library(sqldf) # to use SQL select on data frames





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

tradutor_despesa <-
  readxl::read_excel("data/raw2clean/pof_ibge/input/Tradutor_Despesa_Geral.xls")

if(!file.exists("data/raw2clean/pof_ibge/input/MORADOR.rds")) {

  # read input files
  morador <-
    read.fwf("data/raw2clean/pof_ibge/input/MORADOR.txt"
             , widths = c(2,4,1,9,2,1,2,2,1,2,2,4,3,1,1,
                          1,1,1,2,1,2,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,2,1,1,2,1,1,2,1,1,1,
                          2,1,2,14,14,10,1,1,20,20,20,20)
             , na.strings=c(" ")
             , col.names = c("UF", "ESTRATO_POF", "TIPO_SITUACAO_REG",
                             "COD_UPA", "NUM_DOM", "NUM_UC", "COD_INFORMANTE",
                             "V0306", "V0401", "V04021", "V04022", "V04023",
                             "V0403", "V0404", "V0405", "V0406", "V0407",
                             "V0408", "V0409", "V0410", "V0411", "V0412",
                             "V0413", "V0414", "V0415", "V0416",
                             "V041711", "V041712", "V041721", "V041722",
                             "V041731", "V041732", "V041741", "V041742",
                             "V0418", "V0419", "V0420", "V0421", "V0422",
                             "V0423", "V0424", "V0425", "V0426", "V0427",
                             "V0428", "V0429", "V0430", "ANOS_ESTUDO",
                             "PESO", "PESO_FINAL", "RENDA_TOTAL",
                             "INSTRUCAO", "COMPOSICAO", "PC_RENDA_DISP",
                             "PC_RENDA_MONET", "PC_RENDA_NAO_MONET",
                             "PC_DEDUCAO")
             , dec="."
    )

  # save input in a format that is faster to read in R
  saveRDS(morador, "data/raw2clean/pof_ibge/input/MORADOR.rds")

  } else {

  # read input
    morador <- readRDS("data/raw2clean/pof_ibge/input/MORADOR.rds")
}


if(!file.exists("data/raw2clean/pof_ibge/input/DESPESA_COLETIVA.rds")) {

  despesa_coletiva <-
    read.fwf("data/raw2clean/pof_ibge/input/DESPESA_COLETIVA.txt"
             , widths = c(2,4,1,9,2,1,2,2,7,2,4,10,2,2,1
                          ,10,1,12,10,10,1,1,2,14,14,10,5)
             , na.strings=c(" ")
             , col.names = c("UF", "ESTRATO_POF", "TIPO_SITUACAO_REG",
                             "COD_UPA", "NUM_DOM", "NUM_UC", "QUADRO",
                             "SEQ", "V9001", "V9002", "V9005", "V8000",
                             "V9010", "V9011", "V9012", "V1904",
                             "V1905", "DEFLATOR", "V8000_DEFLA",
                             "V1904_DEFLA", "COD_IMPUT_VALOR",
                             "COD_IMPUT_QUANTIDADE", "FATOR_ANUALIZACAO",
                             "PESO", "PESO_FINAL", "RENDA_TOTAL","V9004")
             , dec="."
    )

    # save input in a format that is faster to read in R
    saveRDS(despesa_coletiva, "data/raw2clean/pof_ibge/input/DESPESA_COLETIVA.rds")

  } else {

  # read input
  despesa_coletiva <- readRDS("data/raw2clean/pof_ibge/input/DESPESA_COLETIVA.rds")
}

if(!file.exists("data/raw2clean/pof_ibge/input/CADERNETA_COLETIVA.rds")) {

  caderneta_coletiva <-
    read.fwf("data/raw2clean/pof_ibge/input/CADERNETA_COLETIVA.txt"
             , widths = c(2,4,1,9,2,1,2,3,7,2,10,12,10,1,2,14,14,10,
                          9,4,5,9,5
             )
             , na.strings=c(" ")
             , col.names = c("UF", "ESTRATO_POF", "TIPO_SITUACAO_REG",
                             "COD_UPA", "NUM_DOM", "NUM_UC", "QUADRO",
                             "SEQ", "V9001", "V9002", "V8000", "DEFLATOR",
                             "V8000_DEFLA", "COD_IMPUT_VALOR",
                             "FATOR_ANUALIZACAO", "PESO", "PESO_FINAL",
                             "RENDA_TOTAL",
                             "V9005", "V9007", "V9009", "QTD_FINAL","V9004")
             , dec="."
    )
  # save input in a format that is faster to read in R
  saveRDS(caderneta_coletiva, "data/raw2clean/pof_ibge/input/CADERNETA_COLETIVA.rds")

  } else {

  # read input
  caderneta_coletiva <- readRDS("data/raw2clean/pof_ibge/input/CADERNETA_COLETIVA.rds")
}


if(!file.exists("data/raw2clean/pof_ibge/input/DESPESA_INDIVIDUAL.rds")) {

  despesa_individual <-
    read.fwf("data/raw2clean/pof_ibge/input/DESPESA_INDIVIDUAL.txt"
             , widths = c(2,4,1,9,2,1,2,2,2,7,2,10,2
                          ,2,1,1,1,12,10,1,2,14,14,10,5)
             , na.strings=c(" ")
             , col.names = c("UF", "ESTRATO_POF", "TIPO_SITUACAO_REG",
                             "COD_UPA", "NUM_DOM", "NUM_UC",
                             "COD_INFORMANTE", "QUADRO", "SEQ", "V9001",
                             "V9002", "V8000", "V9010", "V9011", "V9012",
                             "V4104", "V4105", "DEFLATOR", "V8000_DEFLA",
                             "COD_IMPUT_VALOR", "FATOR_ANUALIZACAO",
                             "PESO", "PESO_FINAL", "RENDA_TOTAL","V9004")
             , dec="."
    )
  # save input in a format that is faster to read in R
  saveRDS(despesa_individual, "data/raw2clean/pof_ibge/input/DESPESA_INDIVIDUAL.rds")

  } else {

  # read input
    despesa_individual <- readRDS("data/raw2clean/pof_ibge/input/DESPESA_INDIVIDUAL.rds")
}

if(!file.exists("data/raw2clean/pof_ibge/input/ALUGUEL_ESTIMADO.rds")) {

  aluguel_estimado <-
    read.fwf("data/raw2clean/pof_ibge/input/ALUGUEL_ESTIMADO.txt"
             , widths = c(2,4,1,9,2,1,2,7,2,10,2,2,12,10,1,2,14,14,10)
             , na.strings=c(" ")
             , col.names = c("UF", "ESTRATO_POF", "TIPO_SITUACAO_REG",
                             "COD_UPA", "NUM_DOM", "NUM_UC", "QUADRO",
                             "V9001", "V9002", "V8000", "V9010", "V9011",
                             "DEFLATOR", "V8000_DEFLA", "COD_IMPUT_VALOR",
                             "FATOR_ANUALIZACAO", "PESO", "PESO_FINAL",
                             "RENDA_TOTAL")
             , dec="."
    )
  # save input in a format that is faster to read in R
  saveRDS(aluguel_estimado, "data/raw2clean/pof_ibge/input/ALUGUEL_ESTIMADO.rds")

  } else {

  # read input
  aluguel_estimado <- readRDS("data/raw2clean/pof_ibge/input/ALUGUEL_ESTIMADO.rds")

}

if(!file.exists("data/raw2clean/pof_ibge/input/RENDIMENTO_TRABALHO.rds")) {

  rendimento_trabalho <-
  read.fwf("data/raw2clean/pof_ibge/input/RENDIMENTO_TRABALHO.txt"
           , widths = c(2,4,1,9,2,1,2,2,1,1,7,1,1,1,1,1,1,7,7,
                        7,7,2,2,3,1,12,10,10,10,10,1,1,14,14,
                        10,4,5)
           , na.strings=c(" ")
           , col.names = c("UF", "ESTRATO_POF", "TIPO_SITUACAO_REG",
                           "COD_UPA", "NUM_DOM", "NUM_UC",
                           "COD_INFORMANTE", "QUADRO", "SUB_QUADRO",
                           "SEQ", "V9001", "V5302", "V53021", "V5303",
                           "V5304", "V5305", "V5307", "V8500", "V531112",
                           "V531122", "V531132", "V9010", "V9011",
                           "V5314", "V5315", "DEFLATOR", "V8500_DEFLA",
                           "V531112_DEFLA", "V531122_DEFLA",
                           "V531132_DEFLA", "COD_IMPUT_VALOR",
                           "FATOR_ANUALIZACAO", "PESO", "PESO_FINAL",
                           "RENDA_TOTAL","V53011","V53061")
           , dec="."
  )
  # save input in a format that is faster to read in R
  saveRDS(rendimento_trabalho, "data/raw2clean/pof_ibge/input/RENDIMENTO_TRABALHO.rds")

  } else {

  # save input in a format that is faster to read in R
  rendimento_trabalho <- readRDS("data/raw2clean/pof_ibge/input/RENDIMENTO_TRABALHO.rds")
}

if(!file.exists("data/raw2clean/pof_ibge/input/OUTROS_RENDIMENTOS.rds")) {

  outros_rendimentos <-
    read.fwf("data/raw2clean/pof_ibge/input/OUTROS_RENDIMENTOS.txt"
             , widths = c(2,4,1,9,2,1,2,2,2,7,10,10,2
                          ,2,12,10,10,1,1,14,14,10
             )
             , na.strings=c(" ")
             , col.names = c("UF", "ESTRATO_POF", "TIPO_SITUACAO_REG",
                             "COD_UPA", "NUM_DOM", "NUM_UC",
                             "COD_INFORMANTE", "QUADRO", "SEQ", "V9001",
                             "V8500", "V8501", "V9010", "V9011",
                             "DEFLATOR", "V8500_DEFLA", "V8501_DEFLA",
                             "COD_IMPUT_VALOR", "FATOR_ANUALIZACAO",
                             "PESO", "PESO_FINAL", "RENDA_TOTAL")
             , dec="."
    )
  # save input in a format that is faster to read in R
  saveRDS(outros_rendimentos, "data/raw2clean/pof_ibge/input/OUTROS_RENDIMENTOS.rds")

  } else {

  # read input
  outros_rendimentos <- readRDS("data/raw2clean/pof_ibge/input/OUTROS_RENDIMENTOS.rds")
}


if(!file.exists("data/raw2clean/pof_ibge/input/DOMICILIO.rds")) {

  domicilio <-
    read.fwf("data/raw2clean/pof_ibge/input/DOMICILIO.txt"
             , widths = c(2,4,1,9,2,1,1,1,1,2,1,1,1,1,1,1,1,1,1,2,
                          1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,14,14,1
             )
             , na.strings=c(" ")
             , col.names = c("UF", "ESTRATO_POF", "TIPO_SITUACAO_REG",
                             "COD_UPA", "NUM_DOM", "V0201", "V0202",
                             "V0203", "V0204", "V0205", "V0206", "V0207",
                             "V0208", "V0209", "V02101", "V02102",
                             "V02103", "V02104", "V02105", "V02111",
                             "V02112", "V02113", "V0212", "V0213",
                             "V02141", "V02142", "V0215", "V02161",
                             "V02162", "V02163", "V02164", "V0217",
                             "V0219", "V0220", "V0221", "PESO",
                             "PESO_FINAL", "V6199")
             , dec="."
    )

  # save input in a format that is faster to read in R
  saveRDS(domicilio, "data/raw2clean/pof_ibge/input/DOMICILIO.rds")

  } else {

  # save input in a format that is faster to read in R
    domicilio <- readRDS("data/raw2clean/pof_ibge/input/DOMICILIO.rds")
}





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# SELECT FAMILIES OF INTEREST (RELATED TO LIVESTOCK)
morador <-
  morador %>%
  dplyr::filter(UF %in% c(11:17, 21, 51)) %>% # select states in the Amazon Biome
  dplyr::left_join(rendimento_trabalho[, c("UF", "ESTRATO_POF", "TIPO_SITUACAO_REG", "COD_UPA", "NUM_DOM", "NUM_UC", "COD_INFORMANTE", "V53061")]) %>%
  dplyr::group_by(UF, ESTRATO_POF, TIPO_SITUACAO_REG, COD_UPA, NUM_DOM, NUM_UC) %>%
  dplyr::mutate(d_livestockRelated = if_else(V53061 %in% c(1201:1209, 1402, 1999), 1, 0), # select activities realted to livestock
                d_livestockRelated = sum(d_livestockRelated)) %>% # sum number of informants in a UC with activity related to livestock
  dplyr::filter(d_livestockRelated > 0) %>% # keep only families with at least one informant in a UC with activity related to livestock
  dplyr::slice(1) %>% # keep only one observation per UC
  dplyr::select(UF, ESTRATO_POF, TIPO_SITUACAO_REG, COD_UPA, NUM_DOM, NUM_UC, PESO_FINAL)

# extract id for each UC (family) of interest
aux.id <- morador %>% tidyr::unite(id, UF, COD_UPA, NUM_DOM, NUM_UC) %>% pull(id) %>% unique()

# SELECT CONSUMPTION EXPENDITURES CODES
aux.consumption.code <- tradutor_despesa %>% filter(Descricao_2 %in% c("Despesas de consumo", "Despesas de Consumo")) %>% pull(Codigo)

# CALCULATE TOTAL CONSUMPTION EXPENDITURE BY UC
aluguel_estimado <-
  aluguel_estimado %>%
  dplyr::mutate(expenditure_code = trunc(V9001/100)) %>%  # create expenditure code to identify consumption expenditures
  dplyr::filter(expenditure_code %in% aux.consumption.code) %>% # keep only consumption expenditures
  tidyr::unite(id, UF, COD_UPA, NUM_DOM, NUM_UC, remove = FALSE) %>%
  dplyr::filter(id %in% aux.id) %>%
  dplyr::mutate(aluguel_value =(V8000_DEFLA*V9011*FATOR_ANUALIZACAO*PESO_FINAL)) %>% # calculate annual value
  dplyr::group_by(UF, COD_UPA, NUM_DOM, NUM_UC) %>%
  dplyr::summarise(value = sum(aluguel_value, na.rm = TRUE)) %>%
  dplyr::select(UF, COD_UPA, NUM_DOM, NUM_UC, value)

despesa_coletiva <-
  despesa_coletiva %>%
  dplyr::mutate(expenditure_code = trunc(V9001/100)) %>%  # create expenditure code to identify consumption expenditures
  dplyr::filter(expenditure_code %in% aux.consumption.code) %>% # keep only consumption expenditures
  tidyr::unite(id, UF, COD_UPA, NUM_DOM, NUM_UC, remove = FALSE) %>%
  dplyr::filter(id %in% aux.id) %>%
  dplyr::mutate(despesa_coletiva_value =ifelse( QUADRO==10|QUADRO==19,
                                                (V8000_DEFLA*V9011*FATOR_ANUALIZACAO*PESO_FINAL),
                                                (V8000_DEFLA*FATOR_ANUALIZACAO*PESO_FINAL))) %>% # calculate annual value
  dplyr::group_by(UF, COD_UPA, NUM_DOM, NUM_UC) %>%
  dplyr::summarise(value = sum(despesa_coletiva_value, na.rm = TRUE)) %>%
  dplyr::select(UF, COD_UPA, NUM_DOM, NUM_UC, value)

caderneta_coletiva <-
  caderneta_coletiva %>%
  dplyr::mutate(expenditure_code = trunc(V9001/100)) %>%  # create expenditure code to identify consumption expenditures
  dplyr::filter(expenditure_code %in% aux.consumption.code) %>% # keep only consumption expenditures
  tidyr::unite(id, UF, COD_UPA, NUM_DOM, NUM_UC, remove = FALSE) %>%
  dplyr::filter(id %in% aux.id) %>%
  dplyr::mutate(caderneta_coletiva_value =(V8000_DEFLA*FATOR_ANUALIZACAO*PESO_FINAL)) %>% # calculate annual value
  dplyr::group_by(UF, COD_UPA, NUM_DOM, NUM_UC) %>%
  dplyr::summarise(value = sum(caderneta_coletiva_value, na.rm = TRUE)) %>%
  dplyr::select(UF, COD_UPA, NUM_DOM, NUM_UC, value)

despesa_individual <-
  despesa_individual %>%
  dplyr::mutate(expenditure_code = trunc(V9001/100)) %>%  # create expenditure code to identify consumption expenditures
  dplyr::filter(expenditure_code %in% aux.consumption.code) %>% # keep only consumption expenditures
  tidyr::unite(id, UF, COD_UPA, NUM_DOM, NUM_UC, remove = FALSE) %>%
  dplyr::filter(id %in% aux.id) %>%
  dplyr::mutate(despesa_individual_value = ifelse( QUADRO==44|QUADRO==47|QUADRO==48|QUADRO==49|QUADRO==50,
                                                   (V8000_DEFLA*V9011*FATOR_ANUALIZACAO*PESO_FINAL),
                                                   (V8000_DEFLA*FATOR_ANUALIZACAO*PESO_FINAL))) %>% # calculate annual value
  dplyr::group_by(UF, COD_UPA, NUM_DOM, NUM_UC) %>%
  dplyr::summarise(value = sum(despesa_individual_value, na.rm = TRUE)) %>%
  dplyr::select(UF, COD_UPA, NUM_DOM, NUM_UC, value)

# combine all consumption expenditures
consumption_expenditure <-
  rbind(aluguel_estimado, despesa_coletiva, caderneta_coletiva,despesa_individual) %>%
  dplyr::group_by(UF, COD_UPA, NUM_DOM, NUM_UC) %>%
  dplyr::summarise(value = sum(value, na.rm = TRUE))

# sum all consumption expenditures by UF
clean.pof <-
  consumption_expenditure %>%
  dplyr::select(UF, value) %>%
  dplyr::group_by(UF) %>%
  dplyr::summarise(pof_consumptionExpenditure_livestockRelated = sum(value, na.rm = TRUE)) %>%
  dplyr::select(UF, pof_consumptionExpenditure_livestockRelated)

# add number of livestock related families by UF
clean.pof <-
  morador %>%
  dplyr::group_by(UF) %>%
  dplyr::summarise(pof_numberFamilies_livestockRelated = sum(PESO_FINAL)) %>%
  dplyr::right_join(clean.pof)

# clean environment
rm(aluguel_estimado, caderneta_coletiva, despesa_coletiva, despesa_individual, domicilio, morador, consumption_expenditure, outros_rendimentos, rendimento_trabalho, tradutor_despesa)




# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(clean.pof$pof_numberFamilies_livestockRelated )         <- "number of families with at least one member working in a livestock related job per state (UF)"
sjlabelled::set_label(clean.pof$pof_consumptionExpenditure_livestockRelated)  <- "total annual consumption expediture of families with at least one member working in a livestock related job per state (UF)"



# POST-TREATMENT OVERVIEW
# summary(clean.pof)
# View(clean.pof)



# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.pof,
     file = file.path("data/raw2clean/pof_ibge/output",
                      "clean_pof.Rdata"))


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------