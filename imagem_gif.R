pacman::p_load(tidyverse, gganimate,magick, refund)

data("DTI")

DTI2 <- DTI %>% 
  drop_na() %>% 
  filter(visit == 1 & case == 1)


gerar_grafico <- function(){

# arrumando os dados para o gráfico 
completo <- DTI2$cca %>% 
  t() %>% 
  as_tibble() %>% 
  janitor::clean_names() %>% 
  mutate(
    tract = 1:93
  )  
#
amostra <- DTI2 %>% 
  slice_sample(n=sample(1:10,1), replace=FALSE)
#
amostra <- amostra$cca %>% 
  t() %>% 
  as_tibble() %>% 
  janitor::clean_names() %>% 
  mutate(
    tract = 1:93
  )  %>% 
  rename_with(~ str_replace(., "_1", "_novo"), contains("_1"))
#
p<- left_join(completo,amostra,by="tract") %>% 
  pivot_longer(
    cols = contains("_1"),
    names_to = "original",
    values_to = "valor_original") %>% 
  pivot_longer(
    cols = contains("novo"),
    names_to = "amostra",
    values_to = "valor_amostra"
  ) %>% 
  group_by(tract) %>% 
  ggplot(aes(tract,valor_original, group = original))+ #gráfico de todos os 66 casos completos
  geom_line(color = "gray")+
  geom_line(aes(tract, valor_amostra, group=amostra,color=amostra), linewidth = 1.2)+
  theme(legend.position = "none")+
  labs(x="tract",
       y ="Anisotropia Fracionária (AF)",
       title = "imagem de tensor de difusão:CCA")

return(p)
}


frames <- vector("list", 20)
for (i in 1:length(frames)) {
  frames[[i]] <- gerar_grafico()
}

# Obter diretório de trabalho atual
dir_atual <- getwd()

# Criar diretório temporário
dir_temp <- tempdir()

# Mudar para o diretório temporário
setwd(dir_temp)

# Salvar gráficos como imagens temporárias
tmp_files <- lapply(seq_along(frames), function(i) {
  tmp_file <- file.path(dir_temp, paste0("frame", i, ".png"))
  ggsave(tmp_file, frames[[i]], dpi = "retina", width = 6, height = 4)
  return(tmp_file)
})

# Corrigir caminhos dos arquivos no Windows
if (.Platform$OS.type == "windows") {
  tmp_files <- file.path(gsub("/", "\\\\", tmp_files))
}

# Criar a animação com as imagens temporárias
animacao <- image_read(tmp_files)
animacao <- image_animate(animacao, fps = 1)

# Salvar a animação como GIF
anim_save(file.path(dir_atual, "animacao.gif"), animacao)

# Voltar para o diretório de trabalho original
setwd(dir_atual)
