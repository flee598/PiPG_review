# Clean up aotea data - species name inconsistencies (typos, mergers, renames)

aotea_raw <- read.csv(file = "../Data/aotea_allyears.csv", row.names = 1) #vegan uses row names, but the tidyverse doesn"t!

aotea_tidy <- aotea_raw %>%
  mutate(Adihis = Adihis + Adihisp) %>% select(-Adihisp) %>%   # Merge and remove...
  mutate(AlseuoSpp = AlsSpp + Alseuo) %>% select(-c( AlsSpp, Alseuo)) %>%
  mutate(Blenov = Blecha + Blenov) %>% select(-Blecha) %>%  
  mutate(CollospSpp = CollSpp + ColSpp) %>% select(-c(CollSpp, ColSpp)) %>%  
  mutate(Dooaus = Doomed + Dooaus) %>% select(-Doomed) %>%  
  mutate(Genlig = Genlig + Genrup) %>% select(-Genrup) %>%  
  mutate(GrammitisSpp = GramSpp + Graspp) %>% select(-c(GramSpp, Graspp)) %>%  
  mutate(EarinaSpp = EarSpp + EarinaSpp) %>% select(-EarSpp) %>%  
  mutate(Kunrob = Kuneri + Kunrob) %>% select(-Kuneri) %>%  
  mutate(Pteesc = Pteesc + Pttesc) %>% select(-Pttesc) %>%
  mutate(Oplhir = Oplhir + Oplimb) %>% select(-Oplimb) %>%
  mutate(NerteraSpp = Nertera + EarinaSpp) %>% select(-Nertera) %>% 
  mutate(Micopus = Phydiv + Micpus) %>% select(-Phydiv) %>%
  mutate(Pipexc = Macexc + Pipexc) %>% select(-Macexc) %>%
  mutate(Metper = Metper + Metperf) %>% select(-Metperf) %>%
  mutate(MeteroSpp = MetSpp + MetUNID) %>% select(-c(MetSpp, MetUNID)) %>%
  mutate(UnciniaSpp = UnciniaSpp + UncSpp) %>% select(-UncSpp) %>%
  mutate(Muecom = Muecom + Meucom) %>% select(-Meucom) %>%
  mutate(Ripsca = Ripsca + Rhisca + Rhisco) %>% select(-c(Rhisca, Rhisco)) %>%
  mutate(PseudopanaxSpp = PseudopanaxSpp + PseSpp) %>% select(-PseSpp)

aotea_tidy <- aotea_tidy[, -grep("Dead|UNID|Exotic|Cladonia", names(aotea_tidy), ignore.case = TRUE)] # Drop unneeded "species"

table(str_sub(rownames(aotea_tidy), 1,3))

# Remove the extra Glenfern sites and Claris fire
sites <- rownames(aotea_tidy)
aotea_tidy <- aotea_tidy[!(str_starts(sites, "CL") | (str_starts(sites, "GF") & ((str_ends(sites, "13") | str_ends(sites, "14"))))), ]

# aotea_chk <- aotea_tidy[grep(focal_sites, rownames(aotea_tidy)),] # Keep a selection of sites
# aotea_chk <- aotea_chk[-grep("CL_", rownames(aotea_chk)),] # Keep a selection of sites

table(str_sub(rownames(aotea_tidy), 1,3))

# Check no empty species (columns)
aotea_tidy <- aotea_tidy[, colSums(aotea_tidy) > 0]

write.csv(aotea_tidy, "perry_2010_aotea.csv")  
# aotea_pa <- decostand(aotea_chk, "pa")

#xy_mds <- metaMDS(aotea_pa, distance = "bray", autotransform = FALSE, wascores = FALSE)
