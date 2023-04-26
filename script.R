# load packages ----

library(googlesheets4)
library(dplyr)
library(janitor)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggh4x)
library(forcats)

# load data and pre-processing ----

samples_legend <- googlesheets4::read_sheet(
  ss = '1xAwGZ2Ppx8sxzECMZTZw5icM5DgEsOHla2ViIKdYqs4',
  sheet = 'mrm_samples'
)

peptides_legend <- googlesheets4::read_sheet(
  ss = '1-sSaPLjMiuEiuVBQ8pif7M95two8nV0hzOoxPvTQTDg',
  sheet = 'peptides'
)

important_columns <- c(
  'peptide_sequence', 'replicate_name', 'ratio_to_standard', 'pass_qc'
)

raw_data <- googlesheets4::read_sheet(
  ss = '1Al1EKZh8haPgoobc-woyb60rrJGi6jJcH_iAA2foZbI',
  sheet = 'raw'
) %>% 
  janitor::clean_names() %>% 
  dplyr::group_by(peptide_sequence) %>% 
  tidyr::fill(pass_qc) %>% 
  dplyr::mutate(
    pass_qc = replace(pass_qc, pass_qc > 0, TRUE),
    pass_qc = replace(pass_qc, is.na(pass_qc), FALSE),
    replicate_name = str_extract(replicate_name, '(?<=_).+(?=_)')
  ) %>% 
  dplyr::select(all_of(important_columns)) %>% 
  merge(
    samples_legend, 
    by = 'replicate_name'
  ) %>% 
  merge(
    peptides_legend,
    by = 'peptide_sequence'
  ) %>% 
  tidyr::uncount(n_protein_groups, .id = "protein_group_id") %>% 
  rowwise() %>% dplyr::mutate(
    protein_name = str_split_i(protein_name, ' & ', protein_group_id),
    aa_subs = str_split_i(aa_subs, ' & ', protein_group_id),
    isoform = replace(isoform, is.na(isoform), '')
  )
dplyr::select(-c('replicate_name'))

# data processing and quality control ----

data_sum <- raw_data %>% 
  rowwise() %>% dplyr::mutate(
    conc = case_when(
      pass_qc == 1 ~ ratio_to_standard * 50,
      pass_qc == 0 ~ ratio_to_standard * 50
    ),
    full_protein = paste(protein_name, aa_subs, isoform, sep = '\n')
  ) %>% 
  dplyr::group_by_at(
    vars(-c('ratio_to_standard', 'conc'))
  ) %>% 
  dplyr::summarise(
    mean_conc = mean(conc, na.rm = TRUE),
    sd_conc = sd(conc),
    cv_conc = sd_conc / mean_conc * 100
  )


data_sum %>% 
  dplyr::mutate(
    replicate_name = fct_reorder(replicate_name, mean_conc)
  ) %>% 
  ggplot2::ggplot(
    aes(weave_factors(as.factor(replicate_name), disease), mean_conc)
  ) +
  geom_point(
    data = dplyr::filter(data_sum, pass_qc == 1, peptide_type != 'common',
                         missed_cleave != TRUE),
    aes(fill = peptide_type),
    position = position_dodge2(0.4), shape = 21,
    size = 4, alpha = 0.7
  ) +
  geom_point(
    data = dplyr::filter(data_sum, pass_qc == 0, peptide_type != 'common',
                         missed_cleave != TRUE),
    aes(color = peptide_type),
    position = position_dodge2(0.4), shape = 13,
    size = 4, stroke = 0.5
  ) +
  facet_wrap(facets = vars(full_protein), scales = 'free', ncol = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  guides(x = "axis_nested")

data_rec <- data_sum %>% 
  dplyr::filter(
    missed_cleave != TRUE,
    peptide_type != 'common',
    pass_qc != 0
  ) %>% 
  dplyr::group_by_at(
    vars(-c('peptide_sequence', 'pass_qc', 'peptide_type', 'mean_conc', 'sd_conc', 'cv_conc'))
  ) %>% 
  dplyr::summarise(
    rec = mean_conc[peptide_type == 'recoded']/(mean_conc[peptide_type == 'recoded'] + mean_conc[peptide_type == 'genomic'])
  )

data_rec %>% 
  ggplot2::ggplot(
    aes(patient, rec, fill = cell_type)
  ) +
  geom_point(position = position_dodge2(0.4), shape = 21,
             size = 4, alpha = 0.7) +
  facet_wrap(facets = vars(full_protein), scales = 'free', ncol = 3)

data_sum_ggplot <- data_sum %>% 
  dplyr::mutate(
    patient = factor(patient, c('FD3S', 'FD5S', 'FF1S', 'FF2S', 'FF4S', 'PO4S', 'BO4S', 'HD76', 'SCA17', 'SCA17M', 'SCA17S')),
    pass_qc = factor(pass_qc, c(0, 1))
  ) %>% 
  dplyr::filter(
    !protein_name %in% c()
  )

data_sum_ggplot %>%
  ggplot(aes(weave_factors(patient, disease), mean_conc)) +
  facet_wrap(vars(protein_name, cell_type), scales = 'free_x', ncol = 2) +
  guides(x = "axis_nested") +
  geom_point(
    aes(shape = pass_qc, color = peptide_type),
    position = position_dodge2(0.1),
    size = 4,
    alpha = 0.5
  )
