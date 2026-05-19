library(NUIMM)
con_mln(
  gene_abundance_file = system.file("extdata", "filtered_reaction_abundance_yashida_filtered_cpm.tsv", package="NUIMM"),
  pathway_abundance_file = system.file("extdata", "filtered_pathway_abundance_yashida_filtered_cpm_cleaned.csv", package="NUIMM"),
  pathway_contribution_file = system.file("extdata", "filtered_pathway_contribution_yashida_filtered_cpm_long_format.tsv", package="NUIMM"),
  metabolite_concentration_file = system.file("extdata", "met_yashida_processed.csv", package="NUIMM"),
  metadata_file = system.file("extdata", "metadata_yashida_cleaned.csv", package="NUIMM"),
  output_directory = "D:/NUIMM/NUIMM_web/test_out",
  ppn_da_method = "simple"
)
iden_hub(
  multi_layered_network_file = "D:/NUIMM/NUIMM_web/test_out/multi_layered_network.csv",
  output_directory = "D:/NUIMM/NUIMM_web/test_out/hub"
)
find_path(
  multi_layered_network_file = "D:/NUIMM/NUIMM_web/test_out/multi_layered_network.csv",
  output_directory = "D:/NUIMM/NUIMM_web/test_out/path"
)
node_prior(
  multi_layered_network_file = "D:/NUIMM/NUIMM_web/test_out/multi_layered_network.csv",
  output_directory = "D:/NUIMM/NUIMM_web/test_out/prior"
)
