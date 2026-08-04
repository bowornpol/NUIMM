devtools::load_all("D:/NUIMM/NUIMM")

out_dir <- "D:/NUIMM/NUIMM_test/Dung_project/results_KO_maaslin3/Multi-strain_delta/mln_final_test"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

NUIMM:::con_mln_int(
  gsea_file = "D:/NUIMM/NUIMM_test/Dung_project/results_KO_maaslin3/Multi-strain_delta/ppn_output/gsea_results_V1_vs_V3.csv",
  mpn_file = "D:/NUIMM/NUIMM_test/Dung_project/results_KO_maaslin3/Multi-strain_delta/mpn_output/filtered_microbe_pathway_edges.csv",
  ppn_file = "D:/NUIMM/NUIMM_test/Dung_project/results_KO_maaslin3/Multi-strain_delta/ppn_output/pathway_pathway_network_V1_vs_V3.csv",
  pmn_file = "D:/NUIMM/NUIMM_test/Dung_project/results_KO_maaslin3/Multi-strain_delta/pmn_output/filtered_pathway_metabolite_edges.csv",
  output_dir = out_dir,
  visualize = TRUE,
  layout_method = "circle",
  node_colors = c("Microbe"="#b4b992", "Pathway"="#d2a9a4", "Metabolite"="#567285"),
  node_shapes = c("Microbe"="square", "Pathway"="dot", "Metabolite"="diamond"),
  base_node_size = 30,
  plot_width = 12,
  plot_height = 8,
  plot_dpi = 300,
  ppn_map_database = "kegg",
  map_file = NULL
)
print("TEST SCRIPT FINISHED")
