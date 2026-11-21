# Export the manuscript panels after running the numbered study scripts.
# Reuses saved plots without repeating model fitting or simulation.
.script_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
.script_file <- if (length(.script_files)) tail(.script_files, 1) else
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
source(file.path(dirname(normalizePath(.script_file)), "Project_Config.R"), encoding = "UTF-8")
project <- PROJECT_ROOT
revision <- project
saved <- file.path(project, "Scripts")
library(ggplot2)
library(dplyr)
library(ggbreak)
out <- file.path(project, "Figures/Manuscript/Panels")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
specs <- list(
  Figure3 = c(A="Study1/figures/Fig1_slope_b1.rds", B="Study1/figures/Fig1_slope_b2.rds"),
  Figure4 = c(A="Study1/figures/Fig1_weights_b1.rds", B="Study1/figures/Fig1_NDOM_b1.rds", C="Study1/figures/Fig1_weights_b2.rds", D="Study1/figures/Fig1_NDOM_b2.rds", E="Study1/figures/Fig1_baseline_weights.rds", F="Study1/figures/Fig1_VDOM.rds"),
  Figure5 = c(A="Study1/figures/Fig2_Slope_Visual.rds", B="Study1/figures/Fig2_Slope_Auditory.rds", C="new:Fig2_ModelFitting.rds", D="Study1/figures/Fig2_confusion_matrix_AIC.rds", E="Study1/figures/Fig2_confusion_matrix_BIC.rds", F="Study1/figures/Fig2_param_recovery.rds", G="new:Fig3_model_validation_cell_means.rds", H="new:Fig3_fit_observed_vs_predicted.rds", I="new:Fig3_weighting_functions_5plq.rds"),
  Figure6 = c(A="Study2/figures/03_ModelComparison.rds", B="Study2/figures/Validation_BestModel.rds", C="Study2/figures/05_Scatter_3p.rds", D="Study2/figures/Fig3_weighting_functions_3p.rds"),
  Figure7 = c(A="Study3/figures/03_ModelComparison.rds", B="Study3/figures/Validation_BestModel.rds", C="Study3/figures/05_Scatter_7p.rds", D="Study3/figures/Fig3_weighting_functions_7p.rds")
)
for (figure in names(specs)) {
  plots <- lapply(specs[[figure]], function(path) {
    filename <- if (startsWith(path,"new:")) file.path(project,"Scripts/Study1/figures",sub("new:","",path)) else file.path(saved,path)
    p <- style_figure(readRDS(filename))
    p <- p + theme(legend.location="plot", legend.key.spacing.x=grid::unit(14,"pt"),
                    legend.box="horizontal", legend.justification="center",
                    legend.margin=margin(0,0,10,0))
    if (is.character(p$labels$x) && nchar(p$labels$x)>24) {
      p <- p + labs(x=sub("Information Criterion", "Information\nCriterion", p$labels$x, fixed=TRUE))
    }
    p <- p + guides(colour=guide_legend(title.position="top", nrow=1))
    if (figure %in% c("Figure5","Figure6","Figure7")) {
      p <- p+theme(axis.title=element_text(size=32),axis.title.x=element_text(size=32),axis.title.y=element_text(size=32))
    }
    if (figure=="Figure5" && any(vapply(p$layers,function(x)inherits(x$geom,"GeomTile"),logical(1)))) {
      value <- if("Correlation" %in% names(p$data)) p$data$Correlation else p$data$prop
      p$data$.matrix_label <- sub("^(-?)0\\.","\\1.",sprintf("%.2f",value))
      fill_scale <- p$scales$get_scales("fill")
      fill_scale$breaks <- if("Correlation" %in% names(p$data)) c(-1,0,1) else c(0,.5,1)
      correlation <- "Correlation" %in% names(p$data)
      p <- p + labs(y=if(correlation) "Recovered\nParameter" else "Recovered\nModel") +
        theme(legend.position="bottom", legend.direction="horizontal",
              legend.location="plot", legend.justification="center",
              legend.title=element_text(size=26, margin=margin(b=8)),
              legend.text=element_text(size=26, margin=margin(t=9)),
              legend.ticks=element_blank(), legend.margin=margin(t=8,b=6),
              legend.box.margin=margin(), legend.box.spacing=grid::unit(12,"pt"),
              axis.title.y=element_text(size=28), plot.margin=margin(8,8,8,8)) +
        guides(fill=guide_colourbar(title.position="top",label.position="bottom",
               display="rectangles",nbin=50,barwidth=grid::unit(145,"pt"),
               barheight=grid::unit(13,"pt")))
      for(j in seq_along(p$layers))if(inherits(p$layers[[j]]$geom,"GeomText")) {
        p$layers[[j]]$aes_params$size<-28*25.4/72.27
        p$layers[[j]]$mapping$label <- aes(label=.matrix_label)$label
      }
    }
    p
  })
  if (figure=="Figure4") {
    for (letter in c("A","C")) plots[[letter]] <- plots[[letter]]+labs(x="Unimodal Valence\nRating")
    grobs <- lapply(plots,ggplotGrob)
    widths <- do.call(grid::unit.pmax,lapply(grobs,function(g)g$widths))
    heights <- do.call(grid::unit.pmax,lapply(grobs,function(g)g$heights))
    for (letter in names(grobs)) {grobs[[letter]]$widths<-widths;grobs[[letter]]$heights<-heights}
    plots <- grobs
  }
  if (figure=="Figure5") {
    # Align the recovery panels without changing their matrix values.
    matrices <- lapply(plots[c("D","E","F")],ggplotGrob)
    widths <- do.call(grid::unit.pmax,lapply(matrices,function(g)g$widths))
    heights <- do.call(grid::unit.pmax,lapply(matrices,function(g)g$heights))
    for (letter in names(matrices)) {
      matrices[[letter]]$widths <- widths
      matrices[[letter]]$heights <- heights
      matrices[[letter]]$respect <- FALSE
      plots[[letter]] <- matrices[[letter]]
    }
  }
  for (letter in names(plots)) {
    export_figure(file.path(out,paste0(figure,"_",letter,".png")),plots[[letter]],
                  width=if(figure=="Figure3")10 else FIGURE_WIDTH_IN,
                  height=if(figure=="Figure3")5 else FIGURE_HEIGHT_IN)
  }
}
# Use a compact count graphic for the four participant-deletion comparisons.
summary_files <- file.path(project,"Results",c("Study2","Study3"),"ParticipantSensitivity/Summary.csv")
combined <- bind_rows(lapply(seq_along(summary_files),function(i) read.csv(summary_files[i]) %>% mutate(Study=paste("Study",i+1))))
p <- plot_selection_frequency(combined)
for (j in seq_along(p$layers)) if(inherits(p$layers[[j]]$geom,"GeomText")) p$layers[[j]]$aes_params$size<-20*25.4/72.27
p <- p+theme(axis.text=element_text(size=22),axis.text.x=element_text(size=22),axis.text.y=element_text(size=22),axis.title=element_text(size=24),axis.title.x=element_text(size=24))
export_figure(file.path(project,"Figures/ParticipantSensitivity/SelectionFrequency.png"),p,11,6.5)
cat("Exported revised figure panels.\n")
