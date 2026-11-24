#!/usr/bin/env Rscript
# Prepare secondary-calibrated screening trees from an IQ-TREE phylogram.
suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
    stop("Usage: 28_make_cafe5_timetrees.R <IQ-TREE phylogram> <output_directory>")
}
input_tree <- normalizePath(args[1], mustWork = TRUE)
output_directory <- args[2]
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

cyclophyllidea <- c(
    "E_granulosus", "E_multilocularis", "H_microstoma",
    "T_multiceps", "T_solium"
)
keep_species <- c(
    "S_solidus", "S_cotti", "Spirometra_Aus1", "E_multilocularis",
    "H_microstoma", "E_granulosus", "T_solium", "T_multiceps"
)

phy <- read.tree(input_tree)
if (inherits(phy, "multiPhylo") || anyDuplicated(phy$tip.label)) {
    stop("Expected one tree with unique tip labels")
}
missing_species <- setdiff(union(cyclophyllidea, keep_species), phy$tip.label)
if (length(missing_species)) {
    stop("Species absent from input tree: ", paste(missing_species, collapse = ", "))
}
check_lengths <- function(tree, stage) {
    lengths <- tree$edge.length
    if (is.null(lengths) || length(lengths) != nrow(tree$edge) ||
        any(!is.finite(lengths)) || any(lengths <= 0)) {
        stop(stage, ": missing, non-finite, or non-positive branch lengths")
    }
}
check_lengths(phy, "Original phylogram")
original_distances <- cophenetic.phylo(phy)
phy$node.label <- NULL
phy$root.edge <- NULL

# ape::root(resolve.root=TRUE) adds a zero-length root edge. Split only
# the selected basal edge into two positive halves before calibration.
# This places a starting root within the chosen edge; it is not dating.
phy <- root(phy, outgroup = cyclophyllidea, resolve.root = TRUE)
root_node <- setdiff(unique(phy$edge[, 1]), phy$edge[, 2])
basal_edges <- which(phy$edge[, 1] == root_node)
if (length(root_node) != 1L || length(basal_edges) != 2L) {
    stop("Expected two branches at the Cyclophyllidea/Diphyllobothriidea root")
}
zero_edges <- which(phy$edge.length == 0)
if (length(zero_edges)) {
    if (length(zero_edges) != 1L || !(zero_edges %in% basal_edges)) {
        stop("Unexpected zero-length branch away from the newly inserted root")
    }
    basal_total <- sum(phy$edge.length[basal_edges])
    if (!is.finite(basal_total) || basal_total <= 0) {
        stop("The selected root edge has no positive length")
    }
    phy$edge.length[basal_edges] <- basal_total / 2
    message("Resolved artificial zero root branch; basal lengths = ",
            paste(format(phy$edge.length[basal_edges], digits = 12), collapse = ", "))
}
check_lengths(phy, "Rerooted phylogram")

# Verify that the root directly separates the two specified taxon sets.
descendant_tips <- function(node, tree) {
    if (node <= Ntip(tree)) return(tree$tip.label[node])
    unlist(lapply(tree$edge[tree$edge[, 1] == node, 2],
                  descendant_tips, tree = tree), use.names = FALSE)
}
root_children <- phy$edge[basal_edges, 2]
if (!any(vapply(root_children, function(node) {
    setequal(descendant_tips(node, phy), cyclophyllidea)
}, logical(1)))) {
    stop("The root does not separate Cyclophyllidea from the remaining taxa")
}

phy <- keep.tip(phy, keep_species)
phy$node.label <- NULL
phy$root.edge <- NULL
if (!is.rooted(phy) || !is.binary(phy)) {
    stop("The pruned tree must be rooted and fully bifurcating")
}
check_lengths(phy, "Eight-species phylogram")
pruned_distances <- cophenetic.phylo(phy)
distance_error <- max(abs(pruned_distances - original_distances[
    rownames(pruned_distances), colnames(pruned_distances)
]))
if (distance_error > 1e-10 * max(1, original_distances)) {
    stop("Rerooting/pruning changed the retained pairwise patristic distances")
}
message("PASS: rerooting and pruning preserve pairwise distances")
write.tree(phy, file.path(output_directory,
                         "screening8_rooted_phylogram_not_dated.nwk"), digits = 15)

# Secondary calibration: T. multiceps genome study, 13-20 Ma estimate
# for the T. multiceps / T. solium split (not a fossil calibration).
# The midpoint below is an analyst-selected screening value.
# Uniform changes to this single fixed age chiefly rescale all times;
# these outputs do not test independent calibration/model uncertainty.
calibration_ages <- c("13Ma" = 13, "16.5Ma" = 16.5, "20Ma" = 20)
summary_rows <- list()
warning_rows <- list()
for (label in names(calibration_ages)) {
    age <- calibration_ages[[label]]
    calibration <- data.frame(
        node = getMRCA(phy, c("T_multiceps", "T_solium")),
        age.min = age, age.max = age, soft.bounds = FALSE
    )
    fit_warnings <- character()
    dated <- withCallingHandlers(
        chronos(phy, lambda = 1, model = "correlated", quiet = TRUE,
                calibration = calibration,
                control = chronos.control(dual.iter.max = 1000,
                                         iter.max = 20000, eval.max = 20000)),
        warning = function(w) {
            fit_warnings <<- c(fit_warnings, conditionMessage(w))
        }
    )
    if (any(grepl("converg|iterations reached|false convergence", fit_warnings,
                  ignore.case = TRUE))) {
        stop("Time-tree optimization needs review for ", label, ": ",
             paste(fit_warnings, collapse = "; "))
    }
    check_lengths(dated, paste("Dated tree", label))
    root_to_tip <- node.depth.edgelength(dated)[seq_len(Ntip(dated))]
    spread <- diff(range(root_to_tip))
    if (spread > 1e-7 * max(1, root_to_tip)) {
        stop("Fitted tree is not ultrametric: ", label)
    }
    dated$node.label <- NULL
    output_file <- file.path(output_directory,
                             paste0("screening8_dated_", label, ".nwk"))
    write.tree(dated, output_file, digits = 15)
    summary_rows[[label]] <- data.frame(
        tree = basename(output_file), calibration_pair = "T_multiceps--T_solium",
        calibration_age_Ma = age, species = Ntip(dated),
        minimum_root_to_tip = min(root_to_tip),
        maximum_root_to_tip = max(root_to_tip), difference = spread,
        fit_warnings = length(fit_warnings)
    )
    warning_rows[[label]] <- data.frame(tree = rep(basename(output_file), length(fit_warnings)),
                                       warning = fit_warnings)
    message("Wrote: ", output_file)
}
write.table(do.call(rbind, summary_rows), file.path(output_directory,
            "timetree_calibration_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, warning_rows), file.path(output_directory,
            "timetree_fit_warnings.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
    paste("Source tree:", input_tree),
    paste("Source tree MD5:", unname(tools::md5sum(input_tree))),
    paste("ape version:", packageVersion("ape")),
    "Method: ape::chronos, correlated model, lambda=1; fixed secondary calibration.",
    "Calibration source: https://pmc.ncbi.nlm.nih.gov/articles/PMC6191302/",
    "13-20 Ma is a published divergence estimate; 16.5 Ma is its arithmetic midpoint.",
    "These are screening time trees; review optimizer warnings before CAFE5.",
    "The three single-age variants mainly test overall time scaling, not independent dating uncertainty.",
    "Source IQ-TREE belongs to the without32 sensitivity run; record that provenance alongside baseline counts."
), file.path(output_directory, "timetree_provenance.txt"))
message("PASS: generated three eight-species screening time trees")
