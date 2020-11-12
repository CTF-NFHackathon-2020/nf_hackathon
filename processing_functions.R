"
Quick function to find orthologous gene names between mouse and human genomes for a given gene set.
gene_name <- Gene symbols
from_species <- Species from which input gene names are derived.
to_species <- Species from which orthologous gene names are to be returned.

Returns:
A data.table with two columns --- 'mouse_gene_name' and 'human_gene_name'. 
Note that there may be one-many mappings possible.
"
get_orthologs <- function(gene_names,from_species,to_species) {
    if (from_species == "mouse") {
        from_mart <- mart_mouse
        to_mart <- mart_human
    } else {
        from_mart <- mart_human
        to_mart <- mart_mouse
    }

    ortholog_dt <- getLDS(attributes=c("external_gene_name"),
            filters="external_gene_name", values=gene_names, mart=from_mart,
            attributesL=c("external_gene_name"), martL=to_mart, verbose=F) %>% data.table(.) %>%
     setnames(.,"Gene.name.1",paste(to_species,"gene_name",sep="_")) %>%
    setnames(.,"Gene.name",paste(from_species,"gene_name",sep="_"))

    return(ortholog_dt)
}

"
The function below loads the E9.5 trunk neural crest data.
In order for this function to run, two sets of files need to be downloaded to the data_dir path ---
1. The .txt.gz read count matrices from GSE129114.
2. 'Table S9.csv' from http://pklab.med.harvard.edu/ruslan/neural_crest/Table%20S9.csv .

This function returns a Seurat object of E9.5 trunk neural crest cells. 
"
load_trunk_neural_crest <- function() {
    data_files <- list.files(file.path(data_dir),full.names = T)

    read_count_dt <- data.table()
    
    "
    The loop below iteratively reads each of the gene expression matrices into a combined gene expression matrix.
    fread actually creates a data.table, which is fast to read from disk.
    "
    for (data_file in data_files) {
        if (grepl("GSE129114",data_file) && grepl(".txt.gz",data_file)) {
            dt <- fread(data_file)

            if (nrow(read_count_dt) == 0) {
                read_count_dt <- dt
            } else {
                read_count_dt <- merge(read_count_dt,dt,by="V1")
            }
        }
    }
    "
    Create a matrix from the read count data.table and create an (unannotated) Seurat object.
    The 'V1' column is actually the gene names. The name V1 is assigned automatically by data.table during the fread()
    call.
    "
    read_count_mat <- as.matrix(read_count_dt[,!c("V1")])
    rownames(read_count_mat) <- read_count_dt$V1
    rm(read_count_dt)
    mouse_ncc_obj <- CreateSeuratObject( read_count_mat )
    mouse_ncc_meta_data_dt <- data.table( mouse_ncc_obj@meta.data, keep.rownames=T ) %>% setnames("rn","cell.name")

    "
    Remember to download Table\ S9.csv from the link specified in point 2 above the function header.
    This file contains cluster assignments for every trunk neural crest cell.
    "
    ncc_cluster_info_dt <- fread(file.path(data_dir,"Table\ S9.csv"),skip = 1)


    "
    Boiler plate code to add annotations from Table S9 into Seurat object.
    "
    meta_data_df <- ncc_cluster_info_dt %>% dplyr::select(-V3) %>% tibble::column_to_rownames("V1") %>%dplyr::rename(cluster=V4)
    mouse_ncc_obj <- AddMetaData( mouse_ncc_obj, meta_data_df ) %>% SetIdent(.,value="V2")
    trunk_ncc_obj <- subset( mouse_ncc_obj, ident="E9.5_trunk_Wnt1") %>% NormalizeData %>% SetIdent(value="cluster")

    "
    The batch_info column added to the meta.data helps correct batch effects during trajectory construction. A
    'cell.type' column is also added since the other Seurat objects contain a cell.type column.
    There is also a single un-annotated cell here (with the cluster label NA). We get rid of that cell.
    "
    trunk_ncc_obj <- AddMetaData( trunk_ncc_obj, 
                                tibble::rownames_to_column(trunk_ncc_obj@meta.data,"cell.name") %>%
                                mutate(batch_info=as.factor(gsub("_[A-Z][0-9]*","",cell.name))) %>%
                                mutate(cell.type=cluster) %>% 
                                 tibble::column_to_rownames("cell.name") ) %>%
                                 subset(subset=cell.type %in% c("mesenchyme","neural tube","autonomic","early migratory","migratory 1","migratory 2","sensory","delaminatory"))

    return(trunk_ncc_obj)
}

"
The function below loads Schwann cell precursors from E12.5 and E13.5 mouse embryos. Data sourced from
https://science.sciencemag.org/content/357/6346/eaal3753.abstract .

This requires two Seurat objects --- E12.5_seurat_annotated.rds and E13.5_seurat_annotated.rds . These objects will be
made avaiable for download shortly.
"
load_sc_precursors <- function() {
    E12.5_seurat_obj <- readRDS(file.path(data_dir,"E12.5_seurat_annotated.rds"))
    E13.5_seurat_obj <- readRDS(file.path(data_dir,"E13.5_seurat_annotated.rds"))

    merged_df <- merge( data.table( E12.5_seurat_obj@meta.data, keep.rownames=T ), 
          tibble::enframe(Idents(E12.5_seurat_obj), name="rn", value="cluster") %>%
          mutate(cell.type=paste("E12.5",cluster,sep="_"),batch_info="E12.5"), by="rn" ) %>%
    tibble::column_to_rownames("rn")
    E12.5_seurat_obj <- AddMetaData( E12.5_seurat_obj, merged_df )
    merged_df <- merge( data.table( E13.5_seurat_obj@meta.data, keep.rownames=T ), 
          tibble::enframe(Idents(E13.5_seurat_obj), name="rn", value="cluster") %>% 
          mutate(cell.type=paste("E13.5",cluster,sep="_"),batch_info="E13.5"), by="rn" ) %>%
    tibble::column_to_rownames("rn")
    E13.5_seurat_obj <- AddMetaData( E13.5_seurat_obj, merged_df )

    "
    The two time-points also contain chromaffin, bridge and sympathoblast cells. We aren't particularly interested 
    in these populations for our current analysis. So we retain only SCPs (Schwann Cell Precursors) for trajectory
    construction.
    "
    merged_obj <- merge( E12.5_seurat_obj, E13.5_seurat_obj) %>% subset(subset=cell.type %in% c("E12.5_SCP","E13.5_SCP"))



    return(merged_obj)
}

"
The function below loads mature schwann cell data. The GEO accession data made available by the authors
unfortunately contained errors. The authors provided us with a Seurat object containing correct read count and
annotation information. This Seurat object was used for our analyses, and is the only input to this function.
"
load_sc_data <- function() { 
    pnas_obj <- readRDS(file.path(data_dir,"PNS_pnas.rds"))

    "
    We add a batch_info column to the meta.data slot so that batch effects can be corrected during Seurat integration.
    "
    pnas_obj <- AddMetaData( pnas_obj, tibble::rownames_to_column(pnas_obj@meta.data,"cell.name") %>%
                             mutate(batch_info="10x") %>% tibble::column_to_rownames("cell.name")) %>%
                           subset(subset=cell.type %in% c("mySC","nmSC"))

    return(pnas_obj)
}

"
Load mouse schwann cell data from plexiform neurofibroma. 

The input is a Seurat object contained Schwann cells and other cells from a PNF-bearing mouse. This object
will be made available for download shortly.
"
load_pnf_scs <- function() {
    mouse_pnf_obj <- readRDS(file.path(data_dir,"mouse_pnf_seurat_object.rds"))

    #Clusters 1 and 3 are Schwann cells. 
    subset_mouse_pnf_obj <- subset( mouse_pnf_obj, subset=RNA_snn_res.0.1 %in% c(1,3))

    return(subset_mouse_pnf_obj)
}


"
A function to deal with the fact that gene names from public sources can contain the dreaded Excel error of converting
some gene names to calendar dates. This brute-force function fixes that.

Input :
1. dt - Data.table object where the gene names are stored in the 'Symbol' column.

Output :
The input data.table with a column 'corrected_gene_name' appended. 

"
fix_gene_names <- function( dt ) {
    rename_dt <- data.table()
    rename_dt$incorrect <- c("1-Sep","2-Sep","3-Sep","4-Sep","5-Sep","6-Sep","7-Sep","8-Sep","9-Sep","10-Sep","11-Sep","12-Sep","13-Sep","14-Sep","15-Sep","2-Oct","2-Oct","3-Oct","4-Oct","6-Oct","6-Oct","7-Oct","9-Oct","11-Oct","1-Dec","1-Mar","2-Mar","3-Mar","4-Mar","5-Mar","6-Mar","7-Mar","8-Mar","9-Mar","10-Mar","11-Mar","1-Feb","2-Feb","5-Feb","6-Feb","7-Feb","9-Feb","10-Feb","2-Apr","3-Apr","1-Apr","1-May","2-Nov","1-Nov")
    rename_dt$correct <- c("SEPT1","SEPT2","SEPT3","SEPT4","SEPT5","SEPT6","SEPT7","SEPT8","SEPT9","SEPT10","SEPT11","SEPT12","SEPT13","SEPT14","SEPT15","POU2F2","SLC22A2","POU5F1","POU5F1","POU3F1","SLC22A16","POU3F2","POU3F4","POU2F3","DEC1","MARCH1","MARCH2","MARCH3","MARCH4","MARCH5","MARCH6","MARCH7","MARCH8","MARCH9","MARCH10","MARCH11","FEB1","FEB2","FEB5","FEB6","FEB7","FEB9","FEB10","FAM215A","ATRAID","MAGEH1","PRKCD","CTGF","C11orf40")
    all_gene_names <- dt$Symbol
    for (idx in 1:nrow(rename_dt)) {
          all_gene_names <- gsub(paste0("^",rename_dt[idx,]$incorrect,"$"),rename_dt[idx,]$correct,all_gene_names)
    }   
    dt$corrected_gene_name <- all_gene_names

    return(dt)
}


"
Data of 16 MPNST samples from https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000792.v1.p1 .

It is currently unclear if the FPKM data from this can be shared freely. We will sort this out and provide this for
download.
"
load_human_mpnst_data <- function() {
    mpnst_gene_exp_dt <- fread(file.path(data_dir,"DiagnosisMergedExpression_medianFPKMFilt_Log2.txt")) %>%
    setnames("V1","Symbol") %>% dplyr::select(c(Symbol,contains("MPNST")))

    scaled_dt <- mpnst_gene_exp_dt[,!c("Symbol")] %>% as.matrix %>% scale %>% as.data.frame %>%
    mutate(Symbol=mpnst_gene_exp_dt$Symbol) %>% data.table
    temp_dt <- fix_gene_names(scaled_dt)[,!c("Symbol")] %>% setnames("corrected_gene_name","Symbol")

    #There's a few gene names that appear multiple times, where each duplicate entry represents 
    #a different transcript.
    merged_mpnst_dt <- data.table(temp_dt)[,lapply(.SD,mean),by=Symbol]


}
"
This function outputs a Monocle3 CDS object with pseudotime values and trajectory topology stored in it.
Inputs :
1. seurat_obj -> A Seurat object.
2. genes_to_use -> An optional argument. A vector of genes to use for assigning pseudotime values to each cell. If set
to NULL, highly variable genes are used by default.


Returns :
A Monocle3 CDS object.
"
compute_pseudotime <- function( seurat_obj, genes_to_use=NULL ) {
    
   meta_data_df <- seurat_obj@meta.data

   #Log-Normalize data and find highly variable genes.
   seurat_obj <- NormalizeData( seurat_obj ) %>% FindVariableFeatures
   expr_mat <- seurat_obj[["RNA"]]@data

   all_genes <- rownames(expr_mat)
   gene_meta_data <- data.frame( gene_short_name=all_genes, row.names=all_genes )

   #Converts a Seurat object into a CDS object.
   cds_obj <- new_cell_data_set( expr_mat, cell_metadata=meta_data_df,
                             gene_metadata = gene_meta_data)
    
    if (is.null(genes_to_use)) {
        genes_to_use <- VariableFeatures(seurat_obj,nfeatures=1000)
    } 

    #norm_method is set to "none" since counts in cds_obj are already log-normalized.
   cds_obj <- preprocess_cds(cds_obj, num_dim = 10, use_genes=genes_to_use,
                            scaling=F, norm_method="none", method="PCA",
                              residual_model_formula_str = "G2M.Score + S.Score")

    cds_obj <- align_cds(cds_obj, alignment_group = "batch_info")
   cds_obj <- reduce_dimension(cds_obj,preprocess_method="PCA")
   cds_obj <- cluster_cells(cds_obj, resolution=1e-4 )

   #The following line is a hack to ensure that all cells get assigned a pseudotime value
   #by ensuring that all clusters belong to the same partition.
    cds_obj@clusters$UMAP$partitions[cds_obj@clusters$UMAP$partitions != "1"] <- "1"
   cds_obj <- learn_graph( cds_obj, use_partition=F )

    #Neural tube cells are given a pseudotime of zero. The code is boiler-plate code from the 
    #Monocle3 manual for non-interactively picking the root cell for the trajectory.
    cell_ids <- WhichCells(seurat_obj,expression = cell.type == "neural tube")
   closest_vertex <- cds_obj@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
   closest_vertex <- as.matrix(closest_vertex[colnames(cds_obj), ])
   root_pr_nodes <-  igraph::V(principal_graph(cds_obj)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]

   cds_obj <- order_cells( cds_obj, root_pr_nodes=c(root_pr_nodes))

    return(cds_obj)
}


"
A function for down-sampling reads from a sub-population of cells to match the UMI/read count distribution
of a target cell population. 

Inputs : 
1. count_mat -- Gene x cell read count matrix.
2. meta_data_dt -- A data.table containing a 'cell.name' column. Typically, this is derived from the @meta.data slot of
the Seurat object with the row names added as a 'cell.name' column.
3. feature_counts_to_subsample -- A vector of read counts to be sub-sampled.
4. feature_counts_reference -- A vector of read counts that serves as the reference read count distribution to which
feature_counts_to_subsample has to be matched.

"
create_resampled_seurat_obj <- function( count_mat, meta_data_dt, feature_counts_to_subsample,
                                        feature_counts_reference ) {

    cell_names_reference <- names(feature_counts_reference)
    num_cells_to_subsample <- length(feature_counts_to_subsample)

    #We make sure that the largest read counts are sub-sampled first.
     target_feature_counts <- sort( sample(feature_counts_reference,size=num_cells_to_subsample,replace=T), decreasing=TRUE )
    feature_counts_to_subsample <- sort(feature_counts_to_subsample,decreasing=TRUE)

    resampled_mat <- count_mat[,cell_names_reference]

    cell_num <- 1
    all_genes <- rownames(count_mat)

    cell_names_to_subsample <- names(feature_counts_to_subsample)
    num_skipped <- 0
    for (cell_name in cell_names_to_subsample) {
        subsampled_feature_count <- target_feature_counts[cell_num]
        if (subsampled_feature_count > feature_counts_to_subsample[cell_num]) {
            #If the cell picked for sub-sampling happens to have fewer read counts than the "target" cell, 
            #then we don't need to sub-sample this cell. We just retain the original read counts of the cell.
            sparse_cell <- as(count_mat[,cell_name],"sparseMatrix" )
            rownames(sparse_cell) <- all_genes
            colnames(sparse_cell) <- cell_name
          resampled_mat <- RowMergeSparseMatrices( resampled_mat, sparse_cell )
            next
        }
        resampled_cell <- SampleUMI( as.matrix(count_mat[,cell_name]), subsampled_feature_count )
        colnames(resampled_cell) <- cell_name
        rownames(resampled_cell) <- all_genes
         resampled_mat <- RowMergeSparseMatrices( resampled_mat, resampled_cell )
        cell_num <- cell_num + 1
    }

    #We now add back any cells that were not re-sampled. 
    remaining_cells <- meta_data_dt[!cell.name %in% c(cell_names_to_subsample,cell_names_reference),cell.name]
    remaining_mat <- count_mat[,remaining_cells]
    resampled_mat <- RowMergeSparseMatrices( remaining_mat, resampled_mat )

    #This ensures that the resampled versions of cells are given the same name as their parent cells.
    cell_name_order <- c(remaining_cells,cell_names_reference,cell_names_to_subsample)
    new_meta_data_df <- as.data.frame( meta_data_dt[match(cell.name,cell_name_order),] )
    rownames(new_meta_data_df) <- new_meta_data_df$cell.name

    #Note that after this merging process, the nFeature_RNA and nCount_RNA values may not make sense. Don't be alarmed.
    resampled_seurat_obj <- CreateSeuratObject( resampled_mat, meta.data=new_meta_data_df )

    return(resampled_seurat_obj)
}


"
This function takes as input the data.frame output by Seurat's FindAllMarkers function, and
creates a list, where each member is a vector of marker genes (with avg_logFC > 0 and p_val_adj < 0.1)

The trouble with FindAllMarkers is that multiple clusters can share an up-regulated marker gene. For our purposes here,
this is problematic. When an up-regulated marker gene is assigned to multiple clusters, we assign the marker gene
to the cluster where it has the highest log-fold change. 

Inputs : 
1. markers_df -- A data.frame output by the FindAllMarkers function (i.e., it needs to contain at least
p_val_adj,gene,cluster and avg_logFC columns)
2. p_val_thresh -- Threshold on the adjusted p-value. Used to remove markers with p_adj > p_val_thresh.
3. num_markers -- Maximum number of up-regulated markers to retain. Set to NULL if all markers are to be used.

Outputs : 
A list, where each member is a vector of marker genes.
"
create_unique_marker_lists <- function(markers_df,p_val_thresh=0.1,num_markers=100) {
    markers_dt <- data.table(markers_df)[avg_logFC > 0 & p_val_adj < p_val_thresh,]  #Data.table is much more convenient to work with than data.frame some times

    #Separate out unique markers and non-unique markers (i.e., markers which are common between multiple clusters)
    non_unique_markers <- markers_dt[,.N,by=gene][N > 1,gene]
    non_unique_markers_df <- markers_dt %>% group_by(gene) %>% top_n(n=1,wt=avg_logFC) %>% ungroup
    unique_markers_df <- markers_dt[!gene %in% non_unique_markers,]

    marker_gene_list <- list()
    for (cluster_ in unique(markers_dt$cluster)) {
        subset_df <- rbind(unique_markers_df %>% dplyr::filter(cluster == cluster_),
                               non_unique_markers_df %>% dplyr::filter(cluster == cluster_)) %>% arrange(-avg_logFC)

        if (!is.null(num_markers)) {
            subset_df <- head(subset_df,num_markers)
        }


        marker_gene_list[[cluster_]] <- subset_df$gene
    }

    return(marker_gene_list)
}


#Cell cycle genes in the human genome
s.genes <- c("MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8")
g2m.genes <- c("HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","HN1","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA")


#BioMart instances for human and mouse genomes
mart_mouse = useMart("ensembl", dataset="mmusculus_gene_ensembl",verbose=F,host="ensembl.org")
mart_human = useMart("ensembl", dataset="hsapiens_gene_ensembl",verbose=F,host="ensembl.org")
#
##Figure out mouse orthologs of S and G2M genes.
temp_dt <- get_orthologs( s.genes,"human","mouse")
mouse_s_genes <- temp_dt$mouse_gene_name
temp_dt <- get_orthologs( g2m.genes,"human","mouse")
mouse_g2m_genes <- temp_dt$mouse_gene_name

