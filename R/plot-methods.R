#############################################
# orthologの結果を表示させる
#############################################
# 系統樹と対応して、orthologの有無の表示
# いつものやつ。aplotで作成する
# 4番目にはplasmidと判定された
#############################################
#' Plot Ortyholog and phylogenetic tree
#'
#' Orthofinderとphylogenetic treeをまとめて書く
#'
#' @param ortho
#' @param tree
#' @param metadata
#' @param plas_ortho plasmidとされるorthologのリスト (3列)
#' @param label
#' @param num_ortho
#' @param delete_genomes
#'
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom aplot insert_top
#' @importFrom aplot insert_left
#' @importFrom ggplot2 ggplot
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr arrange
#' @importFrom dplyr summarise
#' @importFrom dplyr group_by
#' @importFrom dplyr rowid_to_column
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom dplyr unchop
#' @importFrom dplyr filter
#' @importFrom tidyr pivot_longer
#'
#'
plot_ortho <- function(ortho, tree, metadata, plas_ortho = NULL, num_ortho = 10, delete_genomes = 0,
                       tip_lab = T){

  element <- ortho$ortho_count %>%
    column_to_rownames("Orthogroup") # polarisのorthogorupsからorthogroups genecountを抽出
  #あるなしに変更
  # できれば、各Orthologの数の分布をしたにだしたいところ
  element[element > 1] <- 1

  # ゲノムの欠損による可視化ができない場合には、穴抜けになっている遺伝子を削除
  if (delete_genomes != 0){
    genome_num <- ncol(element)
    element <- element[rowSums(element) != (genome_num - delete_genomes),]
  }

  cols <- colnames(element)

  # Orthogroupsの有無によって、グルーピングを行っている
  element_data <- element %>%
    rownames_to_column("OG") %>%
    group_by(across(all_of(cols))) %>%
    summarise(freq = n()) %>%
    ungroup() %>%
    rowid_to_column()

  element_filtered <- element_data %>%
    arrange(desc(freq)) %>%
    filter(freq > num_ortho) # ある一定数以上の表示

  gg_element <- element_filtered %>%
    dplyr::mutate(freq2 = as.factor(freq)) %>%
    transform(freq3 = freq * -1) %>%
    tidyr::pivot_longer(c(-rowid, -freq, -freq2, -freq3), values_to = "presence") %>%
    full_join(metadata, by = c("name" = "genome")) # メタデータ対応、hostが色付けのまま
    # 上記については、なぜか-が.に変換されていることがあるので、注意すること

  if (!is_null(plas_ortho)){ #もしplasmid特有とされるOrthologを指定した場合にplasmidとchromosomeに載っている遺伝子を棒グラフで分けて表示する
    rowid_OG <- element %>%
      rownames_to_column("OG") %>%
      group_by(across(all_of(cols))) %>%
      nest() %>%
      summarise(freq = n(), data = data) %>%
      ungroup() %>%
      rowid_to_column() %>%
      select(rowid, data, freq) %>%
      unnest() %>%
      group_by(rowid) %>%
      mutate(freq = n()) %>%
      ungroup() %>%
      dplyr::filter(freq > num_ortho)

    plasmid_freq <- rowid_OG %>%
      filter(OG %in% unique(plas_ortho$Orthogroup)) %>%
      group_by(rowid, freq) %>%
      summarise(n = n())

    # chromosomeデータの生成
    gg_sum_data <- gg_element %>%
      dplyr::distinct(rowid, freq2, .keep_all = TRUE) %>%
      select(rowid, freq, freq2, freq3) %>%
      mutate(type = "Chromosome")

    # plasmidデータの生成
    gg_sum_plasmid <- tibble(rowid = plasmid_freq$rowid,
                             freq = plasmid_freq$n,
                             freq3 = -plasmid_freq$freq,
                             freq2 = as_factor(plasmid_freq$freq)) %>%
      mutate(type = "Plasmid")

    # このあたりすごく汚いので、書き直しが必要
    gg_sum_data2 <- gg_sum_data %>%
      full_join(gg_sum_plasmid, by = c("rowid" = "rowid")) %>%
      mutate(freq.z =  case_when(!is.na(freq.y) ~ freq.x - freq.y,
                                 is.na(freq.y) ~ freq.x)) %>%
      select(rowid, freq2.x, freq3.x, type.x, freq.z) %>%
      rename(freq2 = freq2.x, freq3 = freq3.x, type = type.x, freq = freq.z)

    # plasmidとchromosomeのbarplotはここで作ってしまう
    ortholog_sum <- gg_sum_data2 %>%
      bind_rows(gg_sum_plasmid) %>%
      ggplot(aes(x = reorder(rowid, freq3), y = freq, fill = type)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label=ifelse(freq != 0, freq, '')), vjust=0) +
      theme_minimal() +
      labs(x= "", y = "The number of Orthogrups") +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank())

  }

  ortholog_presence <- gg_element %>%
    ggplot(aes(x = reorder(rowid, freq3), y = name, fill = host)) + # ここがhostになっている
    geom_point(aes(size=ifelse(presence==0, NA, presence)), #, fill = host),
               shape = 21) +
    labs( x= "Presence of Orthologus", y = "") +
    theme_minimal() +
    guides(size="none") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())

  # 図の上部の棒グラフを作る
  if (is_null(plas_ortho)){ #plasmidの指定を行った時はすでにあるので作らない
    ortholog_sum <- gg_element %>%
      dplyr::distinct(rowid,freq2, .keep_all = TRUE) %>%
      ggplot(aes(x = reorder(rowid, freq3), y = freq)) + # rowidでreorderすることで、数が同じものも表示
      geom_bar(stat = "identity") +
      geom_text(aes(label=freq), vjust=0) +
      theme_minimal() +
      labs(x= "") +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank())
  }

  # 系統樹のラベルを変更
  # これが一般的なエラーなのかを探る必要がある
  tree_tip <- tree$tip.label
  # ここはWSK特異的な変更
  tree$tip.label <- str_replace_all(tree_tip, "-", ".") #なぜかread.treeだと-になるので、.へ変換
  ortholog_tree <- ggtree(tree) + geom_treescale(x=0.05, y=0, offset=2, fontsize = 3) + xlim_tree(0.05)

  if (isTRUE(tip_lab)){
    ortholog_tree <- ortholog_tree + geom_tiplab()
  }

  ortho_tree <- ortholog_presence %>%
    insert_top(ortholog_sum, height = .5) %>%
    insert_left(ortholog_tree)

  return(ortho_tree)
}

#############################################
# 取れたゲノムの系統組成を示す
#############################################
# gtdbtkの結果を基にして、系統組成を作成する
# metadataがある時はmetadataを用いて分ける
#
#############################################
#' Barplot using gtdbtk data with metadata
#'
#' @return
#' @export
#'
#' @examples
plot_gtdbtk <- function(gtdbtk, metadata = NULL, category = NULL,
                        type = c("single", "ani", "pplacer"),
                        taxonomy = c("domain", "phylum", "class", "order", "family", "genus", "species"), ...){

  if(is.null(metadata)){
    stop("metadataが指定されていません")
  }
  if(!("sample" %in% colnames(metadata))){
    stop("medataにsample列がありません")
  }
  if(!(category %in% colnames(metadata))){
    stop("metadataに指定されたカテゴリ列がありません")
  }

  gtdb_taxonomy <- paste(taxonomy, "gtdb", sep = "_")

  # metadataとの統合
  # metadataはgenomeがuser_genome列, その他が適当な列名になっていることが条件
  # sample列が必要
  gtdbtk_meta <- gtdbtk %>%
    left_join(metadata, by = "user_genome") # x軸用のメタデータの取得

  gtdbtk_rate <- gtdbtk_meta %>%
    group_by(!!sym(gtdb_taxonomy), !!sym(category), sample) %>% # categoryやsampleが1種でも問題なし
    summarise(count_genome = n()) %>%
    arrange(desc(count_genome)) %>%
    group_by(sample) %>% # sampleに変更
    mutate(percentage = count_genome / sum(count_genome) * 100) %>%
    mutate(tax = !!sym(gtdb_taxonomy),
           category = !!sym(category)) %>%
    ungroup() %>%
    select(tax, category, percentage, count_genome, sample) %>%
    mutate(tax = case_when(tax == "" ~ "Undetermined",
                           TRUE ~ tax))

  # 表示する色の指定
  fill_order <- gtdbtk_meta %>%
    group_by(!!sym(gtdb_taxonomy)) %>%
    summarise(count_genome = n()) %>%
    arrange(desc(count_genome)) %>%
    mutate(percentage = count_genome / sum(count_genome) * 100) %>%
    mutate(tax = !!sym(gtdb_taxonomy)) %>%
    mutate(tax = case_when(tax == "" ~ "Undetermined", TRUE ~ tax)) %>%
    filter(tax != "Undetermined") %>%
    arrange(percentage) %>%
    select(tax) %>%
    unlist()

  fill_order <- c("Undetermined", fill_order)
  mod_colors <- colors[1:length(fill_order)]
  mod_colors[length(fill_order)] <- "#808080"

  # barplotの作成
  bac_bar <- ggplot(gtdbtk_rate, aes(x = sample, y = percentage,
                                         fill = factor(tax, fill_order))) +
    geom_bar(stat = "identity") + theme_minimal() +
    labs(y = "Relative abundance (%)") +
    scale_fill_manual(values=rev(mod_colors), name = taxonomy) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(fill = guide_legend(reverse = TRUE)) +
    facet_grid(. ~ category, scales = "free_x") # categoryが指定されている場合には分けて

  return(bac_bar)
}


#############################################
# Checkmの結果、いつもの図を作る
#############################################
# checkmの結果を基にしていつもの図を作る
# completenessとcontamination
#
#############################################
#' Title
#'
#' @return
#' @export
#'
#' @examples
plot_checkm <- function(checkm, metadata = NULL, order = NULL, ...){

  quality_gg <- checkm %>%
    ggplot() +
    geom_point(aes(x = reorder(Bin_Id, -Completeness), y = Completeness, colour= "Completeness")) +
    geom_point(aes(x = Bin_Id, y = Contamination, colour= "Contamination")) +
    theme_minimal() +
    scale_colour_manual(values = c("red", "blue")) +
    scale_y_continuous(sec.axis = sec_axis(~.*5, name = "Contamination [%]")) +
    labs(y = "Completeness [%]", x = "The number of SAGs", color = "type")


  print(quality_gg)
}

#############################################
# 系統樹とcheckm, gtdbtkの結果を統合したものを作成
#############################################
# ggtreeExtraを使った可視化をする
#############################################

#' Title
#'
#' @param phylo
#' @param gtdb
#' @param quality_df
#' @param checkm
#' @param metadata
#' @param legend.size
#'
#' @return
#' @export
#'
#' @examples
plot_phylo <- function(phylo, gtdb, checkm, metadata, legend.size = 0.8){

  # ggtree objectの作成
  # 基本設定は円形
  m_p <- ggtree(phylo, layout="fan", branch.length = "none",
                size=0.4, open.angle=10, aes(color=group)) +
    scale_color_manual(values=colors,
                       guide=guide_legend(keywidth = legend.size,
                                          keyheight = legend.size, order=4)) +
    new_scale_color()

  # gtdbtkのannotationを付与
  m_p1 <- m_p %<+%
    gtdb + geom_fruit(data=metadata, geom=geom_tile,
                      mapping=aes(y=user_genome,
                                  fill=spcimen_id),
    ) +
    scale_fill_manual(values=colors,
                      guide=guide_legend(keywidth = legend.size,
                                         keyheight = legend.size, order=4)) +
    new_scale_fill() +
    geom_fruit(data=quality_df, geom=geom_tile,
               mapping=aes(y=Bin_Id, fill=quality),
    ) +
    scale_fill_manual(values=colors,
                      guide=guide_legend(keywidth = legend.size,
                                         keyheight = legend.size, order=4)) +
    new_scale_fill()

  m_p2 <- m_p1 + geom_fruit(data=checkm, geom=geom_bar,
                            mapping=aes(y=Bin_Id, x=Completeness, fill=species_gtdb),
                            orientation="y",
                            stat="identity",
                            axis.params=list(
                              axis       = "x",
                              text.size  = 1,
                              hjust      = -0.5,
                              text.angle  = 315,
                            )
  ) +
    scale_fill_manual(values=colors,
                      guide=guide_legend(keywidth = legend.size,
                                         keyheight = legend.size, order=4))+
    geom_fruit(data=checkm, geom=geom_bar,
               mapping=aes(y=Bin_Id, x=Contamination, fill=species_gtdb),
               orientation="y",
               stat="identity",
               axis.params=list(
                 axis       = "x",
                 text.size  = 1,
                 hjust      = -0.5,
                 text.angle  = 315,
               )
    ) +
    scale_fill_manual(values=colors,
                      guide=guide_legend(keywidth = 0.8,
                                         keyheight = 0.8, order=4))+
    geom_treescale(fontsize=2, linesize=0.3, x=4.9, y=0.1) +
    theme(legend.position=c(0.95, 0.5),
          legend.background=element_rect(fill=NA),
          legend.title=element_text(size=6.5),
          legend.text=element_text(size=4.5),
          legend.spacing.y = unit(0.02, "cm"),
    )
  return(m_p2)
}

#############################################
# ANIの結果を可視化する
#############################################
# ゲノム間のANIの結果を可視化する
# defaultはヒートマップになっている
#
# optionによって、サブクラスごとのANIの分布にすることも可能
#
#############################################
#' Title
#'
#' @param polaris
#' @param plot_type
#' @param subcluster
#' @param ...
#'
#' @return
#' @export
#'
#' @importFrom pheatmap pheatmap
#' @importFrom purrr is_null
#'
#' @examples
plot_ani <- function(ani, genome_list = NULL, summary_df = NULL,
                     completeness = 50, contamination = 10,
                     annotate_df, subcluster = FALSE,
                     subcluster_methods = c("hclust_complete", ""), ...){

  if(!is_null(genome_list)){
    ani <- aai[genome_list,genome_list]
  }

  hclust_obj <- (100 - ani) %>%
    dist() %>%
    hclust(method = "complete")

  order <-  ccsag_se_hclust_obj %>% .$order
  order2 <- rownames(ani[order,])


  ani_heatmap <- pheatmap(ani[order2,order2], show_rownames = F, show_colnames = F,
                       cluster_rows = F, cluster_cols = F, annotation_col = annotate_df)

  return(ani_heatmap)
}

#############################################
# AAIの結果を可視化する
#############################################
# ゲノム間のAAIの結果を可視化する
# defaultはヒートマップになっている
#
# optionによって、サブクラスごとのAAIの分布にすることも可能
#
#############################################
#' Title
#'
#' @param polaris
#' @param plot_type
#' @param subcluster
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_aai <- function(polaris, plot_type = c("heatmap", "distribution"),
                     subcluster = FALSE, ...){

  aai <- polaris # aaiの結果を抽出する

}

#############################################
# ゲノム間の距離を比較する
#############################################
# ゲノム間の距離を比較する
# 比較できるのは
# 系統距離: phyloから
# ANI: aniから
# AAI: aaiから
# gene content: Orthogroupsから
#############################################
# GGallyみたいな比較の図を書きたいね
#############################################

plot_distance <- function(polaris){

  aai <- polaris # aaiの結果を抽出する

}

#############################################
# thresholdを設定した際に、どの様なcluster数になるかplotする
#############################################
#
#############################################

