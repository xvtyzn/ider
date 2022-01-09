#############################################
# orthologの結果を表示させる
#############################################
# 系統樹と対応して、orthologの有無の表示
# いつものやつ。aplotで作成する
#############################################
#' Title
#'
#' @param polaris polaris object
#' @param label
#' @param num_ortho
#' @param
#'
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom aplot
#'
plot_ortho <- function(polaris, label, num_ortho){

  host <- label #メタデータから取得するようにしたいよね

  element <- polaris # polarisのorthogorupsからorthogroups genecountを抽出
  #あるなしに変更
  # できれば、各Orthologの数の分布をしたにだしたいところ
  element[element > 1] <- 1
  cols <- colnames(element)

  element_filtered <- element %>%
    rownames_to_column("OG") %>%
    group_by(across(all_of(cols))) %>%
    summarise(freq = n()) %>%
    ungroup() %>%
    rowid_to_column() %>%
    dplyr::arrange(desc(freq)) %>%
    dplyr::filter(freq > num_ortho) # ある一定数以上の表示

  gg_element <- element_filtered %>%
    dplyr::mutate(freq2 = as.factor(freq)) %>%
    transform(freq3 = freq * -1) %>%
    tidyr::pivot_longer(c(-rowid, -freq, -freq2, -freq3, -Total), values_to = "presence") %>%
    dplyr::mutate(host = rep(host, 43)) # ここをどうにかして、自動化しないと

  gg_element %>%
    ggplot(aes(x = reorder(rowid, freq3), y = name)) +
    geom_point(aes(size=ifelse(presence==0, NA, presence), fill = host),
               shape = 21) +
    labs( x= "Presence of Orthologus", y = "") +
    theme_minimal() +
    guides(size=FALSE) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank()) -> ortholog_presence

  gg_element %>%
    dplyr::distinct(rowid,freq2, .keep_all = TRUE) %>%
    ggplot(aes(x = reorder(rowid, freq3), y = freq)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label=freq), vjust=0) +
    theme_minimal() +
    labs(x= "") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank()) -> ortholog_sum


  ortho_tree <- ortholog_presence %>%
    insert_top(ortholog_sum, height = .5) %>%
    insert_left(ortholog_tree)
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
    group_by(!!sym(gtdb_taxonomy), !!sym(category), sample) %>%
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
  bac_bar <- ggplot(gtdbtk_rate, aes(x = sample, y =percentage,
                                         fill = factor(tax, fill_order))) +
    geom_bar(stat = "identity") + theme_minimal() +
    labs(y = "Relative abundance (%)") +
    scale_fill_manual(values=rev(mod_colors_bac)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(fill = guide_legend(reverse = TRUE)) +
    facet_grid(. ~ category, scales = "free_x")

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

plot_qcphylo <- function(polaris, layout = "fan",
                         legend.size = 0.8,
                         open.angle = 10,
                         fill.element = NULL,
                         bar.element = NULL,
                         phylo_color = NULL,
                         bar_color = NULL,
                         metadata_column = NULL){



}

# 以下のはdraft versionこれを上記の形に変えていく
plot_overview <- function(phylo, gtdb, quality_df, checkm, metadata, legend.size = 0.8){
  m_p <- ggtree(phylo, layout="fan", branch.length = "none",
                size=0.4, open.angle=10, aes(color=group)) +
    scale_color_manual(values=colors,
                       guide=guide_legend(keywidth = legend.size,
                                          keyheight = legend.size, order=4)) +
    new_scale_color()

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
#' @examples
plot_ani <- function(polaris, plot_type = c("heatmap", "distribution"),
                     subcluster = FALSE, ...){

  ani <- polaris # aniのデータを抽出する



  print()

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

