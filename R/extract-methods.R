#############################################
# ortholog解析の結果から特定の分布を示すortholog IDを取得する
#############################################
# ゲノムリストやfrequencyを指定することで、それらのゲノムが
# 特異的に持つOrthogroupsを抽出する
#
#############################################
#' Title
#'
#' @param ortho
#' @param freq
#' @param genome_list
#' @param delete_genomes
#'
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom purrr is_null
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr ungroup
#' @importFrom dplyr select
#' @importFrom tidyr nest
#'
extract_ortho <- function(ortho, freq = NULL, genome_list = NULL, delete_genomes = 0){

  element <- ortho$ortho_count %>%
    column_to_rownames("Orthogroup") # polarisのorthogorupsからorthogroups genecountを抽出
  #あるなしに変更
  element[element > 1] <- 1

  # ゲノムの欠損による可視化ができない場合には、穴抜けになっている遺伝子を削除
  if (delete_genomes != 0){
    genome_num <- ncol(element)
    element <- element[rowSums(element) != (genome_num - delete_genomes),]
  }

  # freqで選択するか genome_listで選択するかの判断が必要
  if(is_null(freq) & is_null(genome_list)){
    stop("freqかgenome_listに値を入れてください")
  } else if(!is_null(freq) & !is_null(genome_list)){
    stop("freqとgenome_list両方の指定はできません")
  }

  if(!is_null(genome_list)){
    ortho_list <- element %>%
      filter_except_zero(genome_list) %>%
      rownames_to_column("OG") %>%
      select(OG)  %>%
      unlist()
  } else if (!is_null(freq)){
    ortho_list <- element %>%
      rownames_to_column("OG") %>%
      group_by(across(all_of(cols))) %>%
      nest() %>%
      mutate(freq = map(data, nrow)) %>%
      filter(freq == freq) %>%
      ungroup() %>%
      select(data)  %>%
      unlist()
  }

  return(ortho_list)
}

#############################################
# 指定したidが1で、それ以外が0となっている要素を取得する
#############################################
# 指定したゲノムが1であり、それ以外が0となっている要素をfilterする
# filter_atは非推奨になってしまったので、使わない
# outputしない
#############################################
#' Title
#'
#' @param dat
#' @param ids
#'
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom dplyr filter
#' @importFrom dplyr if_all
#'
filter_except_zero <- function(dat, ids){ #なぜこれがupdateされないのか

  expected_dat <- filter(dat, if_all(!ids, ~ .x == 0) & if_all(ids, ~ .x != 0))

  return(expected_dat)
}

#############################################
# ortholog解析の結果からOrtholog多い順にゲノムリストを取得
#############################################
# orthlogリストからゲノムのリストを取得する
#
#
#############################################
#' Title
#'
#' @param ortho
#' @param delete_genomes
#'
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom purrr is_null
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr ungroup
#' @importFrom dplyr select
#' @importFrom dplyr across
#' @importFrom dplyr all_of
#' @importFrom dplyr summarise
#' @importFrom dplyr arrange
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr nest
#' @importFrom
#'
extract_genomelist <- function(ortho, delete_genomes = 0, num_ortho = 10){

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
  element_filtered <- element %>%
    rownames_to_column("OG") %>%
    group_by(across(all_of(cols))) %>%
    summarise(freq = n()) %>%
    ungroup() %>%
    rowid_to_column() %>%
    arrange(desc(freq)) %>%
    filter(freq > num_ortho) # ある一定数以上の表示

  gg_element <- element_filtered %>%
    mutate(freq2 = as.factor(freq)) %>%
    transform(freq3 = freq * -1) %>%
    pivot_longer(c(-rowid, -freq, -freq2, -freq3), values_to = "presence")

  gg_element2 <- gg_element %>%
    group_by(rowid) %>%
    nest(c(name)) %>%
    filter(presence == 1) %>%
    select(rowid, freq, data) %>%
    ungroup()

  return(gg_element2)
}

#############################################
# ortholog解析の結果から特定の分布を示すortholog IDを取得する
#############################################
# ゲノムリストやfrequencyを指定することで、それらのゲノムが
# 特異的に持つOrthogroupsを抽出する
#
#############################################
#' Title
#'
#' @param tmp
#' @param ortho_data
#' @param delete_genomes
#'
#' @return
#' @export
#'
#' @examples
extract_ortho2 <- function(tmp, ortho_data, delete_genomes = 0){
  tmp2 <- str_replace(tmp, "\\.", "-")
  aaa <- extract_ortho(ortho_data,
                       genome_list = tmp2,
                       delete_genomes = delete_genomes)
  return(aaa)
}

#############################################
# ortholog解析の結果から特定の分布を示すortholog IDを取得する
#############################################
# ゲノムリストやfrequencyを指定することで、それらのゲノムが
# 特異的に持つOrthogroupsを抽出する
#
#############################################
# Orthologから遺伝子のみを抽出
#' Title
#'
#' @param data
#' @param ortho_data
#'
#' @return
#' @export
#'
#' @examples
extract_gene <- function(data, ortho_data){
  ortho_genes <- ortho_data$ortho_genes %>%
    filter(Orthogroup %in% data) %>%
    select(value) %>%
    arrange(desc(value))
  return(ortho_genes)
}

#############################################
# 遺伝子名でeggnogを抽出し、指定したカテゴリのデータをunchopした状態にする
#############################################
# 遺伝子名のベクトルを得た後に、その遺伝子名で
# eggnongデータからアノテーションを取得する
#############################################
#' Title
#'
#' @param data
#' @param category
#' @param eggnog_data
#'
#' @return
#' @export
#'
#' @examples
extract_eggnog <- function(data, category, eggnog_data){
  eggnog_annot <- eggnog_data %>%
    filter(`#query` %in% data) %>%
    unchop(category)
  return(eggnog_annot)
}
