#############################################
# fastaファイルをまとめて読み込む
#############################################
# 関数をまとめて読みこませるgeneral関数を作った方がいい気がする
#############################################
#' Read fasta files from specific directory
#'
#' fastaファイルをまとめて読み込むための関数
#'
#' @param fastas fastaのファイルパスのリストもしくはディレクトリ
#' @param dir (default: FALSE)
#' @param type 配列のタイプ, DNAかAAか
#' @param extension fastaファイルの拡張子。読み込みの際に指定された拡張子はファイル名から削除されて、ゲノム名となる
#' @param names ゲノム の名前を自分で指定したい場合に指定する
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom stringr str_split
#' @importFrom stringr str_replace_all
#' @importFrom purrr map
#'
#' @return DNAStringSetと名前のリスト
#' @export
#'
#' @examples
#' #fasta_dir <- system.file("inst/extdata/Rickettsiales/fasta", package = "polaris")
#' #fasta_list <- read_fastas(fasta_dir, dir = TRUE, type = "DNA")
#'
#' # fastaのファイル一覧取得。こちらでも動く
#' # fasta_filelist <- list.files(path = fastas, pattern = NULL, full.names = T)
#' # fasta_list <- read_fastas(fasta_filelist)
read_fastas <- function(fastas, dir = FALSE, type =c("DNA", "AA"),
                        extension = "fna", names = NULL, ...){

  # ディレクトリか否か
  if (isTRUE(dir)){
    # ファイルの個々のリストを取得する
    fasta_list <- list.files(path = fastas, pattern = NULL, full.names = T)
  } else {
    fasta_list <- fastas
  }

  # ゲノムの名前を自分で指定するかどうか
  # 要素の数で確認はしていないので、少し挙動が不安
  if(!is_null(names)){
    genome_names <- names
  } else {
    name_list <- fasta_list %>% str_split("/", simplify = T)
    genome_names <- name_list[,ncol(name_list)] %>%
      str_replace_all(pattern=paste0(".",extension), replacement="") %>%
      str_replace_all(pattern=" ", replacement="_")
  }

  # BiostringsのreadDNAまたはreadAAを用いてまとめてmapで取得
  # 中身はBiostringのXstring classのlist
  genomefastas <- fasta_list %>%
    map(readDNAStringSet)

  names(genomefastas) <- genome_names
  genome_list <- list(genome = genomefastas, names = genome_names)
  # fastaファイル自体とそれぞれのgenomeの名前を定義する
  # listの要素名に指定する
  return(genome_list)
}

#############################################
# gffの読み込み
#############################################
#
# gffファイルを読み込む
# 現在は真核生物は非対応。エクソンの処理ができない
# 将来的に対応してもいいがそれよりやることは多い
#
#############################################
#' Read gff files from specific directory
#'
#' gffファイルを一つ読み込むための関数
#'
#' @param file gffファイルのパス
#' @param tibble (default: FALSE) 出力をgffclassではなくtibbleにする
#'
#' @importFrom readr read_tsv
#' @importFrom tidyr separate
#' @importFrom tidyr unnest
#' @importFrom tidyr pivot_wider
#'
#'
#' @return gff gffのクラス
#' @export
#'
#' @examples
#'
read_gff <- function(file, tibble = FALSE, extension = "gff"){

  colnames <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase","attributes")
  coltype <- list("c", "c", "c", "d", "d", "c", "c", "c", "c") # 一文字じゃないといけない

  name <- file %>% str_split("/", simplify = T)
  genome_names <- name[,ncol(name)] %>%
    str_replace(pattern=paste0(".",extension), replacement="") %>%
    str_replace(pattern=" ", replacement="_")

  # 行の頭が#の場合は適切な形で処理する
  # prokkaのgffには末尾にfastaファイルがあるらしい (なんで？？？)
  # これらの処理をするためにna.omitしている。もっといい方法があるように思うが、、、
  gff <- read_tsv(pipe(paste("grep -v \\#", file)), col_names = colnames, col_types = coltype) %>% # 行の頭が#の場合は適切な形で処理する
    na.omit() %>% #warningの原因、fasta部分を削除
    mutate(attributes_list = strsplit(attributes, ";")) %>% # tag情報などを分割して
    unnest(cols = attributes_list) %>% #一行づつの対応にして、
    separate(attributes_list, into = c("attributes_id", "elements"), sep = "=") %>% # 列名となるものと要素を分解し
    pivot_wider(names_from = attributes_id, values_from = elements, values_fill = NA) # 列名と要素をいれる

  gff_name <- gff %>%
    mutate(genome = genome_names)

  # 返り値はgffのデータ
  return(gff_name)
}

#############################################
# gffの読み込み
#############################################
#
# まとめて、gffを読み込む
#
#############################################
#' 複数のgffファイルの読み込み
#'
#' 複数のgffファイルを一つ一つ読み込んで統合することが少々煩雑なので、
#' 最初から複数ファイルを読み込むための関数を提供する
#'
#' @importFrom
#'
#' @return
#' @export
#'
#' @examples
read_gffs <- function(files, dir = FALSE, extension = "gff", names = NULL){

  # ディレクトリか否か
  if (isTRUE(dir)){
    # ファイルの個々のリストを取得する
    file_list <- list.files(path = files, pattern = NULL, full.names = T)
  } else {
    file_list <- files
  }

  # gffの名前を自分で指定するかどうか
  # 要素の数で確認はしていないので、少し挙動が不安
  if(!is_null(names)){
    gff_names <- names
  } else {
    name_list <- file_list %>% str_split("/", simplify = T)
    gff_names <- name_list[,ncol(name_list)] %>%
      str_replace_all(pattern=paste0(".",extension), replacement="") %>%
      str_replace_all(pattern=" ", replacement="_")
  }

  # BiostringsのreadDNAまたはreadAAを用いてまとめてmapで取得
  # 中身はBiostringのXstring classのlist
  gff_output <- file_list %>%
    map(read_gff) %>%
    bind_rows()

  # # gffファイル自体とそれぞれのgenomeの名前を定義する
  # names(gff_list) <- gff_names
  #   # listをgenomefastaと定義し直して、return
  # gff_output <- list(gff = gff_list, names = gff_names)

  return(gff_output)
}

#############################################
# Orthofinderの結果を読み込む
#############################################
# Orthogroups
#############################################
#' Orthofinderの結果を読み込む
#'
#' Orthofinderの結果を読み込む。出力としては3つある。
#' 1つ目はどの遺伝子とどのOrthogroupが対応しているか
#' 2つ目はどのOrthogroupがどのゲノムにどれくらいあるか
#' 3つ目はシングルコピーオルソログはどれか
#'
#' @param ortho_file orthofinderの結果、orthofroups.tsv
#'
#' @return ortho_class ortho
#' @export
#'
#' @examples
read_ortho <- function(ortho_file){
  ortho_tb <- read_tsv(ortho_file)

  # カンマを含んでいる全ての要素はリストに分解したい
  ortho_list <- ortho_tb %>%
    mutate_all(~str_replace_all(.," ", "")) %>% #Use `mutate_at(df, vars(-group_cols()), myoperation)` to silence the message.
    mutate_all(~str_split_list(.,pattern = ",")) %>% #Use `mutate_at(df, vars(-group_cols()), myoperation)` to silence the message.
    unchop(Orthogroup) %>%
    pivot_longer(-Orthogroup) %>%
    unnest(-Orthogroup) %>%
    na.omit()

  ortho_count <- ortho_list %>%
    group_by(Orthogroup, name) %>%
    mutate(count = n()) %>%
    select(-value) %>%
    distinct() %>% #重複を削除
    pivot_wider(id_cols = Orthogroup, names_from = name, values_from = count) %>% #行をOG, 列をgenomeに
    mutate_all(~replace(., is.na(.), 0)) #NAを置換

  singlecopy <- ortho_count %>% #single copy marker geneの抽出が必要
    column_to_rownames("Orthogroup") %>%
    filter_all(all_vars(. == 1)) %>% #全ての要素が1であれば
    rownames()

  ortho <- list(ortho_genes = ortho_list, ortho_count = ortho_count, singlecopy = singlecopy)
  return(ortho)
}

#############################################
# Checkmの結果を読み込む
#############################################
# checkmの結果を読み込む
# tsvがinputであれば、TRUEとしてそうでなければ、
# tsvじゃないものをinput (ちょっと面倒くさい)
#############################################
#############################################

#' read checkm
#'
#' checkmの結果を読み込む
#'
#' @param file
#'
#' @return
#' @export
#'
#' @examples
read_checkm <- function(file){

  # gtdbtkのNAをRのNAへ変換する
  checkm_data <- read_tsv(file) %>%
    rename_all(~str_replace_all(.,"#","")) %>%
    rename_all(~str_replace_all(.,"^ ","_")) %>%
    rename_all(~str_replace_all(.," ","_")) # closest_placement

  return(checkm_data)

}

#############################################
# GTDBtkの結果を読み込む
#############################################
# gtdbtk.bac120.summary.tsvをinputにする
# 基本的にはtaxonomyは分解される
#
#############################################
#' GTDBtkのtsvを読み込む
#'
#' @param file gtdbtk.bac120.summary.tsv
#'
#' @return
#' @export
#'
#' @examples
#' #file <- system.file("inst/extdata/Rickettsiales/gtdbtk/gtdbtk.bac120.summary.tsv", package = "polaris")
#' #gtdb <- read_gtdbtk(file)
#'
#' # polaris %<>%
#' #  add_polaris(gtdb)
#' #
read_gtdbtk <- function(file){

  # gtdbtkのNAをRのNAへ変換する
  gtdb_data <- read_tsv(file) %>%
    na_if("N/A")

  taxonomy_names <- c("domain", "phylum",
                      "class", "order",
                      "family", "genus",
                      "species")

  refined_gtdb <- gtdb_data %>%
    separate(classification, paste(taxonomy_names, "gtdb", sep = "_"), sep = ";") %>% # classification
    separate(fastani_taxonomy, paste(taxonomy_names, "ani", sep = "_"), sep = ";") %>% # fastani_taxonomy
    separate(closest_placement_taxonomy, paste(taxonomy_names, "place", sep = "_"), sep = ";")  %>% #pplacer
    mutate_all(~str_replace_all(.,"^.__","")) # closest_placement

  return(refined_gtdb)
}

#############################################
# GTDBtkの結果を読み込む
#############################################
# gtdbtk.bac120.summary.tsvをinputにする
# 基本的にはtaxonomyは分解される
#
#############################################
#' Eggnong annotationの読み込み
#'
#' Eggnogデータを操作しやすい形で読み込む
#' おそらく、懸念点としては、各cell内でカンマで
#' 分けられてるデータをリスト化すること
#'
#' @param file eggnongの結果、eggnog annotation
#'
#' @return
#' @export
#'
#' @examples
read_eggnog <- function(file, name = NULL, extension = ".annotations", tools = c("eggnog", "atlas")){

  if (tools == "eggnog"){
    eggnog <- read_tsv(file, skip = 4) %>% # versionによって数が変わる
      rename_all(~str_replace_all(.," ","_"))   # スペースを_へ
  } else {
    eggnog <- read_tsv(file)
  }

  if(is_null(name)){
    name <- file %>% str_split("/", simplify = T)
    genome_names <- name[,ncol(name)] %>%
      str_replace(pattern=paste0(".",extension), replacement="") %>%
      str_replace(pattern=" ", replacement="_")
  } else {
    genome_names <- name
  }

  eggnog_splited_ids <- c("GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
                   "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy",
                   "BiGG_Reaction", "eggNOG_OGs", "PFAMs", "COG_category")
  atlas_splited_ids <- c("GO_terms", "EC", "KO", "KEGG_Pathway", "KEGG_Module",
                         "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy",
                         "BiGG_Reaction", "Eggnog", "COG_cat", "PFAMs")
  # 各列において、カンマを含む列を文字列からリストへ
  # acrossで書き換える
  if (tools == "eggnog"){
    refined_eggnong <- eggnog %>%
      mutate(across(eggnog_splited_ids, str_split_list)) %>% # acrossを活用
      mutate(genome = genome_names)
  } else {
    refined_eggnong <- eggnog %>%
      mutate(across(atlas_splited_ids, str_split_list)) %>% # acrossを活用
      mutate(genome = genome_names)
  }

  return(refined_eggnong)
}

# カンマを含む要素をリストへ分解する
str_split_list <- function(cols, pattern = ","){
  cols %>%
    str_split(pattern = pattern) %>%
    as.list()
}

#############################################
# GTDBtkの結果を読み込む
#############################################
# gtdbtk.bac120.summary.tsvをinputにする
# 基本的にはtaxonomyは分解される
#
#############################################
#' Eggnogデータのまとめて読み込み
#'
#' Eggnong-mapperの結果をまとめて読み込む
#'
#' @param eggnog
#' @param dir
#'
#' @return
#' @export
#'
#' @examples
#'
#'
read_eggnogs <- function(eggnog, dir = FALSE, extension = "annotations", names = NULL, tools = c("eggnog", "atlas")){

  # ディレクトリか否か
  if (isTRUE(dir)){
    # ファイルの個々のリストを取得する
    eggnog_files <- list.files(path = eggnog, pattern = NULL, full.names = T)
  } else {
    eggnog_files <- eggnog
  }

  eggnog_list <- eggnog_files %>%
    map(read_eggnog, name = names, extension = extension, tools = tools) %>%
    bind_rows()

  return(eggnog_list)
}

#############################################
# Interproscanの結果を読み込む
#############################################
# Inteproscan
# Protein accession (e.g. P51587)
# Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
# Sequence length (e.g. 3418)
# Analysis (e.g. Pfam / PRINTS / Gene3D)
# Signature accession (e.g. PF09103 / G3DSA:2.40.50.140)
# Signature description (e.g. BRCA2 repeat profile)
# Start location
# Stop location
# Score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
# Status - is the status of the match (T: true)
# Date - is the date of the run
# InterPro annotations - accession (e.g. IPR002093)
# InterPro annotations - description (e.g. BRCA2 repeat)
# (GO annotations (e.g. GO:0005515) - optional column; only displayed if –goterms option is switched on)
# (Pathways annotations (e.g. REACT_71) - optional column; only displayed if –pathways option is switched on)
#############################################
#' interproscanの読み込み
#'
#' 基本的にはtsvを読み込むことにする
#' 複数のゲノムがある場合には、
#'
#' @param file interproscanファイル
#' @param type ファイルのタイプ3パターン
#'
#' @return
#' @export
#'
#' @examples
read_interpro <- function(file, type = c("tsv", "json", "xml")){

  if(type == "tsv") {
    colnames <- c("protein_accession", "md5", "length",
                  "database", "db_accession", "db_annotation",
                  "start", "end", "e_values", "status",
                  "date", "interpro_accession", "interpro_desc",
                  "GO", "Pathway")

    interpro <- read_tsv(file, col_names = colnames) %>% #1行1タンパク質にしたい
      group_nest(protein_accession, md5, length) #何でnestするのかは要検討

  } else if (type == "json") {

  } else if (type == "xml") {

  } else {
    stop("読み込みの形式が不明です。指定してください")
  }

  return(interpro)
}

#############################################
# interproscanのファイルをまとめて読み込み
#############################################
# interproscanの結果をまとめて読み込むための
# ツール。
#############################################
#' Title
#'
#'
#' @return
#' @export
#'
#' @examples
#'
#'
read_interpros <- function(files, dir = FALSE, type = c("tsv", "json", "xml"),
                           names = NULL){

  # ディレクトリか否か
  if (isTRUE(dir)){
    # ファイルの個々のリストを取得する
    file_list <- list.files(path = files, pattern = ".tsv", full.names = T)
  } else {
    file_list <- files
  }

  interpro_data <- file_list %>%
    map(read_interpro, type = type) %>%
    bind_rows()

  return(interpro_data)
}

read_ani <- function(file, extension = "fna", program = c("fastani", "enveomics")){

  ani_colnames <- c("query", "target", "ANI", "tmp1", "tmp2")

  ani <- file %>%
    read_tsv() %>%
    mutate(seq1 = str_replace(seq1, pattern=paste0(".",extension), replacement=""),
           seq2 = str_replace(seq2, pattern=paste0(".",extension), replacement="")) %>%
    pivot_wider(names_from = seq1, values_from = ANI) %>%
    column_to_rownames("seq2")

  return(ani)
}

read_aai <- function(file, extension = "faa"){

  aai <- file %>%
    read_tsv() %>%
    mutate(seq1 = str_replace(seq1, pattern=paste0(".",extension), replacement=""),
           seq2 = str_replace(seq2, pattern=paste0(".",extension), replacement="")) %>%
    pivot_wider(names_from = seq1, values_from = AAI) %>%
    column_to_rownames("seq2")

  return(aai)
}

read_tree <- function(file){

  phylo <- file %>%
    read_tree()

  return(phylo)

}

#############################################
# blastやmmseqs2,diamondの結果を読み込む
#############################################
# 毎回読み込ませるのがめんどくさいので
#############################################
#' Read blast format 6 file
#'
#' @param file
#' @param format
#'
#' @importFrom readr read_tsv
#' @return
#' @export
#'
#' @examples
read_blast <- function(file, format = 6){
  colnames <- c("qseqid", "sseqid", "pident", "length",
                "mismatch", "gapopen", "qstart", "qend",
                "sstart", "send", "evalue", "bitscore")

  blast <- read_tsv(file, col_names = colnames)

  return(blast)
}

#############################################
# サンプルデータの読み込み
#############################################
# サンプルメタデータを読み込ませる
# colnamesは
# 1. genome
# 2. 以下任意 (想定としてはsubclassなどが入る)
#
#############################################
#' metaデータの読み込み
#'
#' ゲノムのメタデータの読み込みを行う
#'
#' @param file tsvまたはcsvファイル
#'
#' @return sample_data
#' @export
#'
#' @examples
#'
#'
read_metadata <- function(file){

  sample_data <- read_tsv(file)

  if (TRUE) { #genome列が必要。列名に含んでいるかどうかを確かめる
    stop("genome列がありません。メタデータにはgenome列が必須です")
  }

  return(sample_data)

}
