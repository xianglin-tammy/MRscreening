#' This is some descriptio of this function.这是这个函数的一些描述。
#' @title query wrapper function 查询包装器功能
#'
#' @description query wrapper function查询包装器功能
#'
#' @details This is a general purpose function to send data request which can be used when there has not been an R equivalent package function to an API endpoint.
#这是一个通用函数，用于发送数据请求，当没有R等效包函数到API端点时可以使用。
#' @importFrom epigraphdb query_epigraphdb
#' @param query A list of parameters associated with the query endpoint.
#' 与查询端点相关联的参数列表。
#' @return a dataframe
#' @export


query_epigraphdb_as_table <- function(query){
  results_subset <- epigraphdb::query_epigraphdb(
    route = "/cypher",
    params = list(query = query),
    method = "POST",
    mode = "table")
}
#’这是一个用于查询EpigraphDB的函数。它接受一个查询作为参数，并返回查询结果的一个子集，以表格的形式展示出来。
#’函数使用epigraphdb包中的query_epigraphdb函数来执行查询操作。查询通过向EpigraphDB的/cypher路由发送POST请求
#’进行处理，并通过设置mode参数为"table"将查询结果以表格的形式返回。
#’将这个功能设置名字为query_epigraphdb_as_table
#' This is some descriptio of this function.
#' @title query wrapper function
#'
#' @description query wrapper function
#'
#' @details query wrapper function
#' @importFrom kableExtra kable kable_styling
#' @importFrom dplyr %>%
#' @param df
#' @return a dataframe
#' @export
#'这是一个用于将数据框转换为漂亮的表格并添加样式的函数。函数名为kable_it，它接受一个数据框(df)作为参数。
#'函数使用了kableExtra包中的kable函数将数据框转换为表格。然后，使用kable_styling函数来为表格添加样式，
#'以使其更具可读性和美观性。
#'最终，函数将返回一个带有样式的漂亮表格。您可以将该表格用于数据展示、报告生成等目的。

kable_it<-function(df){

  df %>%
    kableExtra::kable(.) %>%
    kableExtra::kable_styling()
}

#' This is some descriptio of this function.
#' @title query wrapper function
#'
#' @description make forest for query result为查询结果创建森林图
#'
#' @details make forest for query result
#' @importFrom ggplot2 ggplot aes geom_errorbarh geom_point geom_text theme_light scale_color_manual scale_y_discrete geom_vline element_text theme facet_wrap labs
#' @importFrom stats reorder数据重新排序
#' @importFrom wesanderson wes_palette
#' @param res_sub mr result
#' @param my_title the title name of my plot
#' @return a dataframe
#' @export


makeforest_plot <- function(res_sub,
                            my_title){

  pal<-c(wesanderson::wes_palette("Zissou1"))[c(1,3,5)]
  cols <- c("negative" = pal[1], "overlaps null" = pal[2], "positive" = pal[3])

  p<-ggplot2::ggplot(res_sub,
                     ggplot2::aes(y=stats::reorder(outcome.details, -or), x=or, label=outcome.details, colour=effect_direction)) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin=or_loci, xmax=or_upci), height=.3) +
    ggplot2::geom_point(size=2)+
    ggplot2::geom_text(ggplot2::aes(label=OR_CI),hjust=-0.3, vjust=-0.1, size =3, color = 'darkgrey')+
    ggplot2::theme_light()+
    ggplot2::scale_color_manual(values=cols)+
    ggplot2::scale_y_discrete(position = "left")+
    ggplot2::geom_vline(xintercept=1, linetype='longdash') +
    ggplot2::theme(strip.text = ggplot2::element_text(face = 'bold'))+
    ggplot2::facet_wrap(~outcome, scales = 'free_y', ncol = 1) +
    ggplot2::labs(color = "",y = "", x = "Odds ratio", subtitle="",
                  title=paste0("         ", my_title))+
    ggplot2::theme(legend.position = "none")
  return(p)
}

#' This is some descriptio of this function.
#' @title query wrapper function
#'
#' @description run local mr all outcomes
#'
#' @details run local mr all outcomes
#' @importFrom TwoSampleMR extract_instruments extract_outcome_data harmonise_data mr generate_odds_ratios split_outcome split_exposure
#' @importFrom tidyr separate
#' @importFrom dplyr bind_rows filter rename mutate
#' @param trait_investigated the trait id
#' @param outcomes_to_try the outcome id
#' @return a dataframe
#' @export

run_local_mr_all_outcomes <- function(trait_investigated, outcomes_to_try){

  instruments <- TwoSampleMR::extract_instruments(trait_investigated)

  extract_outcome_data_custom <- function(exposure_dat, outcome_id){
    out <- TwoSampleMR::extract_outcome_data(
      snps = exposure_dat$SNP,
      outcome = outcome_id,
      proxies = TRUE,
      rsq = 0.8,
      maf_threshold = 0.3)

    return(out)
  }

  mr_res_all <-data.frame()

  for (outcome_id in outcomes_to_try){
    outcome_dat <- extract_outcome_data_custom(instruments, outcome_id)
    harmonised_dat <-TwoSampleMR::harmonise_data(exposure_dat = instruments,
                                                 outcome_dat = outcome_dat)
    mr_results <- TwoSampleMR::mr(harmonised_dat,
                                  method_list=c('mr_ivw','mr_egger_regression','mr_weighted_median', 'mr_wald_ratio')) %>%
      TwoSampleMR::split_outcome() %>%
      TwoSampleMR::split_exposure() %>%
      tidyr::separate(outcome, "outcome", sep="[(]") %>%
      TwoSampleMR::generate_odds_ratios()

    mr_res_all<-dplyr::bind_rows(mr_res_all, mr_results)
  }

  mr_res_all_ivw <- mr_res_all %>%
    dplyr::filter(method %in% c("Inverse variance weighted", "Wald ratio")) %>%
    dplyr::rename(outcome.details=id.outcome,
                  or_loci=or_lci95,
                  or_upci=or_uci95) %>%
    dplyr::mutate(OR_CI = paste0(round(or,2), " [",round(or_loci,2) ,":",round(or_upci,2), "]")) %>%
    dplyr::mutate(effect_direction = ifelse(or_loci > 1 & or_upci >= 1, 'positive',
                                            ifelse(or_loci < 1 & or_upci <= 1, 'negative', 'overlaps null')))
  return(mr_res_all_ivw)
}


#' This is some descriptio of this function.
#' @title function to call enrich fro pathways
#'
#' @description function to call enrich fro pathways
#'
#' @details function to call enrich fro pathways
#' @importFrom enrichR enrichr
#' @importFrom dplyr bind_rows filter arrange
#' @importFrom tidyr separate
#' @param gene_list Character vector of gene names or data.
#' @param dbs Character vector of databases to search. See https://maayanlab.cloud/Enrichr/ for available databases
#' @param adjpval_filter =0.05
#' @return a dataframe
#' @export

enrich_dbs<-function(gene_list, dbs, adjpval_filter = 0.05){

  enriched <- enrichR::enrichr(gene_list, dbs)
  # flatten list into a table; handle empty tables
  for (db_name in names(enriched)){

    if( dim(enriched[[db_name]])[1] > 0){

      enriched[[db_name]]$db <- db_name
    } else  {
      # if it's empty, delete it
      enriched[[db_name]] <- NULL
    }
  }

  enriched_df<-dplyr::bind_rows(enriched)

  if (dim(enriched_df)[1] > 0){

    enriched_df<- enriched_df %>%
      dplyr::filter(Adjusted.P.value < adjpval_filter) %>%
      tidyr::separate(Overlap, into = c("found_genes", "total_genes"), sep="/", remove = F)%>%
      dplyr::arrange(Odds.Ratio)

  } else{

    enriched_df<-data.frame()

  }
}

#' This is some descriptio of this function.
#' @title function to call enrich fro pathways
#'
#' @description function to get eQTL by queried snp.
#'
#' @details function to get eQTL by queried snp.
#' @importFrom httr GET content
#' @importFrom jsonlite fromJSON
#' @importFrom tidyr separate
#' @importFrom dplyr select arrange bind_rows
#' @importFrom purrr discard
#' @param variant the SNP list.
#' @return a dataframe
#' @export

get_eQTL_for_snp <- function(variant){

  request = httr::GET(url = "http://www.ebi.ac.uk/eqtl/api/associations",
                      query = list(
                        variant_id = variant,
                        size = 1000,
                        p_upper = 5e-8)
  )
  stopifnot(request$status_code==200)

  response = httr::content(request, as = "text", encoding = "UTF-8")

  variant_assoc = jsonlite::fromJSON(response, flatten = TRUE)$`_embedded`$associations

  if (length(variant_assoc)!=0){

    for (i in 1:length(variant_assoc)) {

      variant_assoc[[i]]<-variant_assoc[[i]] %>%
        purrr::discard(is.null) # drop any null items

    }

    variant_assoc_df <- dplyr::bind_rows(variant_assoc) %>%
      dplyr::select(rsid, gene_id, qtl_group, pvalue, dplyr::everything()) %>%
      dplyr::arrange(pvalue)

  }else{

    return(data.frame())

  }

  nextq = jsonlite::fromJSON(response, flatten = TRUE)$`_links`

  if( any(names(nextq)=="next")){

    stop("API limits 1000 results, but there might be more -- need to investigate")

  }
  return(variant_assoc_df)
}

#' This is some descriptio of this function.
#' @title function to call enrich fro pathways
#'
#' @description function to get eQTL by queried snp.
#'
#' @details function to get eQTL by queried snp.
#' @importFrom biomaRt useEnsembl useDataset getBM
#' @param values Values of the filter, e.g. vector of affy IDs. If multiple filters are specified then the argument should be a list of vectors of which the position of each vector corresponds to the position of the filters in the filters argument.
#' @param GRCh =37 or 38.
#' @return a dataframe
#' @export

map_ensembl_to_genename <- function(values,
                                    GRCh=37){

  #retrieve all genes with their GRCh37 coordinates from biomart

  mart_grch = biomaRt::useEnsembl(biomart="ensembl",GRCh=GRCh)

  mart_grch = biomaRt::useDataset("hsapiens_gene_ensembl", mart_grch)

  # retrieve gene symbols using biomart (eQTL Catalog returns ensembl)

  mart_query = biomaRt::getBM(mart=mart_grch,
                              attributes=c("ensembl_gene_id","hgnc_symbol"),
                              filters= c("ensembl_gene_id"),
                              values=unique(values))

  return(mart_query)
}

#' This is some descriptio of this function.
#' @title function to query mr result based on gived outcome_id in epigraphdb database website.
#'
#' @description function to query mr result based on gived outcome_id in epigraphdb database website.
#'
#' @details function to query mr result based on gived outcome_id in epigraphdb database website.
#' @importFrom TwoSampleMR extract_instruments extract_outcome_data harmonise_data mr generate_odds_ratios split_outcome split_exposure
#' @importFrom tidyr separate
#' @importFrom dplyr bind_rows filter arrange
#' @param outcome_id the mabase id of outcome.
#' @param outcome_trait the outcome trait.
#' @param p1 =1, query all MR results for the outcomes, not restricting by p-value.
#' @param p2 =5e-8, V2 - alternative - more stringent with pval < 5e-8 for all.
#' @param save_path save the result to file path.
#' @return a dataframe
#' @export

mr_eve <- function(outcome_id="ukb-a-100",
                   outcome_trait,
                   p1=1,
                   p2=5e-8,
                   save_path="D:\\课题组的活\\开题\\1"){

  query =
    paste0("
    MATCH (exposure:Gwas)-[mr:MR_EVE_MR]->(outcome:Gwas)
    WHERE outcome.id in ['", paste0(outcome_id="ukb-a-100", collapse = "', '"),"']
    AND  not exposure.id  in ['", paste0(outcome_id="ukb-a-100", collapse = "', '"),"']
    AND (not (toLower(exposure.trait) contains ",paste0("'",as.character(outcome_trait),"'"),")) ",
           "AND mr.pval < ", p1 ,
           " with mr, exposure, outcome
    ORDER BY mr.pval
    RETURN exposure.id, exposure.trait, exposure.sample_size, exposure.sex, exposure.note,
          toInteger(exposure.year) as year, exposure.author as author, exposure.consortium as consortium,
              outcome.id, outcome.sample_size, toInteger(outcome.ncase) as N_case, outcome.year, outcome.nsnp,
              mr.pval, mr.b, mr.se,mr.nsnp,mr.method, mr.moescore
    ")

  full_results <- query_epigraphdb_as_table(query)

  print(dim(full_results))

  paste0("total number of exposure traits connected to outcomes with any result:",length(unique(full_results$exposure.id)))

  if(dim(full_results)[1]==0){

    message(paste0("no found the significant exposure trait for ",outcome_trait))

  }else{

    # calculate CI and get effect direction
    full_results <- full_results %>%
      dplyr::mutate( loci = mr.b - 1.96 * mr.se,
                     upci = mr.b + 1.96 * mr.se,
                     or = exp(mr.b),
                     or_loci = exp(loci),
                     or_upci = exp(upci),
                     OR_CI = paste0(round(or,2), " [",round(or_loci,2) ,":",round(or_upci,2), "]")) %>%
      dplyr::mutate(effect_direction = ifelse(or_loci > 1 & or_upci >= 1, 'positive',
                                              ifelse(or_loci < 1 & or_upci <= 1, 'negative', 'overlaps null'))) %>%
      dplyr::mutate(`MR method and score` = paste0(mr.method," / ", mr.moescore))


    print(full_results)

    # save all for supl data
    full_results_save <- full_results

    readr::write_csv(x = full_results_save,file = paste0(save_path,"/","01_all_mreve_",outcome_id,"_",outcome_trait,"_results.csv"))


    sub_results <- full_results %>%
      dplyr::filter(effect_direction != 'overlaps null')

    print(paste0("unique traits with non-null effect: ",length(unique(sub_results$exposure.id))))

    # now re-extract the full MR results (for all outcomes) for those 1970 traits
    query = paste0("
      MATCH (exposure:Gwas)-[mr:MR_EVE_MR]->(outcome:Gwas)
      WHERE outcome.id in ['", paste0(outcome_id="ukb-a-100", collapse = "', '"),"']
      AND  exposure.id  in ['", paste0(sub_results$exposure.id, collapse = "', '"),"']
      with mr, exposure, outcome
      ORDER BY exposure.trait
      RETURN exposure.id, exposure.trait, exposure.sample_size, exposure.sex, exposure.note,
      toInteger(exposure.year) as year, exposure.author as author, exposure.consortium as consortium,
              outcome.id, outcome.sample_size, toInteger(outcome.ncase) as N_case, outcome.year, outcome.nsnp,
              mr.pval, mr.b, mr.se, mr.nsnp, mr.method, mr.moescore
      ")

    out3 <- query_epigraphdb_as_table(query)

    print(dim(out3))

    print(length(unique(out3$exposure.id)))

    readr::write_tsv(x = out3, file = paste0(save_path,"/02_",outcome_id,"_",outcome_trait,"_all_mr_fromCIs.tsv"))

    ## review CIs vs pval

    full_results %>%
      dplyr::filter(mr.pval < 0.05) %>%
      dplyr::count(effect_direction)

    full_results %>%
      dplyr::filter(effect_direction != 'overlaps null') %>%
      dplyr::count(effect_direction)

    full_results %>%
      dplyr::filter(mr.pval >= 0.05) %>%
      dplyr::count(effect_direction)

    full_results %>%
      dplyr::filter(effect_direction == 'overlaps null') %>%
      dplyr::filter(mr.pval < 0.05) %>%
      dim()

    test <- full_results %>%
      dplyr::filter(mr.pval >= 0.05) %>%
      dplyr::filter(effect_direction != 'overlaps null') %>%
      dplyr::select(exposure.id, exposure.trait,outcome.id,mr.method, mr.pval, or, or_loci, or_upci, OR_CI, effect_direction, mr.b, mr.se,loci, upci  ) %>%
      dplyr::distinct()

    readr::write_csv(x = test, file = paste0(save_path,"/03_",outcome_id,"_",outcome_trait,"_mismatch_pval_ci.csv"))


    # V2 - alternative - more stringent with pval < 5e-8 for all

    # get traits that have pval < 5e-8

    query =
      paste0("
    MATCH (exposure:Gwas)-[mr:MR_EVE_MR]->(outcome:Gwas)
    WHERE outcome.id in ['", paste0(outcome_id="ukb-a-100", collapse = "', '"),"']
    AND  not exposure.id  in ['", paste0(outcome_id="ukb-a-100", collapse = "', '"),"']
    AND (not (toLower(exposure.trait) contains ",paste0("'",as.character(outcome_trait),"'"),")) ",
             "AND mr.pval < ", p2 ,
             " with mr, exposure, outcome
    ORDER BY mr.pval
    RETURN exposure.id, exposure.trait, exposure.sample_size,
            collect(outcome.id) as outcome_ids,
            collect(mr.pval) as MR_pvals, collect(mr.b) as MR_beta
    ")

    out<-query_epigraphdb_as_table(query)

    print(dim(out))

    # get MR results for those traits with outcome trait

    query = paste0("
      MATCH (exposure:Gwas)-[mr:MR_EVE_MR]->(outcome:Gwas)
      WHERE outcome.id in ['", paste0(outcome_id="ukb-a-100", collapse = "', '"),"']
      AND  exposure.id  in ['", paste0(out$exposure.id, collapse = "', '"),"']
      with mr, exposure, outcome
      ORDER BY exposure.trait
      RETURN exposure.id, exposure.trait, exposure.sample_size, exposure.sex, exposure.note,
      toInteger(exposure.year) as year, exposure.author as author, exposure.consortium as consortium,
              outcome.id, outcome.sample_size, toInteger(outcome.ncase) as N_case, outcome.year, outcome.nsnp,
              mr.pval, mr.b, mr.se, mr.method, mr.moescore
      ")

    out2<-query_epigraphdb_as_table(query)

    print(dim(out2))

    print(length(unique(out2$exposure.id)))

    readr::write_tsv(x = out2, file = paste0(save_path,"/04_",outcome_id,"_",outcome_trait,"_all_mr_madewR.tsv"))

    return(out2)

  }

}

#' @details functions output the data with mr-eve confounders, mediators, reverse_mediator or collider query.
#' @importFrom dplyr %>% mutate select filter starts_with
#' @param df the data queried by epigraphdb.
#' @param type ="confounder", "mediator", "reverse_mediator", "collider"
#' @return a dataframe
#' @export

tidy_conf_query_output <- function(df,
                                   type){

  if (dim(df)[1] != 0){

    df_OR <- df %>%
      dplyr::mutate( r1.loci = r1.b - 1.96 * r1.se,
                     r1.upci = r1.b + 1.96 * r1.se,
                     r1.or = exp(r1.b),
                     r1.or_loci = exp(r1.loci),
                     r1.or_upci = exp(r1.upci),
                     r1.OR_CI = paste0(round(r1.or,2), " [",round(r1.or_loci,2) ,":",round(r1.or_upci,2), "]")) %>%
      dplyr::mutate(r1.effect_direction = ifelse(r1.or_loci > 1 & r1.or_upci >= 1, 'positive',
                                                 ifelse(r1.or_loci < 1 & r1.or_upci <= 1, 'negative', 'overlaps null'))) %>%
      dplyr::mutate( r2.loci = r2.b - 1.96 * r2.se,
                     r2.upci = r2.b + 1.96 * r2.se,
                     r2.or = exp(r2.b),
                     r2.or_loci = exp(r2.loci),
                     r2.or_upci = exp(r2.upci),
                     r2.OR_CI = paste0(round(r2.or,2), " [",round(r2.or_loci,2) ,":",round(r2.or_upci,2), "]")) %>%
      dplyr::mutate(r2.effect_direction = ifelse(r2.or_loci > 1 & r2.or_upci >= 1, 'positive',
                                                 ifelse(r2.or_loci < 1 & r2.or_upci <= 1, 'negative', 'overlaps null'))) %>%

      dplyr::mutate( r3.loci = r3.b - 1.96 * r3.se,
                     r3.upci = r3.b + 1.96 * r3.se,
                     r3.or = exp(r3.b),
                     r3.or_loci = exp(r3.loci),
                     r3.or_upci = exp(r3.upci),
                     r3.OR_CI = paste0(round(r3.or,2), " [",round(r3.or_loci,2) ,":",round(r3.or_upci,2), "]")) %>%
      dplyr::mutate(r3.effect_direction = ifelse(r3.or_loci > 1 & r3.or_upci >= 1, 'positive',
                                                 ifelse(r3.or_loci < 1 & r3.or_upci <= 1, 'negative', 'overlaps null'))) %>%

      dplyr::select(-contains("loci"), -contains("upci")) %>%

      dplyr::select("exposure.trait", "exposure.id", "outcome.trait", "outcome.id", "med.trait","med.id" ,
                    dplyr::starts_with("r1"),  dplyr::starts_with("r2"),dplyr::starts_with("r3")) %>%
      dplyr::filter(r1.effect_direction != 'overlaps null' &r2.effect_direction != 'overlaps null' & r3.effect_direction != 'overlaps null') %>%
      dplyr::mutate(type = type)

  } else {

    df_OR <- data.frame()

  }

}

#' @details functions output the data with mr-eve confounders, mediators, reverse_mediator or collider query.
#' @importFrom dplyr bind_rows
#' @param exposure the mrbase id of exposure.
#' @param outcome the mrbase id of outcome.
#' @param outcome_trait the outcome name.
#' @param pval_threshold =0.05.
#' @param pval_threshold_outcome =1.
#' @param mediator = F.
#' @param confounder = F.
#' @param collider = F.
#' @param reverse_mediator = F.
#' @return a dataframe
#' @export

query_and_tidy_conf <- function(exposure,
                                outcome,
                                outcome_trait,
                                pval_threshold=0.05,
                                pval_threshold_outcome = 1,
                                mediator= F,
                                confounder= F,
                                collider= F,
                                reverse_mediator = F){

  res_list_tidy <- list()

  if (mediator){

    print("Querying MR-EvE for mediators ...")
    mediator_query = paste0("
    MATCH (med:Gwas)<-[r1:MR_EVE_MR]- (exposure:Gwas) -[r2:MR_EVE_MR]->(outcome:Gwas) <-[r3:MR_EVE_MR]-(med:Gwas)
    WHERE exposure.id = '", exposure, "' AND outcome.id = '", outcome,"' AND (not (toLower(exposure.trait) contains ",paste0("'",as.character(outcome_trait),"'"),")) ",
                            " AND r1.pval < ", pval_threshold, " AND r3.pval < ", pval_threshold_outcome, "
    AND med.id <> exposure.id AND med.id <> outcome.id AND exposure.id <> outcome.id AND med.trait <> exposure.trait AND med.trait <> outcome.trait AND exposure.trait <> outcome.trait
    RETURN exposure {.id, .trait}, outcome {.id, .trait}, med {.id, .trait}, r1 {.b, .se, .pval, .selection, .method, .moescore}, r2 {.b, .se, .pval, .selection, .method, .moescore}, r3 {.b, .se, .pval, .selection, .method, .moescore} ORDER BY r1.p
          ")

    mediators = query_epigraphdb_as_table(mediator_query)

    mediators =  tidy_conf_query_output(mediators, type = "mediator")

    res_list_tidy$mediaotrs <- mediators
  }

  if (confounder){
    print("Querying MR-EvE for confounders ...")
    confounder_query = paste0("
    MATCH (med:Gwas)-[r1:MR_EVE_MR]-> (exposure:Gwas) -[r2:MR_EVE_MR]->(outcome:Gwas) <-[r3:MR_EVE_MR]-(med:Gwas)
    WHERE exposure.id = '", exposure, "' AND outcome.id = '", outcome,"' AND (not (toLower(exposure.trait) contains ",paste0("'",as.character(outcome_trait),"'"),")) ",
                              " AND r1.pval < ", pval_threshold, " AND r3.pval < ", pval_threshold_outcome, "
    AND med.id <> exposure.id AND med.id <> outcome.id AND exposure.id <> outcome.id AND med.trait <> exposure.trait AND med.trait <> outcome.trait AND exposure.trait <> outcome.trait
    RETURN exposure {.id, .trait}, outcome {.id, .trait}, med {.id, .trait}, r1 {.b, .se, .pval, .selection, .method, .moescore}, r2 {.b, .se, .pval, .selection, .method, .moescore}, r3 {.b, .se, .pval, .selection, .method, .moescore} ORDER BY r1.p
          ")

    confounders = query_epigraphdb_as_table(confounder_query)

    confounders = tidy_conf_query_output(confounders, type = "confounder")

    res_list_tidy$confounders <- confounders
  }

  if (collider){

    print("Querying MR-EvE for colliders ...")
    collider_query = paste0("
    MATCH (med:Gwas)<-[r1:MR_EVE_MR]- (exposure:Gwas) -[r2:MR_EVE_MR]->(outcome:Gwas) -[r3:MR_EVE_MR]->(med:Gwas)
    WHERE exposure.id = '", exposure, "' AND outcome.id = '", outcome,"' AND (not (toLower(exposure.trait) contains ",paste0("'",as.character(outcome_trait),"'"),")) ",
                            " AND r1.pval < ", pval_threshold, " AND r3.pval < ", pval_threshold_outcome, "
    AND med.id <> exposure.id AND med.id <> outcome.id AND exposure.id <> outcome.id AND med.trait <> exposure.trait AND med.trait <> outcome.trait AND exposure.trait <> outcome.trait
    RETURN exposure {.id, .trait}, outcome {.id, .trait}, med {.id, .trait}, r1 {.b, .se, .pval, .selection, .method, .moescore}, r2 {.b, .se, .pval, .selection, .method, .moescore}, r3 {.b, .se, .pval, .selection, .method, .moescore} ORDER BY r1.p
          ")

    colliders = query_epigraphdb_as_table(collider_query)

    colliders = tidy_conf_query_output(colliders, type = "collider")

    res_list_tidy$colliders <- colliders

  }

  if (reverse_mediator){

    print("Querying MR-EvE for reverse mediators ...")

    reverse_mediator_query = paste0("
    MATCH (med:Gwas)-[r1:MR_EVE_MR]-> (exposure:Gwas) -[r2:MR_EVE_MR]->(outcome:Gwas) -[r3:MR_EVE_MR]->(med:Gwas)
    WHERE exposure.id = '", exposure, "' AND outcome.id = '", outcome,"' AND (not (toLower(exposure.trait) contains ",paste0("'",as.character(outcome_trait),"'"),")) ",
                                    " AND r1.pval < ", pval_threshold, " AND r3.pval < ", pval_threshold_outcome, "
    AND med.id <> exposure.id AND med.id <> outcome.id AND exposure.id <> outcome.id AND med.trait <> exposure.trait AND med.trait <> outcome.trait AND exposure.trait <> outcome.trait
    RETURN exposure {.id, .trait}, outcome {.id, .trait}, med {.id, .trait}, r1 {.b, .se, .pval, .selection, .method, .moescore}, r2 {.b, .se, .pval, .selection, .method, .moescore}, r3 {.b, .se, .pval, .selection, .method, .moescore} ORDER BY r1.p
          ")

    rev_mediators = query_epigraphdb_as_table(reverse_mediator_query)

    rev_mediators = tidy_conf_query_output(rev_mediators, type = "reverse_mediator")

    res_list_tidy$rev_mediators <- rev_mediators

  }


  res_list_tidy_df <- dplyr::bind_rows(res_list_tidy)

  return(res_list_tidy_df)

}


#' @details convert MR results to OR.
#' @importFrom dplyr %>% mutate
#' @param dat the mr result.
#' @return a dataframe
#' @export

tidy_display_numbers<- function(dat){

  dat<- dat %>%
    dplyr::mutate(loci = mr.b - 1.96 * mr.se,
                  upci = mr.b + 1.96 * mr.se,
                  or = exp(mr.b),
                  or_loci = exp(loci),
                  or_upci = exp(upci),
                  OR_CI = paste0(round(or,2), " [",round(or_loci,2) ,":",round(or_upci,2), "]")) %>%
    dplyr::mutate(effect_direction = ifelse(or_loci > 1 & or_upci >= 1, 'positive',
                                            ifelse(or_loci < 1 & or_upci <= 1, 'negative', 'overlaps null'))) %>%
    # fix issue with rounding negative effect to 1
    dplyr::mutate(OR_CI = ifelse(effect_direction == 'negative' & OR_CI == '1 [1:1]', "0.99 [0.99:0.99]", OR_CI)) %>%
    # if beta is really really small, let it count as if it overlaps the null
    dplyr::mutate(effect_direction = ifelse(grepl("0 [0:",OR_CI, fixed = T), 'overlaps null' , effect_direction)) %>%
    # convert all super small pval into 1e-15
    dplyr::mutate(pval_truncated = ifelse(mr.pval < 1e-15, 1e-15, mr.pval),
                  #log transform pvals
                  log10pval_trunc = as.integer(-log10(pval_truncated))) %>%
    # round beta for display
    dplyr::mutate(mr.b = round(mr.b, 3)) %>%
    # tidy moe display
    dplyr::mutate(`MR method and score` = paste0(mr.method," / ", mr.moescore)) %>%
    # emply col fro display
    dplyr::mutate(empty_col = " ")

  return(dat)
}


#' @details create categories of exposure traits.
#' @importFrom dplyr %>% mutate case_when
#' @importFrom stringr str_replace
#' @param dat the mr result.
#' @return a dataframe
#' @export

create_exposure_categories <- function(dat){

  # shorten names for some:
  dat <- dat %>%
    dplyr::mutate(exposure = stringr::str_replace(exposure.trait, "Treatment/medication code", "Drug"),
                  exposure = stringr::str_replace(exposure, "Non-cancer illness code  self-reported", "Illness"),
                  exposure = stringr::str_replace(exposure, "Diagnoses - main ICD10:", "Diagnosis"),
                  exposure = stringr::str_replace(exposure, "Mineral and other dietary supplements:", "Supplements"),
                  exposure = stringr::str_replace(exposure, "Vitamin and mineral supplements:", "Vitamins"),
                  exposure = stringr::str_replace(exposure, "Illness  injury  bereavement  stress in last 2 years", "Stress"),
                  exposure = stringr::str_replace(exposure, "Types of transport used (excluding work)", "Transport"),
                  exposure = stringr::str_replace(exposure, "Never eat eggs, dairy, wheat, sugar:", "Never eat:"),
                  exposure = stringr::str_replace(exposure, "Types of physical activity in last 4 weeks", "Physical activity")) %>%
    dplyr::mutate(exposure = gsub("Comparative ", "", exposure)) %>%
    dplyr::mutate(exposure = gsub(" (last menstrual period)", "", exposure, fixed=T)) %>%
    dplyr::mutate(exposure = gsub("hormone-replacement therapy (HRT)", "HRT", exposure, fixed=T)) %>%
    dplyr::mutate(exposure = gsub("\\(eg.*)", "", exposure)) %>%
    dplyr::mutate(exposure = gsub("Other (e.g. ", "(", exposure, fixed = T)) %>%
    dplyr::mutate(exposure = gsub(", because of other reasons", "", exposure)) %>%
    dplyr::mutate(exposure = gsub("or pain relief  constipation  heartburn", "", exposure)) %>%
    dplyr::mutate(exposure = gsub("or pain relief, constipation, heartburn", "", exposure)) %>%
    dplyr::mutate(exposure = gsub("\\(excluding work)", "", exposure)) %>%
    dplyr::mutate(exposure = gsub("\\(not as a means of transport)", "", exposure)) %>%
    dplyr::mutate(exposure = gsub("Number of days/week of", "Days/week of", exposure)) %>%
    dplyr::mutate(exposure = gsub(" eggs, dairy, wheat, sugar", ": ", exposure)) %>%
    dplyr::mutate(exposure = gsub("Number of cigarettes currently smoked daily", " # cigarettes daily ", exposure)) %>%
    dplyr::mutate(exposure = gsub("Pack years adult smoking as proportion of life span exposed to smoking", "Pack years/life span exposed to smoking ", exposure)) %>%
    dplyr::mutate(exposure = dplyr::case_when(grepl("presbyopia", exposure) ~ "Eyesight problems1",
                                              grepl("hypermetropia", exposure) ~ "Eyesight problems2",
                                              TRUE ~ exposure)) %>%



    ### grouping stuff into categories
    dplyr::mutate(exposure_cat = dplyr::case_when(
      grepl("prot", exposure.id) ~ 'Proteins',
      grepl("met", exposure.id) ~ 'Metabolites',
      grepl("Drug|Medication|prescription medications", exposure, ignore.case = T) ~ "Drugs",
      grepl("Diagnos|Illness|cancer|neoplasm|glioma|diagnos|death|disorder|eye|carcino|colitis|disease|diabetes|asthma|sclerosis|infarction|neuroblastoma|arthritis|eczema|cholangitis|pain|hernia|Polyarthropathies|Malignant|Siatica|Gonarthrosis", exposure, ignore.case = T) ~ "Diseases",
      grepl("blood pressure|heart|Pulse|Cardiac|thromboembolism", exposure, ignore.case = T) ~ "CVD_related",
      grepl("operation|operative|Methods of admission|Number of treatments|Spells in hospital|hospital episode", exposure, ignore.case = T) ~ "Medical Procedures",
      grepl("cylindrical|meridian|asymmetry|glasses|hearing|Corneal|ocular|logMAR|teeth|dental|Spherical power", exposure,  ignore.case = T) ~ "eye_hearing_dental",
      grepl("waist|hip c|hip r|obesity|trunk|mass|weight|bmi|body size|height|impedance|fat percentage|body fat|Basal metabolic rate", exposure, ignore.case = T) ~ "Antrophometric",
      grepl("FEV1|FVC|Wheeze|cough", exposure,  ignore.case = T) ~ "lung_related",
      grepl("count|volume|percentage|reticulocyte|Platelet", exposure,  ignore.case = T) ~ "cells_related",
      grepl("alco|wine|spirits|beer", exposure, ignore.case = T) ~ "Alcohol", # must be before diet
      grepl("vitamin|suppl", exposure, ignore.case = T) ~ "Diet and supplements",
      grepl("intake|diet|food|milk|dairy|coffee|cereal|butter|bread|Never eat", exposure, ignore.case = T) ~ "Diet and supplements",
      grepl("age at|age started|parous|contraceptive pill|replacement therapy|HRT|menopause|menarche|live birth|oophorectomy|hysterectomy|menstrual|sexual", exposure, ignore.case = T) ~ "Reproductive",
      grepl("smok|cigar", exposure, ignore.case = T) ~ "Smoking",
      grepl("activi|transport |diy|walking|walked|Time spent|Weekly usage of|stair climbing|walk|spend outdoors", exposure, ignore.case = T) ~ "Physical activity",
      grepl("sleep|Snoring|chronotype|Getting up in morning|Nap during day", exposure, ignore.case = T) ~ "Sleep",
      grepl("Iron|Testosterone|Urate|Urea|Glucose|Sodium", exposure, ignore.case = T) ~ 'Other biomarkers',
      grepl("LDL|HDL|VLDL|cholest|trigl|cholesterol|glyceride|total lipids|Serum total", exposure.trait, ignore.case = T) ~ 'Lipids',
      grepl("Albumin|Apoliprotein|Adiponectin|Lipoprotein|reactive protein|Creatinine|Ferritin|Transferrin|transferase|Haemoglobin|cystatin|SHBG|bilirubin|Total protein|phosphatase|IGF", exposure, ignore.case = T) ~ 'Proteins',
      grepl("Qualifications|GCSE|Townsend|schooling|College|intelligence|arithmetic|education", exposure, ignore.case = T) ~ 'Education',
      grepl("anxiety|feelings|embarrassment|worr|Bulimia|depressed|guilty|Miserableness|mood|Neuroticism|unenthusiasm|tenseness|Loneliness|self-harm|Risk taking|highly strung|ADHD|Drive faster|nerves|irritability|satisfaction", exposure, ignore.case = T) ~ 'psychiatric_or_mental_health',

      TRUE ~ 'other')) %>%
    dplyr::mutate(exposure_cat = ifelse(grepl("LDL|HDL|VLDL|cholest|trigl|cholesterol|glyceride|total lipids|Serum total", exposure.trait, ignore.case = T) & exposure_cat %in% c('Metabolites', 'Other biomarkers'), "Lipids", exposure_cat)) %>%  # recapture those in met-a
    dplyr::mutate(exposure_cat = ifelse(grepl("Average number|Ratio of", exposure.trait, ignore.case = T) , "other", exposure_cat)) %>%  # recapture those in met-a
    dplyr::mutate(exposure_cat = ifelse(grepl("albumin", exposure.trait, ignore.case = T) , "Proteins", exposure_cat))
  return(dat)

}

#' @details create categories of exposure traits.
#' @importFrom dplyr %>% mutate case_when coalesce
#' @importFrom tidyr separate
#' @param dat the mr result.
#' @return a dataframe
#' @export
add_exposure_labels <-function(dat){

  dat<- dat %>%
    # create tidy exposure sample size
    dplyr::mutate(exposure.ss = format(exposure.sample_size, big.mark = ",")) %>%
    tidyr::separate(exposure.ss, into=c('exposure.ss', 'tmp'), sep=",") %>%
    dplyr::mutate(exposure.ss_label = paste0(" [",exposure.ss,"K]"),
                  exposure.ss_label = gsub("   NA", "NA ", exposure.ss_label)) %>%
    dplyr::mutate(exposure.ss_label = ifelse(exposure.sample_size < 1000, paste0(" [",exposure.sample_size,"]"), exposure.ss_label)) %>%

    # add other exposure details
    dplyr::mutate(exposure = dplyr::case_when(
      exposure.sex == 'Females' ~ paste0(exposure, " (F) ", dplyr::coalesce(consortium, author),"/" , year, exposure.ss_label),
      TRUE ~ paste0(exposure, " (M/F) ", dplyr::coalesce(consortium, author),"/" , year, exposure.ss_label))) %>%
    # add note
    dplyr::mutate(exposure = dplyr::case_when(grepl('Adjusted for BMI', exposure.note) ~ paste0(exposure, " AdjBMI"),
                                              TRUE ~ exposure)) %>%
    # create ukb data tag for filtering in the app
    dplyr::mutate(ukb_tag = ifelse(author %in% c('Neale lab', 'Neale'), "Neale lab",
                                   ifelse (consortium == "MRC-IEU", "MRC-IEU", NA)))

  return(dat)
}

#' @details create categories of exposure traits.
#' @importFrom dplyr %>%
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot aes geom_point scale_size facet_grid scale_colour_manual scale_y_discrete theme_light labs theme element_text
#' @param input the mr result.
#' @param font_size =10
#' @return a dataframe
#' @export
plot_bubble_plot <- function(input, font_size = 10){

  input <- input %>%
    # add effect brackets
    create_beta_ranges()

  # create names palette
  getPalette<-grDevices::colorRampPalette((RColorBrewer::brewer.pal(n=7, name = "PiYG")))
  cols<- getPalette(13)
  names(cols) <- levels(input$mr.b_col_verbose)

  p<-ggplot2::ggplot(input, ggplot2::aes(x= outcome, y =stats::reorder(exposure, mr.b), color= mr.b_col_verbose,  size = log10pval_trunc,
                                         text = paste(" ", empty_col,
                                                      '</br>Exposure ID: ', exposure.id,
                                                      '</br>Exposure: ', exposure.trait,
                                                      '</br>Outcome ID: ',  outcome.id,
                                                      '</br>Outcome: ',  outcome,
                                                      '</br>Beta: ', mr.b,
                                                      '</br>P-value: ', mr.pval,
                                                      '</br>P-value (-log10): ', log10pval_trunc,
                                                      '</br>Odds ratio: ', OR_CI,
                                                      '</br>MR method and score: ', `MR method and score`)

  )) +
    ggplot2::geom_point()+
    ggplot2::scale_size(breaks = c(0,3,8,15))+

    ggplot2::scale_colour_manual(values = cols)+

    ggplot2::scale_y_discrete(position = "left")+
    ggplot2::facet_grid(exposure_cat~outcome, scales = 'free') +
    ggplot2::theme_light()+
    #theme_solarized()+
    ggplot2::labs(size="-log10(pval)",
                  color='beta (effect size)', y="Exposures")+
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35, hjust=1),
                   axis.text = ggplot2::element_text(size = font_size),
                   legend.position = 'right')



  return(p)
}

#' @details add effect brackets.
#' @importFrom dplyr %>% mutate case_when
#' @param input the mr result.
#' @return a dataframe
#' @export
create_beta_ranges <- function(input){

  # add effect brackets
  input <- input %>%
    dplyr::mutate(mr.b_col = dplyr::case_when( mr.b > 1 ~ 2,
                                               mr.b > 0.5 ~ 1,
                                               mr.b > 0.25 ~ 0.5,
                                               mr.b > 0.10 ~ 0.25,
                                               mr.b > 0.01 ~ 0.1,
                                               mr.b > 0.001 ~ 0.01,
                                               mr.b < -1 ~ -2,
                                               mr.b < -0.5 ~ -1,
                                               mr.b < -0.25 ~ -0.5,
                                               mr.b < -0.10 ~ -0.25,
                                               mr.b < -0.01 ~ -0.1,
                                               mr.b < -0.001 ~ -0.01,
                                               TRUE ~ 0)) %>%
    dplyr::mutate(mr.b_col_verbose = dplyr::case_when(mr.b > 1 ~    "> 1",
                                                      mr.b > 0.5 ~  "0.50:1",
                                                      mr.b > 0.25 ~ '0.25:0.50',
                                                      mr.b > 0.10 ~ '0.10:0.25',
                                                      mr.b > 0.01 ~ '0.001:0.10',
                                                      mr.b > 0.001 ~'<  0.001',
                                                      mr.b < -1 ~   '< -1',
                                                      mr.b < -0.5 ~  "-0.50:-1",
                                                      mr.b < -0.25 ~ '-0.25:-0.50',
                                                      mr.b < -0.10 ~ '-0.10:-0.25',
                                                      mr.b < -0.01 ~ '-0.001:-0.10',
                                                      mr.b < -0.001 ~'> -0.001',
                                                      TRUE ~ '0')) %>%
    dplyr::mutate(mr.b_col_verbose = factor(mr.b_col_verbose, levels=c("> 1",
                                                                       "0.50:1",
                                                                       '0.25:0.50',
                                                                       '0.10:0.25',
                                                                       '0.001:0.10',
                                                                       '<  0.001',
                                                                       '0',
                                                                       '> -0.001',
                                                                       '-0.001:-0.10',
                                                                       '-0.10:-0.25',
                                                                       '-0.25:-0.50',
                                                                       "-0.50:-1",
                                                                       '< -1')))
  return(input)
}
