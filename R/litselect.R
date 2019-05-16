# palmer litselect api

palmer_api <- function(api_type, terms=NA, terms2=NA, 
                       bin_cut=0.1, go_ratio=1.0, gene_ratio=2.0,
                       no_multi=FALSE){



  # check dependencies
  if(!require("httr", character.only=T, quietly=T)){
    warning("httr package is required to use the API!")
    install.packages("httr", character.only=T)
  }

  if(!require("jsonlite", character.only=T, quietly=T)){ 
    warning("jsonlite package is required to use the API!")
    install.packages("jsonlite", character.only=T)
  }
  
  library("httr", quietly=T)
  library("jsonlite", quietly=T)

  multi_filter <- function(multiple_matches){
    # filter the dataset
    multi <- multiple_matches[multiple_matches$multiple_matches==TRUE, ]
    non_multi <- multiple_matches[multiple_matches$multiple_matches==FALSE, ]
    # if no multiple matches, return
    if(nrow(multi)==0){
      return
    }
    # get a list of unique input terms that have multiple matches
    unique_genes <- unique(multi$input_term)
    # create an empty dataframe
    final_df <- multi[FALSE, ]
    for(i in 1:length(unique_genes)){
      subset <- multi[multi$input_term==tolower(unique_genes[i]),]
      subset_row <- nrow(subset)
      # set the rownames for indexing
      rownames(subset) <- 1:subset_row
      select <- -1
      # to do: add judge if select is an integer here
      # to do: add an option to pass/ignore/choose none
      # check until get the validated input
      while(!(select %in% seq(0, subset_row, 1))){
        print(subset)
        cat("\nPlease choose the genes you want to keep (0 to keep them all):#")
        select <- readLines(file("stdin"), n=1)
        # to do: validate select/ not necessary turn it into an integer
        select <- as.integer(select)
      }
      if(select <= subset_row&&select >= 1){
        new_add <- subset[select, ]
        final_df[nrow(final_df) + 1, ] <- unname(unlist(new_add))      
      }else if(select == 0){
        print(subset)
        final_df <- rbind(final_df, subset)
      }
    }
    all_df <- rbind(non_multi, final_df)
    rownames(all_df) <- 1:nrow(all_df)
    return(all_df)
  }

  # base url
  base_url <- "http://chunglab.io/GAIL/api"
  # query_body 
  query_body <- NA
  # final results
  final_results <- NA
  # check input
  if(is.na(terms)||(api_type == "matrix"&is.na(terms2))){
    stop("Empty input")
  }
  # api function
  if(api_type == "go"){

    #path
    path <- "/select_GO_terms"
    query_url <- paste(base_url, path, sep="")    
    #to do: check terms
    query_body <- list(terms=terms, ratio=go_ratio)
    html_result <- POST(query_url, body=query_body, encode="json")$content
    html_content <- rawToChar(html_result)
    #check 404 error
    if(grepl("Error code: 400", html_content)){
      stop("Error in conncting the API")
    }else if(grepl("Error code: 500", html_content)){
      stop("Server side error")
    }
    text_result <- fromJSON(html_content)
    if(is.null(text_result$output_gene)){
      warning("No results are found")
    }
    #print(text_result)    
    final_result <- text_result$output_gene   

  }else if(api_type == "gene"){

    #path
    path <- "/select_candidate_genes"
    query_url <- paste(base_url, path, sep="")    
    #to do: check terms
    query_body <- list(terms=terms, ratio=gene_ratio)
    html_result <- POST(query_url, body=query_body, encode="json")$content
    html_content <- rawToChar(html_result)
    #check 404 error
    if(grepl("Error code: 400", html_content)){
      stop("Error in conncting the API")
    }else if(grepl("Error code: 500", html_content)){
      stop("Server side error")
    }
    text_result <- fromJSON(html_content)
    if(is.null(text_result$output_gene)){
      warning("No results are found")
    }
    #print(text_result)    
    final_result <- text_result$output_gene  

  }else if(api_type == "matrix"){
    
    #check input
    #path
    path <- "/select_gene_GO_matrix"
    query_url <- paste(base_url, path, sep="")    
    #to do: check terms
    query_body <- list(gene_list_1=terms, gene_list_2=terms2, binary_cutoff=bin_cut, 
                       GO_ratio=go_ratio, gene_ratio=gene_ratio)
    html_result <- POST(query_url, body=query_body, encode="json")$content
    html_content <- rawToChar(html_result)
    #Check 404 error
    if(grepl("Error code: 400", html_content)){
      stop("Error in conncting the API")
    }else if(grepl("Error code: 500", html_content)){
      stop("Server side error")
    }
    text_result <- fromJSON(html_content)
    if(is.null(text_result$matrix)){
      warning("No results are found")
    }
    #print(text_result)    
    final_result <- text_result  

  }else if(api_type == "id_mapper"){

    #path
    path <- "/id_mapper_1"
    query_url <- paste(base_url, path, sep="")
    html_result <- POST(query_url, body=query_body, encode="json")$content
    html_content <- rawToChar(html_result)  
    text_result <- fromJSON(html_content) 
    if(is.null(text_result$content)){
      warning("No results are found")
    }
    gene_df <- text_result$output_gene
    # to do: response error/warning error handling
    if(no_multi){
      final_result <- gene_df[gene_df$multiple_matches==FALSE, ]
    }else{
      if(!require("getPass", character.only=T, quietly=T)){
        warning("getPass package is required to use the API!")
        install.packages("getPass", character.only=T)
      }      
      final_result <- multi_filter(gene_df)
    }

  }else{

    # if other API parameters, pop up the error  
    stop("API not found!")
  
  }
  return(final_result)

}

