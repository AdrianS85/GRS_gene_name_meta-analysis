# library(reutils)
# normalizePath(" ~ ")
# 
# #
# # combine esearch and efetch
# #
# # Download PubMed records that are indexed in MeSH for both 'Chlamydia' and
# # 'genome' and were published in 2013.
# query <- "Chlamydia[mesh] and genome[mesh] and 2013[pdat]"
# # Upload the PMIDs for this search to the History server
# pmids <- esearch(term = query, "pubmed", usehistory = F)
# pmids
# ## Not run:
# # Fetch the records
# articles <- efetch(pmids)
# # Use XPath expressions with the #xmlValue() or #xmlAttr() methods to directly
# # extract specific data from the XML records stored in the 'efetch' object.
# titles <- articles$xmlValue("//ArticleTitle")
# abstracts <- articles$xmlValue("//AbstractText")
# #
# # combine epost with esummary/efetch
# #
# # Download protein records corresponding to a list of GI numbers.
# uid <- c("194680922", "50978626", "28558982", "9507199", "6678417")
# # post the GI numbers to the Entrez history server
# p <- epost(uid, "protein")
# # retrieve docsums with esummary
# docsum <- content(esummary(p, version = "1.0"), "parsed")
# docsum
# # download FASTAs as 'text' with efetch
# prot <- efetch(p, retmode = "text", rettype = "fasta")
# prot
# # retrieve the content from the efetch object
# fasta <- content(prot)
# 
# 
# 
# library(RISmed)
# 
# test <- EUtilsQuery("myeloma[ti] jones[au]")
# 
# 
# res <- EUtilsSummary("myeloma[ti]",retmax=2,reldate=365)
# summary(res)
# fetch <- EUtilsGet(res)
# 
# leftovers <- readr::read_tsv('data_v4_leftovers_for_ncbi.txt')
# 
# leftovers_vector <- as.vector(leftovers$Probe_ID)
# 
# leftovers_vector_test <- leftovers_vector[1:50]
# 
# fetch <- EUtilsGet(x = leftovers_vector_test, type="epost", db="gene")
# 
# 
# db = 'protein';
# $id_list = '194680922,50978626,28558982,9507199,6678417'
# 
# class()

# The Entrez Programming Utilities (E-utilities) are a set of nine server-side programs that provide a stable interface into the Entrez query and database system at the National Center for Biotechnology Information (NCBI).
# The value of tool should be a string with no internal spaces that uniquely identifies the software producing the request. The value of email should be a complete and valid e-mail address of the software developer and not that of a third-party end user. The value of email will be used only to contact developers if NCBI observes requests that violate our policies, and we will attempt such contact prior to blocking access. 
# In order not to overload the E-utility servers, NCBI recommends that users post no more than three URL requests per second and limit large jobs to either weekends or between 9:00 PM and 5:00 AM Eastern time during weekdays.
# E utilities: https://www.ncbi.nlm.nih.gov/books/NBK25497/


# 
# 
# class(testing)
# test_list[[1]]$QueryTranslation
# 
# test_list[[1]]$QueryTranslation
# 
# stringr::str_remove(string = test_list[[1]]$QueryTranslation, pattern = "([All Fields])|[(]|[)]/g")
# 
# stringr::str_remove(string = test_list[[1]]$QueryTranslation, pattern = "(\\[All Fields\\])|\\)|\\(")
# 
# stringr::str_remove(string = test_list[[1]]$QueryTranslation, pattern = "(\\[All Fields\\])")
# stringr::str_remove(string = test_list[[1]]$QueryTranslation, pattern = "\\(")
# stringr::str_remove(string = test_list[[1]]$QueryTranslation, pattern = "\\)")
# 
# gsub(pattern = "(\\[All Fields\\])|\\)|\\(", replacement = '', x = test_list[[1]]$QueryTranslation)
# 
# length(test_list[[1]]$QueryTranslation)
# 
# length(test_list[[5]]$ids)
# 
# 
#  or 
# entrez_dbs()
# entrez_db_searchable('gene')
# 
# xxx <- entrez_search(db="gene", term = 'BQ195046')
# 
# yyy <- entrez_summary(db="gene", id = c('AF123456.1', '16814', '22178', 'BQ195046'))
# yyy <- entrez_summary(db="gene", id = c('BQ195046'))
# 
# zzz <- entrez_link(db="protein", dbfrom="gene", id=test)
# 
# zzz <- entrez_fetch(db="gene", id = c('BQ195046', 'AF123456.1'), rettype = 'xml')
# 
# 
# zzz2 <- XML::xmlToList(zzz)
# rm(zzz2)
# 
# zzz2 <- jsonlite::fromJSON(zzz)
# zzz2 <- rjson::fromJSON(zzz)
# zzz2 <- RJSONIO::fromJSON(zzz)
# RJSONIO::isValidJSON(zzz, asText = T)
# 
