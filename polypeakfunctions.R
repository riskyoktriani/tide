cleanstring <- function(string) {
  string <- gsub("^>.*?\n", "", string)
  string <- toupper(string)
  string <- gsub("[^ACGTRYSWKMBDHVN]", "", string, perl=TRUE)
  return(string)
}
