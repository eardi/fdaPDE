expand.grid.alt <- function(seq1,seq2) {
cbind(rep.int(seq1, length(seq2)),
c(t(matrix(rep.int(seq2, length(seq1)), nrow=length(seq2)))))
}