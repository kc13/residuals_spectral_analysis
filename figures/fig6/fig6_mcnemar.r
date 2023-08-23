#fig6_mcnemar.r
# for alpha beta incidence

M <-
matrix(c(31, 6, 35, 98),
       nrow = 2,
       dimnames = list("res" = c("T", "F"),
                       "shuf" = c("T", "F")))

print("corrected and uncorrected:")
mcnemar.test(M)
mcnemar.test(M,correct = FALSE)





