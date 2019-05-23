setwd("~/git/caisar/output")

a = colMeans(matrix(read.table("Scenario1_time1.txt", sep = ",")[,2], nrow = 20))[c(1,3:10)]
cat(paste(round(a, digits = 2), collapse = " & "))

a = colMeans(matrix(read.table("Scenario1_time3.txt", sep = ",")[,2], nrow = 20))[c(1,3:10)]
cat(paste(round(a, digits = 2), collapse = " & "))

a = colMeans(matrix(read.table("Scenario1_time5.txt", sep = ",")[,2], nrow = 20))[c(1,3:10)]
cat(paste(round(a, digits = 2), collapse = " & "))

a = colMeans(matrix(read.table("Scenario2_time1.txt", sep = ",")[,2], nrow = 20))[c(2,3:10)]
a[1] = a[1] + a[5]
cat(paste(round(a, digits = 2), collapse = " & "))

a = colMeans(matrix(read.table("Scenario2_time3.txt", sep = ",")[,2], nrow = 20))[c(2,3:10)]
a[1] = a[1] + a[5]
cat(paste(round(a, digits = 2), collapse = " & "))

a = colMeans(matrix(read.table("Scenario2_time5.txt", sep = ",")[,2], nrow = 20))[c(2,3:10)]
a[1] = a[1] + a[5]
cat(paste(round(a, digits = 2), collapse = " & "))

a = colMeans(matrix(read.table("Scenario4_time1.txt", sep = ",")[,2], nrow = 20)) * 5.522
cat(paste(round(a, digits = 2), collapse = " & "))

a = colMeans(matrix(read.table("Scenario4_time2.txt", sep = ",")[,2], nrow = 20)) * 1.323
cat(paste(round(a, digits = 2), collapse = " & "))

a = colMeans(matrix(read.table("Scenario4_time3.txt", sep = ",")[,2], nrow = 20)) * 1.595
cat(paste(round(a, digits = 2), collapse = " & "))

a = colMeans(matrix(read.table("Scenario3_time1.txt", sep = ",")[,2], nrow = 20))[2:13] * 0.043
a[1] = a[1] + a[12]
cat(paste(round(a, digits = 3), collapse = " & "))

a = colMeans(matrix(read.table("Scenario3_time2.txt", sep = ",")[,2], nrow = 20))[2:13] * 0.043
a[1] = a[1] + a[12]
cat(paste(round(a, digits = 3), collapse = " & "))
