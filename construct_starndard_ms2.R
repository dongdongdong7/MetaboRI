load("D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/build_package/MetaboLib.ms2/inHouseMs2List.RData")
data("standard_car")
# T3 Pos
standard_car_ms2 <- inHouseMs2List$T3_Pos$ms2[which(inHouseMs2List$T3_Pos$ms2$accession %in% standard_car$id), ]
usethis::use_data(standard_car_ms2)

load("D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/build_package/MetaboLib.ms2/hmdbMs2List.RData")
data("standard_dic")
standard_dic_ms2 <- inHouseMs2List$T3_Neg$ms2[which(inHouseMs2List$T3_Neg$ms2$accession %in% standard_dic$id), ]
usethis::use_data(standard_dic_ms2)
