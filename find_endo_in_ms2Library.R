# This script is used to search for acylcarnitine in ms2 libraries
# HMDB
load("D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/build_package/MetaboLib.ms2/hmdbMs2List.RData")
hmdbMs2List$cmps[hmdbMs2List$cmps$compound_id == "HMDB0000062", ]
hmdbMs2List$ms2[hmdbMs2List$ms2$accession == "HMDB0000062", ]
