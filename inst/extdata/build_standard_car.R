# Build standard car tibble with cmpound and ms2 information
# 250408
# Barry Song

# Carnitine compound information
standard_car <- dplyr::as_tibble(openxlsx::read.xlsx("./inst/extdata/standard_car.xlsx", sheet = 1))
# HMDB ms2
hmdbMsTb <- pubmsR::load_hmdbMsTb(standard = TRUE)
# lipidblast ms2
lipidblastMsTb <- pubmsR::load_lipidblastMsTb(standard = TRUE)
# MoNA ms2
monaMsTb <- pubmsR::load_monaMsTb(standard = TRUE)
# TangLab inhouse ms2
load("D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/build_package/MetaboLib.ms2/inHouseMs2List.RData")
inHouseAmideMsTb <- inHouseMs2List$Amide_Pos$ms2
inHouseT3MsTb <- inHouseMs2List$T3_Pos$ms2

# Prioritize spectra from in-house database, followed by HMDB, MoNA, and LipidBlast in that order.
# Merge multiple spectra
mzList <- list()
intensityList <- list()
# 1. HMDB0000062 (L-Carnitine)
{
  i <- 1
  inchikey <- stringr::str_split(standard_car[i, ]$inchikey, pattern = "-")[[1]][1]
  hmdb_i <- hmdbMsTb[which(stringr::str_detect(hmdbMsTb$inchikey, pattern = inchikey) & hmdbMsTb$polarity == "pos"), ]
  lipidblast_i <- lipidblastMsTb[which(stringr::str_detect(lipidblastMsTb$inchikey, pattern = inchikey) & lipidblastMsTb$polarity == "pos"), ]
  mona_i <- monaMsTb[which(stringr::str_detect(monaMsTb$inchikey, pattern = inchikey) & monaMsTb$polarity == "pos"), ]
  Amide_i <- inHouseAmideMsTb[which(stringr::str_detect(inHouseAmideMsTb$inchikey, pattern = inchikey)), ]
  T3_i <- inHouseT3MsTb[which(stringr::str_detect(inHouseT3MsTb$inchikey, pattern = inchikey)), ]

  standardRows <- rbind(Amide_i, T3_i)
  # standardRows <- hmdb_i
  # standardRows <- mona_i
  # standardRows <- lipidblast_i
  spMatList <- lapply(1:nrow(standardRows), function(i) {
    MetaboSpectra::get_spMat(standardRows[i, ])
  })
  spMat <- MetaboSpectra::merge_spectra(spMatList = spMatList, method = sum, digits = 0)
  MetaboSpectra::plotSpectra(spMat)
  mzList <- append(mzList, list(spMat[, "mz"]))
  intensityList <- append(intensityList, list(spMat[, "intensity"]))
}
# 2. HMDB0000201 (L-Acetylcarnitine)
{
  i <- 2
  inchikey <- stringr::str_split(standard_car[i, ]$inchikey, pattern = "-")[[1]][1]
  hmdb_i <- hmdbMsTb[which(stringr::str_detect(hmdbMsTb$inchikey, pattern = inchikey) & hmdbMsTb$polarity == "pos"), ]
  lipidblast_i <- lipidblastMsTb[which(stringr::str_detect(lipidblastMsTb$inchikey, pattern = inchikey) & lipidblastMsTb$polarity == "pos"), ]
  mona_i <- monaMsTb[which(stringr::str_detect(monaMsTb$inchikey, pattern = inchikey) & monaMsTb$polarity == "pos"), ]
  Amide_i <- inHouseAmideMsTb[which(stringr::str_detect(inHouseAmideMsTb$inchikey, pattern = inchikey)), ]
  T3_i <- inHouseT3MsTb[which(stringr::str_detect(inHouseT3MsTb$inchikey, pattern = inchikey)), ]

  standardRows <- rbind(Amide_i, T3_i)
  # standardRows <- hmdb_i
  # standardRows <- mona_i
  # standardRows <- lipidblast_i
  spMatList <- lapply(1:nrow(standardRows), function(i) {
    MetaboSpectra::get_spMat(standardRows[i, ])
  })
  spMat <- MetaboSpectra::merge_spectra(spMatList = spMatList, method = sum, digits = 0)
  MetaboSpectra::plotSpectra(spMat)
  mzList <- append(mzList, list(spMat[, "mz"]))
  intensityList <- append(intensityList, list(spMat[, "intensity"]))
}
# 3. HMDB0000824 (Propionylcarnitine)
{
  i <- 3
  inchikey <- stringr::str_split(standard_car[i, ]$inchikey, pattern = "-")[[1]][1]
  hmdb_i <- hmdbMsTb[which(stringr::str_detect(hmdbMsTb$inchikey, pattern = inchikey) & hmdbMsTb$polarity == "pos"), ]
  lipidblast_i <- lipidblastMsTb[which(stringr::str_detect(lipidblastMsTb$inchikey, pattern = inchikey) & lipidblastMsTb$polarity == "pos"), ]
  mona_i <- monaMsTb[which(stringr::str_detect(monaMsTb$inchikey, pattern = inchikey) & monaMsTb$polarity == "pos"), ]
  Amide_i <- inHouseAmideMsTb[which(stringr::str_detect(inHouseAmideMsTb$inchikey, pattern = inchikey)), ]
  T3_i <- inHouseT3MsTb[which(stringr::str_detect(inHouseT3MsTb$inchikey, pattern = inchikey)), ]

  standardRows <- rbind(Amide_i, T3_i)
  # standardRows <- hmdb_i
  # standardRows <- mona_i
  # standardRows <- lipidblast_i
  spMatList <- lapply(1:nrow(standardRows), function(i) {
    MetaboSpectra::get_spMat(standardRows[i, ])
  })
  spMat <- MetaboSpectra::merge_spectra(spMatList = spMatList, method = sum, digits = 0)
  MetaboSpectra::plotSpectra(spMat)
  mzList <- append(mzList, list(spMat[, "mz"]))
  intensityList <- append(intensityList, list(spMat[, "intensity"]))
}
# 4. HMDB0002013 (Butyrylcarnitine)
{
  i <- 4
  inchikey <- stringr::str_split(standard_car[i, ]$inchikey, pattern = "-")[[1]][1]
  hmdb_i <- hmdbMsTb[which(stringr::str_detect(hmdbMsTb$inchikey, pattern = inchikey) & hmdbMsTb$polarity == "pos"), ]
  lipidblast_i <- lipidblastMsTb[which(stringr::str_detect(lipidblastMsTb$inchikey, pattern = inchikey) & lipidblastMsTb$polarity == "pos"), ]
  mona_i <- monaMsTb[which(stringr::str_detect(monaMsTb$inchikey, pattern = inchikey) & monaMsTb$polarity == "pos"), ]
  Amide_i <- inHouseAmideMsTb[which(stringr::str_detect(inHouseAmideMsTb$inchikey, pattern = inchikey)), ]
  T3_i <- inHouseT3MsTb[which(stringr::str_detect(inHouseT3MsTb$inchikey, pattern = inchikey)), ]

  standardRows <- rbind(Amide_i, T3_i)
  # standardRows <- hmdb_i
  # standardRows <- mona_i
  # standardRows <- lipidblast_i
  spMatList <- lapply(1:nrow(standardRows), function(i) {
    MetaboSpectra::get_spMat(standardRows[i, ])
  })
  spMat <- MetaboSpectra::merge_spectra(spMatList = spMatList, method = sum, digits = 0)
  MetaboSpectra::plotSpectra(spMat)
  mzList <- append(mzList, list(spMat[, "mz"]))
  intensityList <- append(intensityList, list(spMat[, "intensity"]))
}
# 5. HMDB0013128 (Valerylcarnitine)
{
  i <- 5
  inchikey <- stringr::str_split(standard_car[i, ]$inchikey, pattern = "-")[[1]][1]
  hmdb_i <- hmdbMsTb[which(stringr::str_detect(hmdbMsTb$inchikey, pattern = inchikey) & hmdbMsTb$polarity == "pos"), ]
  lipidblast_i <- lipidblastMsTb[which(stringr::str_detect(lipidblastMsTb$inchikey, pattern = inchikey) & lipidblastMsTb$polarity == "pos"), ]
  mona_i <- monaMsTb[which(stringr::str_detect(monaMsTb$inchikey, pattern = inchikey) & monaMsTb$polarity == "pos"), ]
  Amide_i <- inHouseAmideMsTb[which(stringr::str_detect(inHouseAmideMsTb$inchikey, pattern = inchikey)), ]
  T3_i <- inHouseT3MsTb[which(stringr::str_detect(inHouseT3MsTb$inchikey, pattern = inchikey)), ]

  standardRows <- rbind(Amide_i, T3_i)
  # standardRows <- hmdb_i
  # standardRows <- mona_i
  # standardRows <- lipidblast_i
  spMatList <- lapply(1:nrow(standardRows), function(i) {
    MetaboSpectra::get_spMat(standardRows[i, ])
  })
  spMat <- MetaboSpectra::merge_spectra(spMatList = spMatList, method = sum, digits = 0)
  MetaboSpectra::plotSpectra(spMat)
  mzList <- append(mzList, list(spMat[, "mz"]))
  intensityList <- append(intensityList, list(spMat[, "intensity"]))
}
# 6. HMDB0000756 (Hexanoylcarnitine)
{
  i <- 6
  inchikey <- stringr::str_split(standard_car[i, ]$inchikey, pattern = "-")[[1]][1]
  hmdb_i <- hmdbMsTb[which(stringr::str_detect(hmdbMsTb$inchikey, pattern = inchikey) & hmdbMsTb$polarity == "pos"), ]
  lipidblast_i <- lipidblastMsTb[which(stringr::str_detect(lipidblastMsTb$inchikey, pattern = inchikey) & lipidblastMsTb$polarity == "pos"), ]
  mona_i <- monaMsTb[which(stringr::str_detect(monaMsTb$inchikey, pattern = inchikey) & monaMsTb$polarity == "pos"), ]
  Amide_i <- inHouseAmideMsTb[which(stringr::str_detect(inHouseAmideMsTb$inchikey, pattern = inchikey)), ]
  T3_i <- inHouseT3MsTb[which(stringr::str_detect(inHouseT3MsTb$inchikey, pattern = inchikey)), ]

  standardRows <- rbind(Amide_i, T3_i)
  # standardRows <- hmdb_i
  # standardRows <- mona_i
  # standardRows <- lipidblast_i
  spMatList <- lapply(1:nrow(standardRows), function(i) {
    MetaboSpectra::get_spMat(standardRows[i, ])
  })
  spMat <- MetaboSpectra::merge_spectra(spMatList = spMatList, method = sum, digits = 0)
  MetaboSpectra::plotSpectra(spMat)
  mzList <- append(mzList, list(spMat[, "mz"]))
  intensityList <- append(intensityList, list(spMat[, "intensity"]))
}
# 7. HMDB0000791 (Octanoylcarnitine)
{
  i <- 7
  inchikey <- stringr::str_split(standard_car[i, ]$inchikey, pattern = "-")[[1]][1]
  hmdb_i <- hmdbMsTb[which(stringr::str_detect(hmdbMsTb$inchikey, pattern = inchikey) & hmdbMsTb$polarity == "pos"), ]
  lipidblast_i <- lipidblastMsTb[which(stringr::str_detect(lipidblastMsTb$inchikey, pattern = inchikey) & lipidblastMsTb$polarity == "pos"), ]
  mona_i <- monaMsTb[which(stringr::str_detect(monaMsTb$inchikey, pattern = inchikey) & monaMsTb$polarity == "pos"), ]
  Amide_i <- inHouseAmideMsTb[which(stringr::str_detect(inHouseAmideMsTb$inchikey, pattern = inchikey)), ]
  T3_i <- inHouseT3MsTb[which(stringr::str_detect(inHouseT3MsTb$inchikey, pattern = inchikey)), ]

  standardRows <- rbind(Amide_i, T3_i)
  # standardRows <- hmdb_i
  # standardRows <- mona_i
  # standardRows <- lipidblast_i
  spMatList <- lapply(1:nrow(standardRows), function(i) {
    MetaboSpectra::get_spMat(standardRows[i, ])
  })
  spMat <- MetaboSpectra::merge_spectra(spMatList = spMatList, method = sum, digits = 0)
  MetaboSpectra::plotSpectra(spMat)
  mzList <- append(mzList, list(spMat[, "mz"]))
  intensityList <- append(intensityList, list(spMat[, "intensity"]))
}
# 8. HMDB0000651 (Decanoylcarnitine)
{
  i <- 8
  inchikey <- stringr::str_split(standard_car[i, ]$inchikey, pattern = "-")[[1]][1]
  hmdb_i <- hmdbMsTb[which(stringr::str_detect(hmdbMsTb$inchikey, pattern = inchikey) & hmdbMsTb$polarity == "pos"), ]
  lipidblast_i <- lipidblastMsTb[which(stringr::str_detect(lipidblastMsTb$inchikey, pattern = inchikey) & lipidblastMsTb$polarity == "pos"), ]
  mona_i <- monaMsTb[which(stringr::str_detect(monaMsTb$inchikey, pattern = inchikey) & monaMsTb$polarity == "pos"), ]
  Amide_i <- inHouseAmideMsTb[which(stringr::str_detect(inHouseAmideMsTb$inchikey, pattern = inchikey)), ]
  T3_i <- inHouseT3MsTb[which(stringr::str_detect(inHouseT3MsTb$inchikey, pattern = inchikey)), ]

  standardRows <- rbind(Amide_i, T3_i)
  # standardRows <- hmdb_i
  # standardRows <- mona_i
  # standardRows <- lipidblast_i
  spMatList <- lapply(1:nrow(standardRows), function(i) {
    MetaboSpectra::get_spMat(standardRows[i, ])
  })
  spMat <- MetaboSpectra::merge_spectra(spMatList = spMatList, method = sum, digits = 0)
  MetaboSpectra::plotSpectra(spMat)
  mzList <- append(mzList, list(spMat[, "mz"]))
  intensityList <- append(intensityList, list(spMat[, "intensity"]))
}
# 9. HMDB0002250 (Dodecanoylcarnitine)
{
  i <- 9
  inchikey <- stringr::str_split(standard_car[i, ]$inchikey, pattern = "-")[[1]][1]
  hmdb_i <- hmdbMsTb[which(stringr::str_detect(hmdbMsTb$inchikey, pattern = inchikey) & hmdbMsTb$polarity == "pos"), ]
  lipidblast_i <- lipidblastMsTb[which(stringr::str_detect(lipidblastMsTb$inchikey, pattern = inchikey) & lipidblastMsTb$polarity == "pos"), ]
  mona_i <- monaMsTb[which(stringr::str_detect(monaMsTb$inchikey, pattern = inchikey) & monaMsTb$polarity == "pos"), ]
  Amide_i <- inHouseAmideMsTb[which(stringr::str_detect(inHouseAmideMsTb$inchikey, pattern = inchikey)), ]
  T3_i <- inHouseT3MsTb[which(stringr::str_detect(inHouseT3MsTb$inchikey, pattern = inchikey)), ]

  standardRows <- rbind(Amide_i, T3_i)
  # standardRows <- hmdb_i
  # standardRows <- mona_i
  # standardRows <- lipidblast_i
  spMatList <- lapply(1:nrow(standardRows), function(i) {
    MetaboSpectra::get_spMat(standardRows[i, ])
  })
  spMat <- MetaboSpectra::merge_spectra(spMatList = spMatList, method = sum, digits = 0)
  MetaboSpectra::plotSpectra(spMat)
  mzList <- append(mzList, list(spMat[, "mz"]))
  intensityList <- append(intensityList, list(spMat[, "intensity"]))
}
# 10. HMDB0005066 (Tetradecanoylcarnitine)
{
  i <- 10
  inchikey <- stringr::str_split(standard_car[i, ]$inchikey, pattern = "-")[[1]][1]
  hmdb_i <- hmdbMsTb[which(stringr::str_detect(hmdbMsTb$inchikey, pattern = inchikey) & hmdbMsTb$polarity == "pos"), ]
  lipidblast_i <- lipidblastMsTb[which(stringr::str_detect(lipidblastMsTb$inchikey, pattern = inchikey) & lipidblastMsTb$polarity == "pos"), ]
  mona_i <- monaMsTb[which(stringr::str_detect(monaMsTb$inchikey, pattern = inchikey) & monaMsTb$polarity == "pos"), ]
  Amide_i <- inHouseAmideMsTb[which(stringr::str_detect(inHouseAmideMsTb$inchikey, pattern = inchikey)), ]
  T3_i <- inHouseT3MsTb[which(stringr::str_detect(inHouseT3MsTb$inchikey, pattern = inchikey)), ]

  standardRows <- rbind(Amide_i, T3_i)
  # standardRows <- hmdb_i
  # standardRows <- mona_i
  # standardRows <- lipidblast_i
  spMatList <- lapply(1:nrow(standardRows), function(i) {
    MetaboSpectra::get_spMat(standardRows[i, ])
  })
  spMat <- MetaboSpectra::merge_spectra(spMatList = spMatList, method = sum, digits = 0)
  MetaboSpectra::plotSpectra(spMat)
  mzList <- append(mzList, list(spMat[, "mz"]))
  intensityList <- append(intensityList, list(spMat[, "intensity"]))
}
# 11. HMDB0000222 (Palmitoylcarnitine)
{
  i <- 11
  inchikey <- stringr::str_split(standard_car[i, ]$inchikey, pattern = "-")[[1]][1]
  hmdb_i <- hmdbMsTb[which(stringr::str_detect(hmdbMsTb$inchikey, pattern = inchikey) & hmdbMsTb$polarity == "pos"), ]
  lipidblast_i <- lipidblastMsTb[which(stringr::str_detect(lipidblastMsTb$inchikey, pattern = inchikey) & lipidblastMsTb$polarity == "pos"), ]
  mona_i <- monaMsTb[which(stringr::str_detect(monaMsTb$inchikey, pattern = inchikey) & monaMsTb$polarity == "pos"), ]
  Amide_i <- inHouseAmideMsTb[which(stringr::str_detect(inHouseAmideMsTb$inchikey, pattern = inchikey)), ]
  T3_i <- inHouseT3MsTb[which(stringr::str_detect(inHouseT3MsTb$inchikey, pattern = inchikey)), ]

  standardRows <- rbind(Amide_i, T3_i)
  # standardRows <- hmdb_i
  # standardRows <- mona_i
  # standardRows <- lipidblast_i
  spMatList <- lapply(1:nrow(standardRows), function(i) {
    MetaboSpectra::get_spMat(standardRows[i, ])
  })
  spMat <- MetaboSpectra::merge_spectra(spMatList = spMatList, method = sum, digits = 0)
  MetaboSpectra::plotSpectra(spMat)
  mzList <- append(mzList, list(spMat[, "mz"]))
  intensityList <- append(intensityList, list(spMat[, "intensity"]))
}
# 12. HMDB0000848 (Stearoylcarnitine)
{
  i <- 12
  inchikey <- stringr::str_split(standard_car[i, ]$inchikey, pattern = "-")[[1]][1]
  hmdb_i <- hmdbMsTb[which(stringr::str_detect(hmdbMsTb$inchikey, pattern = inchikey) & hmdbMsTb$polarity == "pos"), ]
  lipidblast_i <- lipidblastMsTb[which(stringr::str_detect(lipidblastMsTb$inchikey, pattern = inchikey) & lipidblastMsTb$polarity == "pos"), ]
  mona_i <- monaMsTb[which(stringr::str_detect(monaMsTb$inchikey, pattern = inchikey) & monaMsTb$polarity == "pos"), ]
  Amide_i <- inHouseAmideMsTb[which(stringr::str_detect(inHouseAmideMsTb$inchikey, pattern = inchikey)), ]
  T3_i <- inHouseT3MsTb[which(stringr::str_detect(inHouseT3MsTb$inchikey, pattern = inchikey)), ]

  standardRows <- rbind(Amide_i, T3_i)
  # standardRows <- hmdb_i
  # standardRows <- mona_i
  # standardRows <- lipidblast_i
  spMatList <- lapply(1:nrow(standardRows), function(i) {
    MetaboSpectra::get_spMat(standardRows[i, ])
  })
  spMat <- MetaboSpectra::merge_spectra(spMatList = spMatList, method = sum, digits = 0)
  MetaboSpectra::plotSpectra(spMat)
  mzList <- append(mzList, list(spMat[, "mz"]))
  intensityList <- append(intensityList, list(spMat[, "intensity"]))
}
# 13. HMDB0006460 (Arachidyl carnitine)
{
  i <- 13
  inchikey <- stringr::str_split(standard_car[i, ]$inchikey, pattern = "-")[[1]][1]
  hmdb_i <- hmdbMsTb[which(stringr::str_detect(hmdbMsTb$inchikey, pattern = inchikey) & hmdbMsTb$polarity == "pos"), ]
  lipidblast_i <- lipidblastMsTb[which(stringr::str_detect(lipidblastMsTb$inchikey, pattern = inchikey) & lipidblastMsTb$polarity == "pos"), ]
  mona_i <- monaMsTb[which(stringr::str_detect(monaMsTb$inchikey, pattern = inchikey) & monaMsTb$polarity == "pos"), ]
  Amide_i <- inHouseAmideMsTb[which(stringr::str_detect(inHouseAmideMsTb$inchikey, pattern = inchikey)), ]
  T3_i <- inHouseT3MsTb[which(stringr::str_detect(inHouseT3MsTb$inchikey, pattern = inchikey)), ]

  standardRows <- rbind(Amide_i, T3_i)
  # standardRows <- hmdb_i
  # standardRows <- mona_i
  # standardRows <- lipidblast_i
  spMatList <- lapply(1:nrow(standardRows), function(i) {
    MetaboSpectra::get_spMat(standardRows[i, ])
  })
  spMat <- MetaboSpectra::merge_spectra(spMatList = spMatList, method = sum, digits = 0)
  MetaboSpectra::plotSpectra(spMat)
  mzList <- append(mzList, list(spMat[, "mz"]))
  intensityList <- append(intensityList, list(spMat[, "intensity"]))
}
# 14. HMDB0062468 (Docosanoylcarnitine)
{
  i <- 14
  inchikey <- stringr::str_split(standard_car[i, ]$inchikey, pattern = "-")[[1]][1]
  hmdb_i <- hmdbMsTb[which(stringr::str_detect(hmdbMsTb$inchikey, pattern = inchikey) & hmdbMsTb$polarity == "pos"), ]
  lipidblast_i <- lipidblastMsTb[which(stringr::str_detect(lipidblastMsTb$inchikey, pattern = inchikey) & lipidblastMsTb$polarity == "pos"), ]
  mona_i <- monaMsTb[which(stringr::str_detect(monaMsTb$inchikey, pattern = inchikey) & monaMsTb$polarity == "pos"), ]
  Amide_i <- inHouseAmideMsTb[which(stringr::str_detect(inHouseAmideMsTb$inchikey, pattern = inchikey)), ]
  T3_i <- inHouseT3MsTb[which(stringr::str_detect(inHouseT3MsTb$inchikey, pattern = inchikey)), ]

  # standardRows <- rbind(Amide_i, T3_i)
  # standardRows <- hmdb_i
  # standardRows <- mona_i
  standardRows <- lipidblast_i
  spMatList <- lapply(1:nrow(standardRows), function(i) {
    MetaboSpectra::get_spMat(standardRows[i, ])
  })
  spMat <- MetaboSpectra::merge_spectra(spMatList = spMatList, method = sum, digits = 0)
  MetaboSpectra::plotSpectra(spMat)
  mzList <- append(mzList, list(spMat[, "mz"]))
  intensityList <- append(intensityList, list(spMat[, "intensity"]))
}
# 15. HMDB0240665 (Lignoceroylcarnitine)
{
  i <- 15
  inchikey <- stringr::str_split(standard_car[i, ]$inchikey, pattern = "-")[[1]][1]
  hmdb_i <- hmdbMsTb[which(stringr::str_detect(hmdbMsTb$inchikey, pattern = inchikey) & hmdbMsTb$polarity == "pos"), ]
  lipidblast_i <- lipidblastMsTb[which(stringr::str_detect(lipidblastMsTb$inchikey, pattern = inchikey) & lipidblastMsTb$polarity == "pos"), ]
  mona_i <- monaMsTb[which(stringr::str_detect(monaMsTb$inchikey, pattern = inchikey) & monaMsTb$polarity == "pos"), ]
  Amide_i <- inHouseAmideMsTb[which(stringr::str_detect(inHouseAmideMsTb$inchikey, pattern = inchikey)), ]
  T3_i <- inHouseT3MsTb[which(stringr::str_detect(inHouseT3MsTb$inchikey, pattern = inchikey)), ]

  # standardRows <- rbind(Amide_i, T3_i)
  # standardRows <- hmdb_i
  # standardRows <- mona_i
  standardRows <- lipidblast_i
  spMatList <- lapply(1:nrow(standardRows), function(i) {
    MetaboSpectra::get_spMat(standardRows[i, ])
  })
  spMat <- MetaboSpectra::merge_spectra(spMatList = spMatList, method = sum, digits = 0)
  MetaboSpectra::plotSpectra(spMat)
  mzList <- append(mzList, list(spMat[, "mz"]))
  intensityList <- append(intensityList, list(spMat[, "intensity"]))
}
# 16. HMDB0006347 (Hexacosanoyl carnitine)
{
  i <- 16
  inchikey <- stringr::str_split(standard_car[i, ]$inchikey, pattern = "-")[[1]][1]
  hmdb_i <- hmdbMsTb[which(stringr::str_detect(hmdbMsTb$inchikey, pattern = inchikey) & hmdbMsTb$polarity == "pos"), ]
  lipidblast_i <- lipidblastMsTb[which(stringr::str_detect(lipidblastMsTb$inchikey, pattern = inchikey) & lipidblastMsTb$polarity == "pos"), ]
  mona_i <- monaMsTb[which(stringr::str_detect(monaMsTb$inchikey, pattern = inchikey) & monaMsTb$polarity == "pos"), ]
  Amide_i <- inHouseAmideMsTb[which(stringr::str_detect(inHouseAmideMsTb$inchikey, pattern = inchikey)), ]
  T3_i <- inHouseT3MsTb[which(stringr::str_detect(inHouseT3MsTb$inchikey, pattern = inchikey)), ]

  standardRows <- rbind(Amide_i, T3_i)
  # standardRows <- hmdb_i
  # standardRows <- mona_i
  # standardRows <- lipidblast_i
  spMatList <- lapply(1:nrow(standardRows), function(i) {
    MetaboSpectra::get_spMat(standardRows[i, ])
  })
  spMat <- MetaboSpectra::merge_spectra(spMatList = spMatList, method = sum, digits = 0)
  MetaboSpectra::plotSpectra(spMat)
  mzList <- append(mzList, list(spMat[, "mz"]))
  intensityList <- append(intensityList, list(spMat[, "intensity"]))
}

MetaboSpectra::plotSpectra(MetaboSpectra::get_spMat(hmdb_i[1, ]))
MetaboSpectra::plotSpectra(MetaboSpectra::get_spMat(lipidblast_i[1, ]))
MetaboSpectra::plotSpectra(MetaboSpectra::get_spMat(mona_i[3, ]))
MetaboSpectra::plotSpectra(MetaboSpectra::get_spMat(Amide_i[1, ]))
MetaboSpectra::plotSpectra(MetaboSpectra::get_spMat(T3_i[1, ]))

standard_car$mz <- mzList
standard_car$intensity <- intensityList

usethis::use_data(standard_car)
