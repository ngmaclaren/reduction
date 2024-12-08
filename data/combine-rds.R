## Do not run this file from start to finish without checking what's happening.

                                        # Need to combine marker_cafe and prosper

network <- "prosper"
filepattern <- paste0("fullstate-", network, "-")#, "[:alpha:]+", ".rds")
outfilename <- paste0("fullstate-", network, ".rds")

filenames <- list.files()[grep(filepattern, list.files())]

                                        # HARD-CODED!!!
if(any(grepl("first", filenames))) {
    a <- readRDS(grep("first", filenames, value = TRUE))
    b <- readRDS(grep("last", filenames, value = TRUE))
    c <- rbind(a$mutualistic, b$mutualistic)
    rownames(c) <- seq(nrow(c))
    d <- list(mutualistic = c)

    newfilename <- paste0("fullstate-", network, "-mutualistic.rds")
    saveRDS(d, newfilename)
    
    filenames <- filenames[-grep("first|last", filenames)]
    filenames <- c(filenames, newfilename)
}

rdss <- lapply(filenames, readRDS)
rdss <- unlist(rdss, recursive = FALSE)
rdss <- rdss[c("doublewell", "SIS", "genereg", "mutualistic")]

saveRDS(rdss, outfilename)


                                        # check
fs <- readRDS(outfilename)
other <- readRDS("fullstate-celegans.rds")

all.equal(names(fs), names(other))

sapply(fs, dim)
sapply(other, dim)

all.equal(fs, readRDS(paste0("./old/", outfilename)))
