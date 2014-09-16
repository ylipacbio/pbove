#rootDir = "C:/Users/yli/Documents/R/HGAP/redo_C4_Baseline"
#outDir = "out_plot"
#all_out_fn = "all_out.txt"

require("ggplot2")
require("gridExtra")

main <- function(rootDir, all_out_fn, outDir) {
    setwd(rootDir)
    dir.create(file.path(rootDir, outDir), showWarnings = FALSE)
    dat = getTPFPData_from_all_csv_in_fn(all_out_fn)
    plot_all_items(dat, outDir)
    plot_all_groups(dat, outDir)
    plot_avg_groups(dat, outDir)
}


PrintGGPlot <- function(p, outFile="", type="") {
    if (outFile == "") {
        print (p)
    } else {
        if (type == "pdf" || type == "") {
            pdf(outFile, paper="a4")
            print(p)
            dev.off()
        } else if (type == "jpg") {
            jpeg(outFile, width=8,height=6,units="in",res=400, pointsize=1)
            print(p)
            dev.off()
        } else if (type == "png") {
            png(outFile, width=8,height=6,units="in",res=400, pointsize=1)
            print(p)
            dev.off()
        }
    }
}

plot_all_items <- function(dat, outDir) {
    # Simply plot all items (runs) in one plot
    p = ggplot(dat, aes(FDR, sensitivity, color=name)) +
        geom_line() + geom_point() +
        xlab ("False Discovery Rate") +
        ylab ("Sensitivity") +
        ggtitle("Performance in Preassembly Overlap Detection") +
        theme(plot.title=element_text(lineheight=.8, face="bold"))

    PrintGGPlot(p,
                paste(outDir, "/", "all", "sensitivityFDR.jpg", sep=""),
                "jpg")

    PrintGGPlot(p,
                paste(outDir, "/", "all", "sensitivityFDR.pdf", sep=""),
                "pdf")
}

plot_all_groups <- function(dat, outDir) {
    # Group items by 'group' and plot items of each group in a figure
    for (i in 1:dim(array(levels(dat$group)))[1]) {
        by_group = levels(dat$group)[i]
        s_dat = subset(dat, dat$group==by_group)
        title = paste("Performance in Preassembly Overlap Detection\n within Group: ", 
                      by_group, sep="")
        p = ggplot(s_dat, aes(FDR, sensitivity, color=name)) +
            geom_line() + geom_point() +
            xlab ("False Discovery Rate") +
            ylab ("Sensitivity") + 
            ggtitle(title) + 
            theme(plot.title=element_text(lineheight=.8, face="bold"))

        PrintGGPlot(p,
                    paste(outDir, "/", by_group, "_sensitivityFDR.jpg", sep=""),
                    "jpg")

        PrintGGPlot(p,
                    paste(outDir, "/", by_group, "_sensitivityFDR.pdf", sep=""),
                    "pdf")
    }
}

plot_avg_groups <- function(dat, outDir, sensitivity_step_sz=0.02) {
    # Group items by 'group' and compute the mean FDR at each
    # sensitivity range, e.g., FDR of all runs of group "P5C3"
    # at sensitivity range (0.0, 0.02).
    init_avgFDRs = 0
    for (i in 1:dim(array(levels(dat$group)))[1]) {
        s_dat = subset(dat, dat$group==levels(dat$group)[i])
        by_group = toString(levels(dat$group)[i])
        print (paste("group=", by_group))
        for (sensitivity_lb in seq(0, max(s_dat$sensitivity), sensitivity_step_sz)) {
            sensitivity_ub = sensitivity_lb + sensitivity_step_sz
            ss_dat = subset(s_dat,
                            s_dat$sensitivity >= sensitivity_lb & s_dat$sensitivity < sensitivity_ub)
            num_name = 0 # number of runs which have data point within the
            # specified sensitivity range [lb, ub)
            init_res = 0
            print (paste('   ', sensitivity_lb, sensitivity_ub, ' '))

            for (j in 1:dim(array(levels(ss_dat$name)))[1]) {
                sss_dat = subset(ss_dat, ss_dat$name==levels(ss_dat$name)[j])
                name = toString(levels(ss_dat$name)[j])
                # sss_dat, rows of which name is the specified name
                if (nrow(sss_dat) > 0) {
                    num_name = num_name + 1
                    # Compute average FDR within specified sensitivity range
                    mean_FDR = mean(as.numeric(sss_dat$FDR))
                    if (init_res == 0) {
                        res = cbind(name, mean_FDR)
                        init_res = 1
                    } else {
                        res = rbind(res, cbind(name, mean_FDR))
                    }
                } else {
                    # No data point of run with name 'name' falls into
                    # the specified sensitivity range, pass
                }
            }
            # At this point, compute average of (averge FDR of a run with
            # name 'name' at specified sensitivity range) as average
            # FDR at specified sensitivity range
            avgFDR_at_sensitivity = mean(as.numeric(res[, 'mean_FDR']))
            x = cbind(by_group, sensitivity_lb, sensitivity_ub, avgFDR_at_sensitivity)
            if (init_avgFDRs == 0) {
                avgFDRs = x
                init_avgFDRs = 1
            } else {
                avgFDRs = rbind(avgFDRs, x)
            }
        }
    }
    avgFDRs = data.frame(avgFDRs)
    names(avgFDRs) = c("group", 'sensitivity_lb', 'sensitivity_ub', 'avgFDR')
    avgFDRs$sensitivity_range = as.factor(paste(avgFDRs$sensitivity_lb, "_", avgFDRs$sensitivity_ub, sep=""))
}

getTPFPData_from_all_csv_in_fn <- function(all_out_fn) {
    # all_out_fn is a string representing a file name.
    # The all_out_fn file should contain many csv results from pbove,
    # each line of which should have three fields:
    # out_csv, name, group
    all_out = read.table(all_out_fn, header=T)
    for (i in 1:dim(all_out)[1]) {
        if (i == 1) {
            dat = getTPFPData(all_out[i,])
        } else {
            dat = rbind(dat, getTPFPData(all_out[i,]))
        }
    }
    dat
}


getTPFPData <- function(all_out_i) {
    # all_out_i has three columns: out_csv, name, group
    out_csv = all_out_i$out_csv
    name    = all_out_i$name
    group   = all_out_i$group

    tb = read.table(toString(out_csv), header=T)
    tb$name = rep(name, nrow(tb))
    tb$group = rep(group, nrow(tb))
    tb$sensitivity <- tb$numTP/(tb$numTP + tb$numFN)
    tb$specificity <- tb$numTN/(tb$numFP + tb$numTN)
    tb$FPR <- 1.0 - tb$specificity
    tb$FDR <- tb$numFP/(tb$numFP + tb$numTP)
# [1] "ScoreCutoff"       "numTP"             "numFP"
# [4] "numFN"             "numTN"             "numPW"
# [7] "numNW"             "numGTPos"          "numGTNeg"
#[10] "numGTWeak"         "numPredPos"        "numPredNeg"
#[13] "numUnmappableAlns" "numMappableAlns"   "numAlns"
#[16] "name"              "group"             "sensitivity"
#[19] "specificity"       "FPR"               "FDR"
    tb
}

main(rootDir, all_out_fn, outDir)
