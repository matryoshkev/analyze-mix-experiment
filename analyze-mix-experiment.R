# analyze-mix-experiment.R v0.4
# Script for use with the R software environment for statistical computing
# jeff smith and R. Fredrick Inglis 2019

# TO DO: 
# 
# Open questions
# - should output include calculated initial and final states?
# - should plots include zero/Inf data? 
# 
# Plots
# - labels every 3 logs when necessary (esp within) ex: 
# - breaks = 10^seq(-6, 6, by = 2)
# - breaks = 10^seq(-3, 3, by = 1)
# - breaks = c(0.05, 0.1, 0.5, 1, 5, 10, 50)
# - breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10)
# - axis limits: extend to nearest 10-fold? (initial.ratio, fitness, fitness.ratio)
# - make.linear.labels on fitness scale 
# - easily-modified variables for point colors
# - larger points when multiple zero counts? 
# 
# Error warnings: 
# - negative abundance
# - number.X > number.total
# - initial or final proportion not in [0, 1]
# - zeroes/Inf in outcomes
# - more than one strain.A or more than one strain.B
# - stop script execution if warning?
# 
# Calculations
# - format numbers for results-fitness.txt, not formatC()
# 
# Code cleanup and streamlining:
# - revise description and instructions
# - revise output annotation: script version number? more description? 
# - include license (GPL-- GNU General Public License?)
# 

# 
# WHAT THIS SCRIPT DOES
# 
# This script takes data from microbial "mix experiments", calculates some recommended fitness 
# measures, and plots them against the initial abundance of each strain in the mixture.  
# For more details, see smith & Inglis [ref here] 
# 
# 
# HOW TO USE THIS SCRIPT
# 
# 1. Obtain a copy of the R free software environment for statistical computing
#    from https://www.r-project.org/
# 2. If you don't already have the ggplot2 or gtable package, run this command in R: 
#    install.packages(c("ggplot2", "gtable"))
# 3. Format your data to match the script's expected input (described below)
# 4. Change the filename here to your data: 

	data.filename <- "data-tmp.txt"

# 5. Run the script in R by selecting "Edit/Source Document" from the top menu (cmd-E in macOS) 
# 6. Calculated fitness measures will be output to the file "results-fitness.txt". 
#    Fitness plots will be output to files "results-plot1.pdf" and "results-plot2.pdf"
# 
# 
# EXPECTED DATA FORMAT
# 
# The script expects data to be a tab-delimited text file containing a table with a header row. 
# An easy way to generate such a file is to copy and paste cells from MS Excel into a text 
# editor like NotePad or TextEdit. 
# 
# Every row must include columns named "strain.A" and "strain.B" that include the names 
# of the two strains in the mix experiment. Strain A is assumed to be the primary strain of 
# interest. 
# 
# The data table must also include columns with sufficient data to identify the initial and 
# final abundance af each strain. The possible columns are: 
# 
#  initial.number.A, initial.number.B, initial.number.total, 
#  initial.proportion.A, initial.proportion.B, 
#  final.number.A, final.number.B, final.number.total, 
#  final.proportion.A, final.proportion.B
# 
# The data table need not include all of these columns, but it must include two initial and 
# two final quantities from which the rest can be calculated. Ideally, you'd include whatever 
# most closely matches what you actually measured. The columns do not need to be in any 
# particular order. 
# 
# When only one strain is present, mark the other strain's abundance as zero (do not leave 
# empty or mark NA). Density (cells/ml) can be used for number with no difficulty. The script 
# assumes that each row is data from an independent biological replicate, so any repeated 
# measurements from the same replicate must but dealt with (averaged) beforehand. 
# 
# For example, a simple dataset might look like this: 
# 
#  # You can include metadata and other comments like this
#  # strain descriptions, experimental treatments, author, date, etc.
#  # 
#  strain.A	strain.B	initial.number.A	initial.number.B	final.number.A	final.number.B
#  GVB206.3	GJV10	5.00E+08	0.00E+00	6.00E+01	0.00E+00
#  GVB206.3	GJV10	2.50E+08	2.50E+08	1.07E+05	3.90E+05
#  GVB206.3	GJV10	0.00E+00	5.00E+08	0.00E+00	1.76E+07
# 

# =================================================
# YOU DON'T NEED TO CHANGE ANYTHING BELOW THIS LINE
# =================================================

# Tested with R version 3.6.1 (2019-07-05)
library("ggplot2")  # For plots, tested with version 3.2.1
library("gtable")   # For multiple plots in one window, tested with version 0.3.0

analyze.mix.expt <- function(input.filename) {
	# Wrapping everything in a function to keep namespace and memory clean

	# Load data
	dataset <- read.table(input.filename, header = TRUE, sep = "\t", as.is = TRUE)
	dataset.tmp <- dataset  # Copy for calculations

	# 
	# Calculate fitness
	# 

	dataset.tmp <- within(dataset.tmp, {
	
		# Calculate initial population state
		if (exists("initial.number.A") & exists("initial.number.B")) {
			initial.number.total <- initial.number.A + initial.number.B
			initial.proportion.A <- initial.number.A / initial.number.total
			initial.proportion.B <- initial.number.B / initial.number.total
			initial.ratio.A      <- initial.number.A / initial.number.B
		} else if (exists("initial.number.total") & exists("initial.proportion.A")) {
			initial.proportion.B <- 1 - initial.proportion.A
			initial.number.A     <- initial.number.total * initial.proportion.A
			initial.number.B     <- initial.number.total * initial.proportion.B
			initial.ratio.A      <- initial.number.A / initial.number.B
		} else if (exists("initial.number.total") & exists("initial.proportion.B")) {
			initial.proportion.A <- 1 - initial.proportion.B
			initial.number.A     <- initial.number.total * initial.proportion.A
			initial.number.B     <- initial.number.total * initial.proportion.B
			initial.ratio.A      <- initial.number.A / initial.number.B
		} else if (exists("initial.number.total") & exists("initial.number.A")) {
			initial.number.B     <- initial.number.total - initial.number.A
			initial.proportion.A <- initial.number.A / initial.number.total
			initial.proportion.B <- initial.number.B / initial.number.total
			initial.ratio.A      <- initial.number.A / initial.number.B
		} else if (exists("initial.number.total") & exists("initial.number.B")) {
			initial.number.A     <- initial.number.total - initial.number.B
			initial.proportion.A <- initial.number.A / initial.number.total
			initial.proportion.B <- initial.number.B / initial.number.total
			initial.ratio.A      <- initial.number.A / initial.number.B
		} else {
			warning("Cannot determine initial population state from data")
		}

		# Calculate final population state
		if (exists("final.number.A") & exists("final.number.B")) {
			final.number.total <- final.number.A + final.number.B
			final.proportion.A <- final.number.A / final.number.total
			final.proportion.B <- final.number.B / final.number.total
		} else if (exists("final.number.total") & exists("final.proportion.A")) {
			final.proportion.B <- 1 - final.proportion.A
			final.number.A     <- final.number.total * final.proportion.A
			final.number.B     <- final.number.total * final.proportion.B
		} else if (exists("final.number.total") & exists("final.proportion.B")) {
			final.proportion.A <- 1 - final.proportion.B
			final.number.A     <- final.number.total * final.proportion.A
			final.number.B     <- final.number.total * final.proportion.B
		} else if (exists("final.number.total") & exists("final.number.A")) {
			final.number.B     <- final.number.total - final.number.A
			final.proportion.A <- final.number.A / final.number.total
			final.proportion.B <- final.number.B / final.number.total
		} else if (exists("final.number.total") & exists("final.number.B")) {
			final.number.A     <- final.number.total - final.number.B
			final.proportion.A <- final.number.A / final.number.total
			final.proportion.B <- final.number.B / final.number.total
		} else {
			warning("Cannot determine final population state from data")
		}

		# Calculate fitness measures
		fitness.A       <- final.number.A / initial.number.A
		fitness.B       <- final.number.B / initial.number.B
		fitness.total   <- final.number.total / initial.number.total
		fitness.ratio.A <- fitness.A / fitness.B

		# Format results (scientific notation)
		# initial.ratio.A <- formatC(initial.ratio.A, format = "e", digits = 6)
		# fitness.A       <- formatC(fitness.A, format = "e", digits = 6)
		# fitness.B       <- formatC(fitness.B, format = "e", digits = 6)
		# fitness.total   <- formatC(fitness.total, format = "e", digits = 6)
		# fitness.ratio.A <- formatC(fitness.ratio.A, format = "e", digits = 6)

		# NA results for single-strain populations
		fitness.A[initial.number.A == 0] <- NA
		fitness.B[initial.number.B == 0] <- NA
		fitness.ratio.A[(initial.number.A == 0) | (initial.number.B == 0)] <- NA

	})

	# Add results to input data frame
	dataset$initial.proportion.A <- dataset.tmp$initial.proportion.A
	dataset$initial.ratio.A      <- dataset.tmp$initial.ratio.A
	dataset$fitness.A       <- dataset.tmp$fitness.A
	dataset$fitness.B       <- dataset.tmp$fitness.B
	dataset$fitness.total   <- dataset.tmp$fitness.total
	dataset$fitness.ratio.A <- dataset.tmp$fitness.ratio.A

	# Write results to file
	my.connection <- file("results-fitness.txt", open = "wt")
	writeLines(
		c(
			"# Calculated fitness values for microbial mix experiment", 
			paste("# Created:", Sys.time(), "by analyze-mix-experiment.R"), 
			paste("# Input data:", input.filename), 
			"# ", 
			"# Wrightian fitness w = final.number / initial.number", 
			"# Fitness ratio = w_A / w_B", 
			"# "
		), my.connection
	)
	suppressWarnings(write.table(
		dataset, my.connection, quote = FALSE, sep = "\t", row.names = FALSE, append = TRUE
	))
	close(my.connection)

	# 
	# Make plots
	# 

	# Format data for plotting
	name.strain.A <- unique(dataset$strain.A)[1]
	name.strain.B <- unique(dataset$strain.B)[1]
	data.for.plot.long <- reshape(
		subset(dataset, select = c(
			initial.proportion.A, initial.ratio.A, fitness.A, fitness.B, fitness.total
		)), 
		direction = "long", 
		varying = c("fitness.A", "fitness.B", "fitness.total"), v.names = c("fitness"), 
		times = c(name.strain.A, name.strain.B, "Total group"), timevar = "strain"
	)
	data.for.plot.long <- within(data.for.plot.long, {
		strain <- factor(strain, levels = c(name.strain.A, name.strain.B, "Total group"))
		my.facet <- (strain == "Total group")
	})
	data.for.plot.long <- subset(data.for.plot.long, !is.na(fitness))
	data.for.plot.wide <- subset(dataset, !is.na(fitness.ratio.A), 
		select = c(initial.proportion.A, initial.ratio.A, fitness.ratio.A)
	)

	# Define fitness axis limits
	limits.fitness <- with(data.for.plot.long, 
		range(1, fitness[fitness > 0], na.rm = TRUE)
	)
	limits.fitness.ratio <- with(data.for.plot.wide, 
		range(1, fitness.ratio.A[(fitness.ratio.A > 0) & is.finite(fitness.ratio.A)], na.rm = TRUE)
	)
	range.fitness       <- log10(limits.fitness[2]/limits.fitness[1])
	range.fitness.ratio <- log10(limits.fitness.ratio[2]/limits.fitness.ratio[1])

	# Minimum 10-fold fitness range
	if (range.fitness < 1) {
		limits.fitness <- c(
			10^(log10(limits.fitness[1]) - (1 - range.fitness)/2), 
			10^(log10(limits.fitness[2]) + (1 - range.fitness)/2)
		)
		range.fitness <- log10(limits.fitness[2]/limits.fitness[1])
	}
	if (range.fitness.ratio < 1) {
		limits.fitness.ratio <- c(
			10^(log10(limits.fitness.ratio[1]) - (1 - range.fitness.ratio)/2), 
			10^(log10(limits.fitness.ratio[2]) + (1 - range.fitness.ratio)/2)
		)
		range.fitness.ratio <- log10(limits.fitness.ratio[2]/limits.fitness.ratio[1])
	}

	# Shared fitness scale
	if (range.fitness > range.fitness.ratio) {
		limits.fitness.ratio <- c(
			10^(log10(limits.fitness.ratio[1]) - (range.fitness - range.fitness.ratio)/2), 
			10^(log10(limits.fitness.ratio[2]) + (range.fitness - range.fitness.ratio)/2) 
		)
		range.fitness.ratio <- log10(limits.fitness.ratio[2]/limits.fitness.ratio[1])
	} else if (range.fitness.ratio > range.fitness) {
		limits.fitness<- c(
			10^(log10(limits.fitness[1]) - (range.fitness.ratio - range.fitness)/2), 
			10^(log10(limits.fitness[2]) + (range.fitness.ratio - range.fitness)/2) 
		)
		range.fitness <- log10(limits.fitness[2]/limits.fitness[1])
	}
	# print(c(log10(limits.fitness), range.fitness))
	# print(c(log10(limits.fitness.ratio), range.fitness.ratio))

	# Define fitness breaks and labels
	make.log.labels <- function(x) {
		if (x == 1) { 
			as.expression(1) 
		} else {
			as.expression(bquote(10^.(log10(x))))
		}
	}
	make.linear.labels <- function(x) {
		if (x == 1) { 
			as.expression(1) 
		} else if (x == 0) { 
			as.expression(0)
		} else { 
			as.expression(x)
		}
	}
	if (range.fitness > 4) {
		breaks.fitness       <- 10^seq(-10, 10, by = 2)
		breaks.fitness.minor <- 10^seq(-10, 10, by = 1)
		labels.fitness       <- sapply(breaks.fitness, make.log.labels)
	} else if (range.fitness > 2) {
		breaks.fitness       <- 10^c(-10:10)
		breaks.fitness.minor <- 10^c(-10:10)
		labels.fitness       <- sapply(breaks.fitness, make.log.labels)
	} else if (range.fitness > 1) {
		breaks.fitness       <- c(0.05, 0.1, 0.5, 1, 5, 10, 50)
		breaks.fitness.minor <- c(0.01*1:9, 0.1*1:9, 1:9, 10*1:9)
		labels.fitness       <- c(0.05, 0.1, 0.5, 1, 5, 10, 50)
	} else {
		breaks.fitness       <- c(0.1, 0.2, 0.5, 1, 2, 5, 10)
		breaks.fitness.minor <- c(0.1*1:9, 1:9)
		labels.fitness       <- c(0.1, 0.2, 0.5, 1, 2, 5, 10)
	}

	# Define initial ratio axis: limits, breaks, labels
	limits.initial.ratio <- range(
		1e-1, 1e1,  # Minimum range
		with(data.for.plot.long, 
			initial.ratio.A[(initial.proportion.A > 0) & (initial.proportion.A < 1)]
		), 
		na.rm = TRUE
	)
	breaks.initial.ratio <- 10^c(-10:10)
	if (log10(limits.initial.ratio[2]/limits.initial.ratio[1]) < 3) {
		breaks.initial.ratio.minor <- 5 * 10^c(-10:10)
	} else {
		breaks.initial.ratio.minor <- 10^c(-10:10)
	}
	labels.initial.ratio <- sapply(breaks.initial.ratio, make.log.labels)
	labels.initial.proportion <- sapply(seq(0, 1, by = 0.2), make.linear.labels)

	# Define shared plot elements
	plot.base <- ggplot() + 
		scale_fill_manual(values  = c("tan",  "lightsteelblue",  gray(0.65))) + 
		scale_color_manual(values = c("tan4", "lightsteelblue4", gray(0.1))) + 
		geom_hline(yintercept = 1, color = "white", size = 0.8) + 
		ggtitle("") +
		theme(
			text                 = element_text(size = 8), 
			legend.title         = element_blank(), 
			legend.background    = element_blank(), 
			legend.direction     = "horizontal", 
			legend.justification = c(0.5, 0.15), 
			legend.position      = c(0.5, 1),
			strip.text           = element_blank(),
			strip.background     = element_blank()
		)
	scale.initial.proportion <- scale_x_continuous(
		name         = paste("Initial proportion", name.strain.A), 
		limits       = c(0, 1), 
		breaks       = seq(0, 1, by = 0.2), 
		minor_breaks = seq(0, 1, by = 0.1), 
		labels       = labels.initial.proportion
	)
	scale.initial.ratio <- scale_x_log10(
		name         = bquote("Initial ratio" ~ .(name.strain.A) / .(name.strain.B)), 
		limits       = limits.initial.ratio, 
		breaks       = breaks.initial.ratio, 
		minor_breaks = breaks.initial.ratio.minor,
		labels       = labels.initial.ratio
	)
	scale.fitness <- scale_y_log10(
		name         = "Wrightian fitness\n (final no. / initial no.)", 
		limits       = limits.fitness, 
		breaks       = breaks.fitness,
		minor_breaks = breaks.fitness.minor, 
		labels       = labels.fitness
	)
	scale.fitness.ratio <- scale_y_log10(
		name         = bquote("Fitness ratio" ~ .(name.strain.A)/.(name.strain.B)),
		limits       = limits.fitness.ratio, 
		breaks       = breaks.fitness,
		minor_breaks = breaks.fitness.minor, 
		labels       = labels.fitness
	)

	# Create individual plots
	plot.fitness <- plot.base %+% data.for.plot.long + 
		aes(x = initial.proportion.A) + scale.initial.proportion + 
		aes(y = fitness) + scale.fitness + 
		aes(color = strain, fill = strain) + facet_wrap(~ my.facet, nrow = 1) + 
		geom_point(shape = 21) + geom_point(shape = 1)
	plot.fitness.logratio <- plot.base %+% data.for.plot.long + 
		aes(x = initial.ratio.A) + scale.initial.ratio + 
		aes(y = fitness) + scale.fitness + 
		aes(color = strain, fill = strain) + facet_wrap(~ my.facet, nrow = 1) + 
		geom_point(shape = 21) + geom_point(shape = 1)
	plot.fitness.ratio <-  plot.base %+% data.for.plot.wide + 
		aes(x = initial.proportion.A) + scale.initial.proportion + 
		aes(y = fitness.ratio.A) + scale.fitness.ratio + 
		geom_point(color = gray(0.65)) + geom_point(shape = 1)
	plot.fitness.ratio.logratio <-  plot.base %+% data.for.plot.wide + 
		aes(x = initial.ratio.A) + scale.initial.ratio + 
		aes(y = fitness.ratio.A) + scale.fitness.ratio + 
		geom_point(color = gray(0.65)) + geom_point(shape = 1)

	# Combine plots
	pdf(file = NULL)  
		# workaround so that ggplotGrob() doesn't create a blank plot window
		# see https://github.com/tidyverse/ggplot2/issues/809
	plot.proportion <- gtable_add_grob(
		gtable(widths = unit(rep(1, 25), "null"), heights = unit(rep(1, 1), "null")), 
		list(ggplotGrob(plot.fitness), ggplotGrob(plot.fitness.ratio)), 
		l = c(1, 17), r = c(16, 25), t = c(1, 1), b = c(1, 1)
	)
	plot.logratio <- gtable_add_grob(
		gtable(widths = unit(rep(1, 25), "null"), heights = unit(rep(1, 1), "null")), 
		list(ggplotGrob(plot.fitness.logratio), ggplotGrob(plot.fitness.ratio.logratio)), 
		l = c(1, 17), r = c(16, 25), t = c(1, 1), b = c(1, 1)
	)
	dev.off()  # end of workaround

	# Save plots to pdf
	ggsave("results-plot1.pdf", plot.proportion, device = "pdf", width = 6.25, height = 2.25)
	ggsave("results-plot2.pdf", plot.logratio,   device = "pdf", width = 6.25, height = 2.25)
}
analyze.mix.expt(data.filename)
# tmp <- analyze.mix.expt(data.filename)

