# analyze-mix-experiments.R v0.3
# jeff smith & Fredrick R. Inglis, 2015-2016

# TO DO LIST: 
# 
# - Make 10^0 labels "1"
# - Extend limits to nearest major tick?  (If < half segment needed)
# - More space between top and bottom rows?
# - BUG: Sometimes wierd shit happening with first plot.  
#        Has to do with dev.new() being called in block, I think. 
# 
# - Stop script execution if warning("Cannot determine final strain abundances from data")
# - Annotate results.txt: explain Wrightian fitness
# - HOW TO: Explain how script keeps "extra columns"
# - "Number" vs "density:" what do?
# 


# 
# WHAT THIS SCRIPT DOES
# 
# This script takes data from microbial "mix experiments", calculates a variety of fitness 
# measures, and plots them against the initial abundance of each strain in the mixture.  
# For more details, see [full reference here]. 
# 
# 
# HOW TO USE THIS SCRIPT
# 
# 1. Obtain a copy of the R free software environment for statistical computing and 
#    graphics (http://www.r-project.org/)
# 2. If you don't have the ggplot2 and gtable packages, run this command in R: 
#    install.packages(c("ggplot2", "gtable"))
# 3. Format your data to match the script's expected input (described below)
# 4. Change the filename here to your data: 
# 
     # data.filename <- "example-data.txt"
     data.filename <- "data-Fiegna-2006-OC-PX.txt"
# 
# 5. Run the script in R by selecting "Edit/Source Document" from the top menu (cmd-E in OS X) 
# 
# 
# 
# EXPECTED DATA FORMAT
# 
# The script expects data to be a tab-delimited text file containing a table with a header 
# row. An easy way to generate such a file is to copy and paste cells from MS Excel into a 
# text editor like NotePad or TextEdit. 
# 
# Every row must include columns named "strain.A" and "strain.B" that include the names 
# of the two strains in the mix experiment. Strain A is assumed to be the primary strain of 
# interest. 
# 
# The data table must also include columns with sufficient data to identify the initial and 
# final abundance af each strain. The possible columns are: 
# 
# 	initial.number.A, initial.number.B, final.number.A, final.number.B, 
# 	initial.number.total, initial.proportion.A, initial.proportion.B, 
# 	final.number.total, final.proportion.A, and final.proportion.B
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
# For example, a simple data set might look like this: 
# 
#   strain.A	strain.B	initial.number.A	initial.number.B	final.number.A	final.number.B
#   GVB206.3	GJV10	5.00E+08	0.00E+00	6.00E+01	0.00E+00
#   GVB206.3	GJV10	2.50E+08	2.50E+08	1.07E+05	3.90E+05
#   GVB206.3	GJV10	0.00E+00	5.00E+08	0.00E+00	1.76E+07
# 
# 
# 
# LICENSE
# 
# This work is licensed under a Creative Commons Attribution 4.0 International License
# (https://creativecommons.org/licenses/by/4.0/). 
# 
# If you use this script (or derivatives of this script), please cite: 
# [full reference here]
# 


# =================================================
# You don't need to change anything below this line
# =================================================

## 
## LOAD PACKAGES AND DATA
## 

library("ggplot2")  # for graphics, tested with v2.0.0
library("gtable")   # for multiple plots in single window, tested with v0.1.2
# library("scales")    # for better log-scale tick labels, tested with v0.2.4
# Script tested with R v3.2.3

my.data <- read.table(data.filename, header = TRUE, sep = "\t", as.is = TRUE)


## 
## FORMAT DATA AND CALCULATE FITNESS MEASURES
## 

# Extract names of strains in experiment
name.strain.A <- unique(my.data$strain.A)[1]
name.strain.B <- unique(my.data$strain.B)[1]

# Calculate initial strain abundances from what's given
my.data <- within(my.data, {
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
		warning("Cannot determine initial strain abundances from data")
	}
})

# Remove Inf/0 values from initial ratio 
# my.data$initial.ratio.A[is.infinite(my.data$initial.ratio.A)] <- NA
# my.data$initial.ratio.A[my.data$initial.ratio.A == 0]         <- NA

# Calculate final strain abundances from what's given
my.data <- within(my.data, {
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
		warning("Cannot determine final strain abundances from data")
	}
})

# Calculate fitness measures
my.data <- within(my.data, {
	fitness.A       <- final.number.A / initial.number.A
	fitness.B       <- final.number.B / initial.number.B
	fitness.total   <- final.number.total / initial.number.total
	fitness.ratio.A <- fitness.A / fitness.B
	fitness.ratio.B <- fitness.B / fitness.A
	# initial.ratio.1p <- (initial.number.A + 1) / (initial.number.B + 1)
})

# Remove Inf/0 values from fitness ratio 
my.data$fitness.ratio.A[is.infinite(my.data$fitness.ratio.A)] <- NA
my.data$fitness.ratio.B[is.infinite(my.data$fitness.ratio.B)] <- NA
my.data$fitness.ratio.A[my.data$fitness.ratio.A == 0]         <- NA
my.data$fitness.ratio.B[my.data$fitness.ratio.B == 0]         <- NA
my.data$fitness.ratio.A[is.nan(my.data$fitness.ratio.A)]      <- NA
my.data$fitness.ratio.B[is.nan(my.data$fitness.ratio.B)]      <- NA

# head(my.data)

# Write calculated fitness measures to file
my.connection <- file("results.txt", open = "wt")
writeLines(paste("# Test annotation"), my.connection)
writeLines(paste("#", data.filename), my.connection)
writeLines(paste("#", Sys.time()), my.connection)
suppressWarnings(write.table(
	my.data, my.connection, 
	quote = FALSE, sep = "\t", row.names = FALSE, append = TRUE
))
close(my.connection)

# Make long-form data frame for plotting with color scale for strain
tmp.data.A <- subset(my.data, select = c(
	strain.A, 
	initial.number.A, final.number.A, fitness.A, 
	initial.proportion.A, initial.ratio.A 
	#, initial.ratio.1p
))
tmp.data.B <- subset(my.data, select = c(
	strain.B, 
	initial.number.B, final.number.B, fitness.B, 
	initial.proportion.A, initial.ratio.A
	# , initial.ratio.1p
))
tmp.data.total <- subset(my.data, select = c(
	strain.A, 
	initial.number.total, final.number.total, fitness.total, 
	initial.proportion.A, initial.ratio.A
	# , initial.ratio.1p
))
colnames(tmp.data.A) <- c(
	"strain", "initial.number", "final.number", "fitness", 
	"initial.proportion.A", "initial.ratio.A"
	# , "initial.ratio.1p"
)
colnames(tmp.data.B)     <- colnames(tmp.data.A)
colnames(tmp.data.total) <- colnames(tmp.data.A)
tmp.data.total$strain    <- rep("Total group", nrow(tmp.data.total))
my.data.long <- rbind(tmp.data.A, tmp.data.B, tmp.data.total)
my.data.long <- subset(my.data.long, (initial.number > 0) & !(is.na(final.number)) )
my.data.long$strain <- factor(my.data.long$strain, levels = c(name.strain.A, name.strain.B, "Total group"))

# head(my.data.long); tail(my.data.long)


##
## DETERMINE BEST SCALES FOR PLOTS
##

# Fitness y-axis

fitness.limits <- range(
	1, 
	my.data.long$fitness[my.data.long$fitness > 0], # Includes total group
	na.rm = TRUE
)
fitness.range.log10 <- log10(fitness.limits[2]/fitness.limits[1])

if ( fitness.range.log10 < 1 ) {	
	# Minimum 10-fold range
	fitness.midpoint    <- log10(fitness.limits[1]) + fitness.range.log10/2
	fitness.limits[1]   <- 10^(fitness.midpoint - 0.5)
	fitness.limits[2]   <- 10^(fitness.midpoint + 0.5)
	fitness.range.log10 <- 1
}

if (fitness.range.log10 > 4) {
	fitness.breaks       <- 10^seq(-16, 16, by = 2)
	fitness.breaks.minor <- 10^seq(-16, 16, by = 1)
	fitness.labels       <- sapply(
		fitness.breaks, function(i) as.expression(bquote( 10^ .(log10(i)) ))
	)
} else if (fitness.range.log10 > 1.5) {
	fitness.breaks       <- 10^c(-5:5)
	fitness.breaks.minor <- 0.5 * 10^c(-5:5)
	fitness.labels       <- sapply(
		fitness.breaks, function(i) as.expression(bquote( 10^ .(log10(i)) ))
	)
} else {
	fitness.breaks       <- c(0.1, 0.2, 0.5, 1, 2, 5, 10)
	fitness.breaks.minor <- c(0.1 * 1:9, 1:9)
	fitness.labels       <- c(0.1, 0.2, 0.5, 1, 2, 5, 10)
}

# Fitness ratio y-axis
limits.fitness.ratio      <- range(1, my.data$fitness.ratio.A, na.rm = TRUE)
fitness.ratio.range.log10 <- log10(limits.fitness.ratio[2]/limits.fitness.ratio[1])
midpoint.fitness.ratio    <- log10(limits.fitness.ratio[1]) + fitness.ratio.range.log10/2
limits.fitness.ratio[1]   <- 10^(midpoint.fitness.ratio - fitness.range.log10/2) # Use fitness scale
limits.fitness.ratio[2]   <- 10^(midpoint.fitness.ratio + fitness.range.log10/2) # 

# Initial ratio x-axis
initial.ratio.limits <- range(
	1e-1, 1e1,
	my.data.long$initial.ratio.A[
		(my.data.long$initial.ratio.A > 0) & !(is.infinite(my.data.long$initial.ratio.A))
	], 
	na.rm = TRUE
)
initial.ratio.breaks <- 10^c(-8:8)
initial.ratio.labels <- sapply(
	initial.ratio.breaks, function(i) as.expression(bquote( 10^ .(log10(i)) ))
)
initial.ratio.range.log10 <- log10(initial.ratio.limits[2]/initial.ratio.limits[1])
if (initial.ratio.range.log10 < 3) {
	initial.ratio.breaks.minor <- 5 * 10^c(-4:4)
} else {
	initial.ratio.breaks.minor <- 10^c(-8:8)
}



## 
## MAKE SHARED PLOT ELEMENTS
## 

# Make base (empty) plot
base.plot <- ggplot(data = my.data) + 
	scale_fill_manual( values = c("tan",  "lightsteelblue", gray(0.65))) + 
	scale_color_manual(values = c("tan4", "lightsteelblue4", gray(0.1))) + 
	geom_hline(yintercept = 1, color = "white", size = 1.2) + 
	theme(
		text         = element_text(size = 8), 
		# plot.title   = element_text(hjust = -0.12, vjust = 1, face = "bold"),  # Panel labels
		# axis.title.x = element_text(vjust = 0),     # More space betweeen axis labels and tick labels
		# axis.title.y = element_text(vjust = 0.75), 
		legend.title         = element_blank(), 
		legend.background    = element_blank(), 
		legend.direction     = "horizontal", 
		legend.justification = c(0.5, 0.15), 
		legend.position      = c(0.5, 1),
		strip.text           = element_blank(),
		strip.background     = element_blank()
	)  

# Make shared axes
scale.initial.proportion <- scale_x_continuous(
	name         = paste("Initial proportion", name.strain.A), 
	limits       = c(0, 1), 
	breaks       = seq(0, 1, by = 0.2), 
	minor_breaks = seq(0, 1, by = 0.1)
)
scale.initial.ratio <- scale_x_log10(
	name   = bquote("Initial ratio" ~ .(name.strain.A) / .(name.strain.B)), 
	limits = initial.ratio.limits, 
	breaks = initial.ratio.breaks, 
	minor_breaks = initial.ratio.breaks.minor,
	labels = initial.ratio.labels
	# breaks = trans_breaks("log10", function(x) 10^x),
	# labels = trans_format("log10", math_format(10^.x))
)
# scale.initial.ratio.1p <- scale_x_log10(
	# name   = bquote("Initial ratio" ~ .(name.strain.A) / .(name.strain.B)), 
	# # limits = initial.ratio.limits, 
	# breaks = initial.ratio.breaks, 
	# minor_breaks = initial.ratio.breaks.minor,
	# labels = initial.ratio.labels
# )
scale.fitness <- scale_y_log10(
	name         = "Fitness\n (final num / initial num)", 
	limits       = fitness.limits, 
	breaks       = fitness.breaks,
	minor_breaks = fitness.breaks.minor, 
	labels       = fitness.labels
)
scale.fitness.ratio <- scale_y_log10(
	name         = bquote("Fitness ratio" ~ .(name.strain.A)/.(name.strain.B)),
	limits       = limits.fitness.ratio,
	breaks       = fitness.breaks, 
	minor_breaks = fitness.breaks.minor, 
	labels       = fitness.labels
)


## 
## MAKE INDIVIDUAL PLOTS
## 

# Plot fitness
data.for.plot <- my.data.long
data.for.plot <- within(data.for.plot, {
	my.facet <- !(strain %in% c(name.strain.A, name.strain.B))
})
plot.fitness <- base.plot %+% data.for.plot + 
	aes(y = fitness) + scale.fitness + 
	geom_point(shape = 21) + geom_point(shape = 1) +
	aes(color = strain, fill = strain) + facet_wrap(~ my.facet, nrow = 1)
plot.fitness.logit <- plot.fitness +
	aes(x = initial.ratio.A) + scale.initial.ratio
# plot.fitness.log1p <- plot.fitness +
	# aes(x = initial.ratio.1p) + scale.initial.ratio.1p
plot.fitness <- plot.fitness + 
	aes(x = initial.proportion.A) + scale.initial.proportion 

# Plot within-group fitness ratio
plot.fitness.ratio <- base.plot %+% 
	subset(my.data, (initial.proportion.A > 0) & (initial.proportion.B > 0)) +
	aes(x = initial.proportion.A, y = fitness.ratio.A) + 
	geom_point(color = gray(0.65)) + geom_point(shape = 1) + 
	scale.initial.proportion + scale.fitness.ratio
plot.fitness.ratio.logit <- base.plot %+% 
	subset(my.data, (initial.proportion.A > 0) & (initial.proportion.B > 0)) +
	aes(x = initial.ratio.A, y = fitness.ratio.A) + 
	geom_point(color = gray(0.65)) + geom_point(shape = 1) + 
	scale.initial.ratio + scale.fitness.ratio
# plot.fitness.ratio.log1p <- base.plot %+% 
	# subset(my.data, (initial.proportion.A > 0) & (initial.proportion.B > 0)) +
	# aes(x = initial.ratio.1p, y = fitness.ratio.A) + 
	# geom_point(color = gray(0.5)) + geom_point(shape = 1) + 
	# scale.initial.ratio.1p + scale.fitness.ratio


## 
## DRAW PAGE OF PLOTS
##

dev.new(width = 6.5, height = 4.5)
my.plots <- gtable(widths = unit(rep(1, 3), "null"), heights = unit(rep(1, 2), "null"))
my.plot.list <- list(
	ggplotGrob(plot.fitness + ggtitle(" ")), 
	ggplotGrob(plot.fitness.ratio + ggtitle(" ")), 
	ggplotGrob(plot.fitness.logit + ggtitle(" ")), 
	ggplotGrob(plot.fitness.ratio.logit + ggtitle(" "))
)
my.plots <- gtable_add_grob(
	my.plots, my.plot.list, 
	l = c(1, 3, 1, 3),  # left extents
	r = c(2, 3, 2, 3),  # right extents
	t = c(1, 1, 2, 2),  # top extents
	b = c(1, 1, 2, 2)   # bottom extents
)
plot(my.plots)

pdf(file = "results.pdf", width = 6.5, height = 4.5)
plot(my.plots)
dev.off()
