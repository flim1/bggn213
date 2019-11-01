#' ---
#' title: "Class 5 Data Exploration and visualization in R"
#' output: github_document
#' ---

# Class 5 Data visualization

x <- rnorm(1000)

# How many things are in x
length(x)

mean(x)
sd(x)

summary(x)

boxplot(x)
hist(x)
rug(x)

hist(x, breaks = 3)
hist(x, breaks = 30)

# Section 2
weight <- read.table(file="bimm143_05_rstats/weight_chart.txt", 
                     header = TRUE)

plot(weight$Age, weight$Weight, 
     type="o", 
     pch=15, 
     cex=1.5, 
     lwd=2, 
     ylim=c(2,10), 
     xlab="Age (months)", 
     ylab="Weight (kg)", 
     main="Baby weight with age")

mouse <- read.table(file="bimm143_05_rstats/feature_counts.txt", 
           header=TRUE,
           sep="\t")

par(mar=c(3,9,3,5))
barplot(mouse$Count,
        horiz=TRUE,
        names.arg=mouse$Feature,
        main="Number of features in the mouse GRCm38 genome",
        las=1,
        xlim=c(0,80000),
        cex.axis=0.7,
        cex.name=0.6)


# Section 3. Using color
mf <- read.delim("bimm143_05_rstats/male_female_counts.txt")

par(mar=c(6,6,2,2))
barplot(mf$Count,
        names.arg = mf$Sample,
        ylab="Counts",
        cex.names=0.7,
        las=2,
        col=c("blue2","red2"))

genes <- read.delim("bimm143_05_rstats/up_down_expression.txt")
nrow(genes)
table(genes$State)

plot(genes$Condition1, genes$Condition2,
     xlab="Expression condition 2",
     ylab="Expression condition 1",
     col=genes$State)


