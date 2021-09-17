## MSB7102 Mini-project, Semester I, 2021


**1.Import the data described above into R, provide descriptive summaries of the subject data (using appropriate graphics and statistical summary measures) given in the diabimmune_16s_t1d_metadata.csv file. In addition, use appropriate test(s) to check for association/independency between disease status and other variables (delivery mode, gender and age).**
  
  ```{r}
data <- read.csv("taxa_table_.csv");
#View(data)
```

```{r}
summary(data)
```

```{r}
head(data)
```

```{r}
dim(data)
```

```{r}
class(data)
```


**Find the number of cases and controls, males and females and the delivery routes.**
  
  ```{r}
table(data$Case_Control)
table(data$Gender)
table(data$Delivery_Route)
```


**Plot bargraphs showing the different associations.**
  **Load the ggplot2 package.**
  
  ```{r}
library(ggplot2)
```

```{r}
qplot(data$Case_Control, fill = data$Gender) + geom_bar() + labs(title = "A bargraph showing the Cases and Controls against Gender", x = "Case_Control", y= "Frequency", fill = "Gender")
```


```{r}
qplot(data$Case_Control, fill = data$Delivery_Route) + geom_bar()+ labs(title = "A bargraph showing the Cases and Controls against Delivery Route", x = "Case_Control", y= "Frequency", fill = "Delivey Route")
```


```{r}
qplot(data$Age_at_Collection) + geom_bar() + labs(title = "A bargraph showing the Age of participants", x = "Age (Days)", y= "Frequency")
```


**Check for association/independency between disease status and age, gender and delivery mode.**
  
  **Generate a contingency table and perform a Chi-square statistics.**
  
  ```{r}
tbAge <- table(data$Case_Control,data$Age_at_Collection)
chisq.test(tbAge)
```

**From the statistical analysis, the P value is 0.6115, which is greater than 0.05, therefore we accept the null hypothesis indicating that there is no association between the disease state and the ages of the participants at collection.**
  
  
  ```{r}
tbGender <- table(data$Case_Control, data$Gender) 
chisq.test(tbGender)
```

**From the statistical analysis, the P value is 0.5796, which is greater than 0.05, therefore we accept the null hypothesis indicating that there is no association between the disease state and the gender.**
  
  
  ```{r}
tbDelivery_Route <- table(data$Case_Control, data$Delivery_Route)
chisq.test(tbDelivery_Route)
```

**From the statistical analysis, the P value is 3.949e-09, which is less than 0.05, therefore we accept the alternate hypothesis indicating that there is a strong association between the disease state and the delivery route.**
  
  
  
  **2.	Using phyloseq, create a phyloseq object. This will comprise the OTU abundance, taxonomy (provided in the .txt file) and sample data (provided in the .csv file).**
  
  **Import the OTU table.**
  ```{r}
otuTable <- read.table("otu_table")  #Importing the OTU table
head(otuTable, n=1)
```


**Import the taxonomy table.**
  ```{r}
taxaTable <- read.table("taxa_table") #Importing the Taxa table
head(taxaTable, n=1)
```


**Check the dimension and class of each table.**
  
  ```{r}
dim(otuTable)
```

```{r}
dim(taxaTable)
```

```{r}
class(otuTable)
```
```{r}
class(taxaTable)
```


**Convert the tables from dataframe to matrix. This is because, phyloseq works better with a matrix.**
  
  ```{r}
#Converting the Taxa and OTU Table into a Matrix
mtaxaTable <- as.matrix(taxaTable)
motuTable <- as.matrix(otuTable)
```


**Check the class of the matrices created.**
  
  ```{r}
class(motuTable)
```
```{r}
class(mtaxaTable)
```


**Data Cleansing.**
  
  This is done to have consistent data across all the matrices. It involves making sure that the OTU/taxa row names match. Currently they don't as taxa have a trailing ";"

#rownames(mtaxaTable)[rownames(mtaxaTable) == "4333897;"] = "4333897"

```{r}
head(mtaxaTable, n=1)
Rnames <- rownames(mtaxaTable)  #Extract rownames from the matrix
NRnames <- gsub(x = Rnames, pattern = ";", replacement = "")  #Remove the ; from the extracted rownames
#NRnames
rownames(mtaxaTable) <- NRnames  #Set the new rownames
head(mtaxaTable, n=1) #Check to confirm that changes have been made

```


**Load the phyloseq package.**
```{r}
library(phyloseq)
```


**Create a phyloseq object**
```{r}
#Tell phyloseq to load them into a phyloseq object
OTU = otu_table(motuTable, taxa_are_rows = TRUE)
TAX = tax_table(mtaxaTable)

#OTU
#TAX

#Generating the phyloseq object and view it
physeq = phyloseq(OTU, TAX)
physeq

#Plotting the phyloseq
plot_bar(physeq, fill = "Family.")
```



**3. Generate Alpha diversity plots and ordination plots. Examine any observed patterns by delivery mode, gender and disease status.**


**Generate a default plot from the plot_richness function.**

```{r}
plot_richness(physeq)  #Default plot produced by the plot_richness function
```


**Merge the sample data into the phyloseq object to observe patterns across variables.**

```{r}
library(tibble)
s_data <- column_to_rownames(data, var="ï..Sample_ID")
s_data <- sample_data(s_data)
#s_data
mergedPhyseq <- merge_phyloseq(physeq, s_data)
mergedPhyseq
```


**Alpha diversity comparison between the Delivery route in cases and controls.**
```{r}
plot_richness(mergedPhyseq, x ="Case_Control", color="Delivery_Route")
```


**Alpha diversity comparison between the Gender in Cases and Controls.**

```{r}
plot_richness(mergedPhyseq, x ="Case_Control", color="Gender")

```


**Alpha diversity comparison between the Age in Cases and Controls.**

```{r}
plot_richness(mergedPhyseq, x ="Case_Control", color="Age_at_Collection")
```



**4. Perform differential abundance using DEseq2.**


**Load the DESeq2 package.**

```{r}
library("DESeq2")
```


**Convert the phyloseq object to deseq.**

```{r}
cacos = phyloseq_to_deseq2(mergedPhyseq, ~ Case_Control)
```


**Calculate geometric mean to estimate size factor.**

```{r}
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(cacos), 1, gm_mean)
cacos = estimateSizeFactors(cacos, geoMeans = geoMeans)
cacos = DESeq(cacos, fitType="local")
```


**Perform the test using the DESeq function.**

```{r}
cacos = DESeq(cacos, test="Wald", fitType="parametric")
```


**Investigate the test result table.**

```{r}
res = results(cacos, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(mergedPhyseq)[rownames(sigtab), ], "matrix"))
head(sigtab)
```

```{r}
dim(sigtab)
```


**Visualize the result.**

```{r}
library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
```

