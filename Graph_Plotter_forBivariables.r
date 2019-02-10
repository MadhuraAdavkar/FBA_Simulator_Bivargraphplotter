##################################################################################################
#@Auther: Madhura Adavkar
#@Logic: This is R script for running FBA simulation for small FBA model for Demo		#
#This plotter generating dotplots and pairwise plot		#
##################################################################################################

print("****************WELCOME TO BivariateGRAPH_PLOTTER***************************")

###Calling R LIM simulator
limfile_name = "TCA-analysis.lim"

#Output file TCA-analysis.fxr file
sp <- c(strsplit(as.character(limfile_name),"\\."))
filename<- grep(".+[^lim]", sp[[1]], value=TRUE)
fxrfilename <-as.character(filename)
ff <- c(fxrfilename,".fxr")

out1 <-paste(ff, collapse = '')			#This is string concatenation in R

###Calling R LIM simulator and capturing output in .fxr file
print("========================================================================")

print("*********Estimating optimal solution using linear programming*************")
require(LIM)
LIMEcoli <- Setup(limfile_name)
pars <- Ldei(LIMEcoli)
LP <- Linp(LIMEcoli)
data.frame(t(LP$X))
out<-capture.output(data.frame(t(LP$X)))
#setwd(newdir)				#This is for changing directory and depositing .fxr file in that directory
cat(out,file=out1,sep="\n",append=TRUE)

print(paste("Output .fxr file Generated ", as.character(out1)))

print("========================================================================")
print("**********Estimating ranges of all unknown reaction rates****************")
xr <- Xranges(LIMEcoli)
data.frame(simplest = pars$X, optimal = LP$X, xr)

print("**********Generating Dotchart for showing optimal solution within range*************")

sub("^\\s+", "", fxrfilename)
pdf_name <- paste(fxrfilename,"_","Dotchart_with_optimal_solution_within_Range",".pdf")
pdf(pdf_name)

print(paste("Generated dotplot : ",pdf_name))

print("========================================================================")

#Extracting column names
Rxn_names <- colnames(LP$X)

#For printing range and giving optimal solution
par(mfrow = c(1, 2))
nr <- LIMEcoli$NUnknowns

ii <- 1:(nr/2)
dotchart(LP$X[ii], xlim = range(xr), pch = 16, cex = 0.8,labels = Rxn_names[ii])

segments(xr[ii,1], 1:nr, xr[ii,2], 1:nr)

iii <- (nr/2+1):nr
dotchart(LP$X[iii], xlim = range(xr), pch = 16, cex = 0.8,labels=Rxn_names[iii])
segments(xr[iii,1], 1:nr, xr[iii,2], 1:nr)
mtext(side = 3, cex = 1.5, outer = TRUE, line = -1.5,"E coli Core Metabolism, optimal solution and ranges")
# It has genereated plot and saved as .pdf file
dev.off()

print("===========================================================================")

# user defined input number of iterations tobe run
#num_itr <- readline(prompt="You want to run simulation for how many number of iterations: ") 
#num_itr_ii <-c("Number of random samples considered for analysis is: ",num_itr)
#print(num_itr_ii)
num_itr = 500
xs <- Xsample(LIMEcoli, iter = as.numeric(num_itr), type = "mirror", test = TRUE)
#pairs(xs[, c(1,4)], pch =".", cex = 2, gap = 0, upper.panel = NULL)

print("===========================================================================")

## put histograms on the diagonal
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5))
    h <- hist(x, plot = FALSE)    
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
    # rect(breaks[-nB], 0, breaks[-1], y, col="blue", ...)
}

color_vector <- c("light blue", "yellow","green","red","black","brown","gold","dark cyan")

print (paste("Number of reaction fluxes in model: ",nr))

print(paste("List of reaction variables:"))
print(colnames(xs[,c(1:as.numeric(nr))]))


#Plotting first user defined two variables
pdf_name2 <- paste(fxrfilename,"_","Pairwise_plot","_for_two_variables",".pdf")
pdf(pdf_name2)	
	
#rxn1 <- readline(prompt="Please enter two reaction variable names one by one: Rx1:") 
#rxn2 <- readline(prompt="Please enter two reaction variable names one by one: Rx2:") 
rxn1 <- "MEMBRANETRASPORT12"
rxn2 <- "GLYEMP12"
print(paste("Number of reaction considered for analysis is: ",rxn1,rxn2))
	#sp <- array(c(strsplit(as.character(two_rxns),",")))
		
pairs(xs[,c(as.character(rxn1),as.character(rxn2))], pch =21, bg = color_vector[1:2], cex = 1, gap = 1, upper.panel = NULL, panel = panel.smooth , diag.panel = panel.hist )
	
par(xpd=TRUE)
legend("topright", colnames(xs[,c(as.character(rxn1),as.character(rxn2))]) ,fill = color_vector[1:2] ,inset=.001, title="Data sets")
		
dev.off()
	
print("Generated Pairwise plot for two reaction variables. ")
print ("*****************************************************")
	






