# Debug plots for C flow data preprocessing

# Clean up & initialize
rm (list = ls ())
setwd ('~/Work/Code/C/2012/mpiStochEst')

# Set dimensions
M = 80000
K = 5
P = 1

# Read data from the C program
cdata <- scan ('test.txt')
cdata <- matrix (cdata, K*P, M)
cdata <- t (cdata)

# Change working dir
setwd ('/Volumes/lillaccg$/Bio Data/Flow/Lac-GFP expt 3/Gated TXT')

# Read and process data of overnight culture
ovn <- read.delim ('overnight.txt')
ovn <- ovn[,7]
ovn <- ovn[ovn>1]
ovn <- log10(ovn)

# Read and process data of background culture
bg <- read.delim ('CTRL overnight.txt')
bg <- bg[,7]
bg <- bg[bg>1]
bg <- log10(bg)

# Resample background population
sambg = sample (bg, length (ovn), replace = T)

# Perform the deconvolution
dec = ovn - sambg

# Set up plot
quartz ()
par (mar = c (4, 4, 1, 0) + 0.1, mfrow = c (1,5))

for (i in seq (1, K*P))
{
	plot (density (cdata[,i]), main = NA, xlab = NA, xlim = c(-1,4))
	if (i==1)
	{
		lines (density (dec), col = 'red')
	}
}