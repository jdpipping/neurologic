################
### PACKAGES ###
################

# install.packages("haven")
library(haven)
# install.packages("tidyverse)
library(tidyverse)

#################
### LOAD DATA ###
#################

# DEMOGRAPHICS #

# load
demo_2011_2012 = read_xpt("data/nhanes/demographics/DEMO_G.xpt.txt")
demo_2013_2014 = read_xpt("data/nhanes/demographics/DEMO_H.xpt.txt")
# preview
head(demo_2011_2012)
head(demo_2013_2014)

# HIQ #

# load
hiq_2011_2012 = read_xpt("data/nhanes/hiq-huq/HIQ_G.xpt.txt")
hiq_2013_2014 = read_xpt("data/nhanes/hiq-huq/HIQ_H.xpt.txt")
# preview
head(hiq_2011_2012)
head(hiq_2013_2014)

# HUQ #

# load
huq_2011_2012 = read_xpt("data/nhanes/hiq-huq/HUQ_G.xpt.txt")
huq_2013_2014 = read_xpt("data/nhanes/hiq-huq/HUQ_H.xpt.txt")
# preview
head(huq_2011_2012)
head(huq_2013_2014)

# STROKE #

# load
stroke_2011_2012 = read_xpt("data/nhanes/stroke/MCQ_G.xpt.txt")
stroke_2013_2014 = read_xpt("data/nhanes/stroke/MCQ_H.xpt.txt")
# preview
head(stroke_2011_2012)
head(stroke_2013_2014)

# TBI #

# load
tbi_2011_2012 = read_xpt("data/nhanes/tbi/CSQ_G.xpt.txt")
tbi_2013_2014 = read_xpt("data/nhanes/tbi/CSQ_H.xpt.txt")
# preview
head(tbi_2011_2012)
head(tbi_2013_2014)