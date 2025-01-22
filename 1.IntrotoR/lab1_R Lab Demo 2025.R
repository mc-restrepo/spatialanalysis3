################################################################
### Introduction to the R Statistical Computing Environment  ###
################################################################

### Welcome to R! This source file is intended to provide a brief introduction 
### to this software and how best to navigate the computing environment.

### Remember that you will not be expected to write your own code for this class,
### however you will need to run, comprehend, and at times manipulate the code 
### that has been provided to you.


### For those of you using Rstudio take a moment to familiarize yourself with the layout. 
### There will be 4 panes, Source, Console, Environment et al., & Files et al.
### These panes can be changed using Preferences (Macs) or Tools -> Global options (PCs).
### Global options -> Appearance will also allow you to change the font size, window color,
### and text color.

### The first thing you need to know is how to run a command line in R
### This is done in PCs by pressing Ctrl+Enter and in Macs by pressing Command+Enter
### If you are using R studio you can also run the entire source file by pressing 
### the "Run" button located at the top-right of the source page, but I wouldn't recommended 
### that for this course. You can also select parts of lines or multiple lines to run if you'd like.

#---- Now let's perform some basic mathematical calculations in R #----

1 + 1   # Notice how the answer appears in the console window,
#  Your source window remains the same.

2 ^ 3

sqrt(81)  # Using one the many "functions" that can manipulate numbers

#Here's a fun one your high school math teacher probably gave you...
sin(pi / 6) * (sqrt(4)) + factorial(3) ^ 2 - (46656) ^ (1 / 3) 
# ... much easier then typing this into a calculator!

### Hopefully you've noticed that if you try running anything with a "#" sign before it 
### R will not see this as code. This is so you can annotate code and is extremely 
### useful both for this class and for your own coding. R can speed up this process for you
### by "commenting out" sections of selected code by pressing CTRL+SHIFT+C

### While R can act as a calculator, it is really meant to be an object-oriented
### computing environment. This means that you can assign large pieces of 
### information, such as vectors, to simple symbols... 

#for example:

x <- c(1,2,3,4,5,6,7,8,9,10)
y<-c(3,2,5,4,6,5,7,6,12,8)

### Note that R accepts both lines above - the spaces in between arguments are ignored.
### R will auto-format and clean up your code to insert spaces and tabs if you'd like, by
### selecting a chunk of code and pressing CTRL+SHIFT+A.


### In Rstudio you will see your objects stored under the "environment" window

x  # Now running a line that calls the first object "x" will display the data in the Console

y  # The same is true for "y"
 
X  # Objects are case-sensitive, so be careful!
 
length(x)  # Using the length function on the first object, there are many functions you can use on objects

# The following functions give various descriptive statistics.  All should be self-explanatory
mean(y)
 
median(y)

var(y)

stdy<-sqrt(var(y))  # Calculating the standard deviation as the square root
                    # of the variance and assigning it to a variable called "stdy"
stdy

summary(y)  # A particularly useful summary function.

z <- c(x, y)  # Concatenating two vectors into one long vector
z


# Now let's make a plot of x versus y.  R will automatically open a graphics window.
# and R studio will use the "Plots" tab.
plot(x, y)

#  Let's plot again, but  re-label the axes, you can change the names to whatever you like
plot(x,
     y,
     xlab = "X Variable",
     ylab = "Y Variable",
     type = "s")

# This is an example of producing a histogram with your own title
hist(y, main = "This is a Histogram of Y")

#  In Rstudio you can use the arrows to look at all the plots you've made since you opened the software.
#  You can also zoom in and adjust the window size of your plot.


xymatrix <- cbind(x, y)  # The cbind function combines columns of data into a matrix
                         # Note that Rstudio can allow you to visually explore matrix/dataframe objects:
                         # click on the table icon associated with the xymatrix in the Environment tab to explore

xymatrix  # Type it and hit return to see what it looks like

dim(xymatrix)  # Prints the dimensions of this matrix

### Square brackets "[]" allow us to index within objects (vectors, dataframes, lists, etc.). For single-
### dimensional objects (like vectors), we can just specify the position we're interested in.
y[5]  # Prints the 5th item in y

### For 2-dimensional objects (like dataframes), we need to specify the row and column if we want to pull a
### specific value. Since x and y have been bound as columns in xymatrix, we can type
xymatrix[5, 2]  # Prints the item in the 5th row in the 2nd column. Indexing a 2-dimensional object always lists
                # the row first, then the column, separated by a comma.

rm(xymatrix)  # This removes the matrix from the database
xymatrix  # The object should now be gone, you will see it disappear from the Environment pane as well

# If you want to remove all objects from your Environment, you can also select the broom icon in the Environment pane 
# I would only recommend if you want to start over with your analysis or are moving to a different project or source file


### R can be used to simulate probability distributions, 
### This is why it is a "statistical" computing environment.

e <- rnorm(1000, 0, 500)  # This generates Normal data with mean of zero and stdev 500
hist(e)                   # looks pretty normal

(x <- 1:1000)  # x is now re-assigned the integers 1 to 1000.
               # The parentheses prints the output while assigning it, 
               # it's as if you had afterwards typed "x" and ran the line

y <- 200 + 4 * x + e  # Assignment to the variable y, actually simulating a simple
                      # linear regression model with outcome y, predictor x,
                      # intercept=200, and slope=4. See it?

plot(x, y)  # Plot of x versus y
cor(x, y)  # The correlation coefficient between x and y

modelfit <- lm(y ~ x)  # This runs a linear regression with outcome y and predictor x
summary(modelfit)  # This generates a summary of the fit for your regression

class(modelfit)  # The linear regression function, lm, actually creates a 
                 # complicated object  (class "lm"). 
                 # In R studio it will be described as a class "List" with 12 items

names(modelfit)  # The names function provides the titles for 
                 # each list item, and each of these titles has it's own 
                 # information ... as you can see, objects in R can be very complex!

abline(modelfit$coefficients, lwd = 2, col = "red")  # The function "abline" adds a line to your last generated plot.  
                                                     # "lwd" controls thickness and "col" controls the color.  
                                                     # Here we are just adding the fitted regression line.


# We use the "$" to access different items from the list. (remember from "names"?)
modelfit$coefficients  # The fitted coefficients, we saw these in "summary" as well

plot(x, modelfit$residuals)  # Plot the residuals
abline(h = 0, lwd = 2, col = "blue")  # Add a horizontal line at zero

help(lm)    # You'll get to know this function well
            # It should open up another window with the help for lm
?lm         # I prefer this short hand version!


#---- Now let's practice reading in a file and creating a very simple map #----

### First you need to set the working directory for your current R session

## If you are using a Mac it will look something like this...
#setwd("/Users/bjdavis/OneDrive - Johns Hopkins University/SAGIS_1819/SAGIS3")

## If you are using a PC it will look more like this
#setwd("C:/Users/aunek/Dropbox/Hopkins/TA/140.698 - Spatial Analysis III/2025/Lab 1")

### If you are using Rstudio you can also use the dropdown "Session" and set working directory.
### If you want to set it to where your source code is you can select "To Source Files Location"
### Otherwise you can choose directory. You can also set your working directory through the Files
### pane under the "More" tab.

### Make sure you set it to whatever folder you saved the ozone data for this lab
### (ideally in your class folder).

#I'm using the here package instead

library(here)

# To check the current working directory
getwd()


# If you're in the correct working directory you can just read in the file using this 
ozone <- read.csv(here("1.IntrotoR", "Attribute_Ozone2007.csv"))

# The following command prompts you to navigate R to the file if you're feeling lazy
# and this is how we will read the ozone file in for this lab demo as it is stored
# on the desktop
ozone <- read.csv(file = file.choose())

### Also if you are using Rstudio, and you can click "Import Dataset" 
### in the Environment tab to bring in the csv file directly.
### For the ozone dataset, select "From Text (base)..." change the name to "ozone", 
### and make sure for header you select "yes".


# Let's explore this data file we just read in
names(ozone)
dim(ozone)
class(ozone)  # A dataframe is a special R object, the one we use most for basic regression analyses

# We can also use indexing to find out the classes of our columns
class(ozone[, 1])  # The type of class for each column
class(ozone[, 2])
class(ozone[, 4])

# Useful functions to get a quick glimpse at what our data look like
head(ozone)  # All columns of the first rows (usually 5)
tail(ozone)  # All columns of the last rows 


# A particularly useful function is "str()" to look at the structure of the object. 
# Let's try it on this data.frame

str(ozone)  # Similar to class(ozone[, 2]), but a lot more efficient, right?

# It basically summarizes all of the functions we just ran!


# In Rstudio you can also click the dropdown button on the object to get a similar structure summary.
# Again, you can also further explore the dataframe by clicking the table icon in the Environment tab


#---- The R software you downloaded is actually just a base set of packages. There are #----
### close to 10,000 packages developed by users that have their own functions and 
### unique object classes. Think of R as statistical computing 2.0.
### For this course we will use a variety of packages, 
### let's practice by installing the "maps" package

# You can use this code below, or in R studio click Tools, then "Install Packages", or from the 
# Files pane, click on the Packages tab, and then the Install button.
#install.packages("maps")

# This may take a moment...

### You only need to install packages once (don't forget to update them every 
### now and then!). Every time you start up R you need to load the package library:

library(maps)  # This loads the "maps" packages for this R session

# In Rstudio you can alternatively go to the Packages tab to select/deselect installed packages


# Let's try out our new package

map('italy')
map('usa')
map('state')
map('county')

# We can now add points using the ozone monitoring station dataset we brought in
points(ozone$Long,  # The x position of our points
       ozone$Lat,  # The y position of our points
       pch = 16,  # A reference to the point shape (see http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r)
       col = "blue",  # Selecting the color; R has some presets, but you can also use HTML color codes (e.g. "#03fc49")
       cex = 0.6)  # Scaling the point size 
title("Ozone Monitoring Stations")


## Hopefully by now you're getting the hang of using R!

## To save a graphic in Rstudio click "Export" under the plots tab.
## You can save as an image or pdf and choose the dimensions and file location.

#-------------------------------------------------------------------------------------------------------------

# For this first lab save the map of the US counties with the ozone monitoring locations added (in pdf format)
# This was the last plot we made using the commands above

# Hand in this map by 10:00pm EDT Friday January 24 to the CoursePlus dropbox 

#-------------------------------------------------------------------------------------------------------------



# The commands below are optional for you to work through and you are encouraged to do so
# Some of the analysis below is related to Question 1 in Problem Set 1
# we will start covering this material in lecture on Thursday January 26.



# -------------------------------------------------------------------------
# Now lets do some spatial analysis using the Swiss Rainfall dataset in R
# -------------------------------------------------------------------------

library(geoR)	# This loads the package geoR into your current R session
  # First make sure you have installed that package onto your machine
  install.packages("geoR")
  # If you have used R before make sure all your packages are updated


data(SIC)  # this loads the Swiss rainfall dataset into your workspace
?SIC # see what different data objects you have
names(sic.all)  # component names in the list sic.all

# The Swiss rainfall data is called the SIC data because it was used in a
# Spatial Interpolation Comparison (SIC) project back in the late 1990's.

# What is the dimension of the coordinates component?
dim(sic.all$coords)

plot(sic.all$coords, pch = 16)  # R by default plots the first column versus
                                # the second column when given a matrix.
                                # In case you didn't know this the command
                                # below does the same thing, just with
                                # different axis labels. Note I also changed
                                # the plotting symbol using pch to be a
                                # filled in circle.
plot(sic.all$coords[, 1], sic.all$coords[, 2], pch = 16)	# "[, 1]" is shorthand for saying use
									                                        # all the rows and just the column 1
									                                        # Similarly for "[,2]"

# Can you guess what is going on in these 4 commands below? (Remember the 
# help command to learn more about functions you aren't familiar with)
summary(sic.all$data)

points.geodata(sic.all,
               pt.divide = "quartiles",
               borders = sic.borders,
               x.leg = 0,
               y.leg = 290)
  
hist(sic.all$data, nclass = 20, xlab = "Rainfall")

plot(sic.all)  # Look at the nice set of plots R produces for this spatial data

# The plots in the top left and bottom right are self explanatory, although it would
# be quite helpful to have a legend automatically provided with the top left plot. The
# other two plots, top right and bottom left, are plotting the data versus each coordinate.
# In the top right plot the data is actually on the x-axis and the Y Coordinate of the
# aquifer points on the y-axis. This is so you can connect to the plot to its immediate
# left since the y-axes are the same.  Similarly for the plot in the bottom left. The  
# x-axis for this one and the plot above it are the same.


### Now to create a semivariogram for the dataset

rain.vario <- variog(sic.all) # Enter ?variog in the console to find out more about the object you created
names(rain.vario)
rain.vario$max.dist
rain.vario <- variog(sic.all, max.dist = 335.7076 / 2) # Restricting to half the actual maximum distance
# You can also use variables in functions
max.d <- summary(sic.all)$distances.summary["max"] # Storing the maximum distance as a variable
rain.vario <- variog(sic.all, max.dist = max.d / 2) # Using our new variable to create the semivariogram
plot(rain.vario)
title("Semivariogram of Swiss Rainfall Data")

rainfall <- sic.all$data # Storing rainfall as a standalone vectore
altitude <- sic.all$covariate$altitude # Repeating for altitude

plot(altitude, rainfall, pch= 16) # Is altitude a good predictor of rainfall?

modelfit2 <- lm(rainfall ~ altitude)

summary(modelfit2)

abline(modelfit2$coeff, lwd = 2, col = "red")

# Let's check out the distribution of the linear regression residuals
plot(altitude, modelfit2$residuals)
abline(h = 0, lwd = 2, col = "blue")

# The package geoR likes to have its geostatistical data sets all conform to
# be of a specific type and structure, for example, coordinates, followed by data, etc.
# So there is a function called as.geodata which makes this type of object.  The SIC data
# was already internally saved as this type of object.

resid.geo <- as.geodata(cbind(sic.all$coords, modelfit2$residuals)) # Making a geodata object out of
                                                                    # the residuals by concatenating
                                                                    # the coordinates with the residuals
names(resid.geo)

plot.geodata(resid.geo)  

# The following command does the same thing, R is smart
plot(resid.geo)

# These next two commands do the same thing.  There is a special function
# called points.geodata that just works on geodata objects.  But R is smart
# and is an object oriented programming language so when we just use the points
# command it sees that we gave as input a geodata object and so it knows to
# internally call the points.geodata function.  The same thing happened above
# when I used plot(sic.all) to really represent plot.geodata(sic.all)

points.geodata(resid.geo, pt.divide = "quartiles", borders = sic.borders)
points(resid.geo)  # R is smart 
lines(sic.borders) # Country boundary is provided with SIC

resid.vario <- variog(resid.geo)

resid.vario$max.dist	# Want to find the max inter-point distance

resid.vario <- variog(resid.geo, max.dist = 335.7076 / 2)  	# Re-estimate the variogram using
								                                          	# half the max inter-point distance

plot(resid.vario,
     pch = 16,
     xlab = "Distance (Km)", # We know the units are km from ?sic
     ylab = "Semivariogram")
title("Residual Semivariogram")

par(mfrow = c(1, 2))		# Set up the graphics window as a 1 by 2 matrix (you can turn this off 
                        # with 'dev.off()' or 'par(mfrow = c(1, 1))')

plot(rain.vario,
     pch = 16,
     xlab = "Distance (Km)",
     ylab = "Semivariogram")
title("Rainfall Semivariogram", cex = 0.3)
plot(resid.vario,
     pch = 16,
     xlab = "Distance (Km)",
     ylab = "Semivariogram")
title("Residual Semivariogram", cex = 0.3)

# Note the two plots look the same.  However they are very close but not exactly the same.  
# Here is a way to see the actual estimates.

temp <- data.frame(Distance = rain.vario$u,
                   RainVario = rain.vario$v,
                   ResidVario = resid.vario$v)

temp  # Type the newly created object temp and hit return

# Do you think the variable altitude accounted for much of the spatial variation
# in rainfall?  Explain.

par(mfrow = c(1, 1))		# Resets the plotting region to be for one graphic

plot(temp$Distance, temp$RainVario, pch = 16, col = "black")
points(temp$Distance,
       temp$ResidVario,
       pch = 1,
       col = "blue")


# Here are a set of commands to first "eyeball" a set of parameter estimates and then tell
# R to draw that corresponding semivariogram function on the current plot.  Below I am
# suggesting a spherical model with range of 85, partial sill of 14,000 and nugget of 1000.

plot(rain.vario,
     pch = 16,
     xlab = "Distance (Km)",
     ylab = "Semivariogram")
title("Rainfall Semivariogram", cex = 0.3)
rang <- 85
psill <- 14000
nug <- 1000
lines.variomodel(cov.model = "spherical",
                 cov.pars = c(psill, rang),
                 nugget = nug,
                 max,
                 max.dist = 335.7076 / 2)

# For practice change the values in the above set of commands for the range, partial sill, and nugget
# (rang, psill, nug) to see the different lines that get drawn.  If you then run the commands starting
# with "rang<-....." command then previous lines will remain on the plot.  Start from "plot(rain.vario..."
# instead to make one clean version with your final decided upon model.


# -----------------------------------------------------------------------------
# Some more spatial stuff using the Wolfcamp Aquifer data which is also in geoR
# -----------------------------------------------------------------------------

data(wolfcamp)

?wolfcamp  # This will provide some very basic and very brief background for this data

summary(wolfcamp)

plot(variog(wolfcamp, max.dist = 436.2 / 2))  # Variogram estimation and plotting in
title("Semivariogram of Pressure")  # The same command, however the computed
							                      # variogram estimate is not saved in any object
							                      # Did you see where I got the max dist from?

plot(wolfcamp)

# Do you see evidence of any large scale spatial trend as a function of 
# one or both coordinates?  The plots that get generated from the above command
# (4 plots on one page) are quite useful for looking at trends in coordinates.
# The upper left plot is mapping the quartiles of the data (a legend would be nice).
# The plots in the upper right and bottom left are looking at trends in both coordinates,
# the upper right is the data (x-axis) plotted against the Y or north coordinate (y-axis),
# the axes are switch like this so you can read the upper left and upper right plots
# together, see how the y-axis in both of these plots match up.  That is you can look at 
# the map (upper left) and see the trend moving north to south.  The plot in the bottom left 
# is similar, now plotting the X or east coordinate (x-axis) against the data (y-axis), and see 
# how now the two x-axes line up so this plot can be read with the map above it.  The plot in 
# the bottom right is a histogram of the data.


