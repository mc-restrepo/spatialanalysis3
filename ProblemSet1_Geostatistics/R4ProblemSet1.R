# ------------------------------
# Problem Set 1: Geostatistics
#
# R Commands for Question 1
# ------------------------------


# Note: The commands below are the same as those that appear in the
# second half of the R Spatial Lab Demo command file

library('geoR')	# This loads the package geoR into your current R session
# First make sure you have installed that package onto your machine, use annotated code below or Rstudio dropdown
# install.packages("geoR")
# If you have used R before make sure all your packages are updated (the latest version of RStudio will
# automatically check this for you - if this didn't happen, consider updating RStudio as well)

data(SIC)  # this loads the Swiss rainfall dataset from the geoR package into your workspace
? SIC # A description of the included data objects
names(sic.all)  # component names in the sic.all geodataset

# The Swiss rainfall data is called the SIC data because it was used in a
# Spatial Interpolation Comparison (SIC) project back in the late 1990's.

# what are the dimensions of the coordinates component?
dim(sic.all$coords)

plot(sic.all$coords, pch = 16)
# R by default plots the first column versus the second column when given a matrix. In case you
# didn't know, this the command below does the same thing, just with different axis labels.  Note
# I also changed the plotting symbol using pch to be a filled in circle.
plot(sic.all$coords[, 1], sic.all$coords[, 2], pch = 16)
# "[,1]" is shorthand for saying use all the rows and just the column 1. Similarly for "[,2]".
plot(sic.all$coords[, "V2"], sic.all$coords[, "V3"], pch = 16)
# "[, "V2"] can also be used to call the "V2" column by name instead of its index number. Similarly
# for "[, "V3"]. Make sure not to forget the comma (the first object in square brackets refers to cells,
# the second object refers to columns).

# Can you guess what is going on in these 4 commands below?

summary(sic.all$data)

points.geodata(sic.all,
               pt.divide = "quartiles",
               borders = sic.borders,
               x.leg = 0,
               y.leg = 290)

hist(sic.all$data, nclass = 20, xlab = "Rainfall")

plot(sic.all)  # look at the nice set of plots R produces for this spatial data

# The plots in the top left and bottom right are self explanatory, although it would
# be quite helpful to have a legend automatically provided with the top left plot. The
# other two plots, top right and bottom left, are plotting the data versus each coordinate.
# In the top right plot the data is actually on the x-axis and the Y Coordinate of the
# aquifer points on the y-axis.  This is so you can connect to the plot to its immediate
# left since the y-axes are the same.  Similarly for the plot in the bottom left.  The
# x-axis for this one and the plot above it are the same.


### Now to create a semivariogram for the dataset

rain.vario <- variog(sic.all) # enter ?variog in the console to find out more about the object you created
names(rain.vario)
rain.vario$max.dist
rain.vario <- variog(sic.all, max.dist = 335.7076 / 2) # Restricting to half the actual maximum distance
plot(rain.vario)
title("Semivariogram of Swiss Rainfall Data")

rainfall <- sic.all$data
altitude <- sic.all$covariate$altitude

plot(altitude, rainfall, pch = 16)   # Is altitude a good predictor of rainfall?

modelfit2 <- lm(rainfall ~ altitude) # Simple linear regression model with altitude regressed on rainfall

summary(modelfit2)

abline(modelfit2$coeff, lwd = 2, col = "red")

plot(altitude, modelfit2$residuals) # Plot the residuals to check your linear regression assumptions and model fit
abline(h = 0, lwd = 2, col = "blue")

# The package geoR requires to have its geostatistical dataset objects all conform to
# a specific type and structure (class geodata), for example, coordinates, followed by data, etc.
# So there is a function called as.geodata which makes this type of object.  The SIC data
# was already internally saved as a geodata object.

resid.geo <- as.geodata(cbind(sic.all$coords, modelfit2$residuals))
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

resid.vario <- variog(resid.geo) # Creating a semivariogram for the

resid.vario$max.dist	# Want to find the max inter-point distance

resid.vario <- variog(resid.geo, max.dist = 335.7076 / 2)
# re-estimate the variogram using half the max inter-point distance

plot(resid.vario,
     pch = 16,
     xlab = "Distance (Km)",
     ylab = "Semivariogram")
title("Residual Semivariogram")

par(mfrow = c(1, 2))		# set up the graphics window as a 1 by 2 matrix

plot(rain.vario,
     pch = 16,
     xlab = "Distance (Km)",
     ylab = "Semivariogram")
title("Rainfall Semivariogram", cex = .3)
plot(resid.vario,
     pch = 16,
     xlab = "Distance (Km)",
     ylab = "Semivariogram")
title("Residual Semivariogram", cex = .3)

# Note the two plots look identical.  However, they are very close but not exactly the same.
# Here is a way to see the actual estimates.

temp <-data.frame(Distance = rain.vario$u,
                  RainVario = rain.vario$v,
                  ResidVario = resid.vario$v)

temp # The newly created object temp

# Do you think the variable altitude accounted for much of the spatial variation
# in rainfall? Explain.

par(mfrow = c(1, 1)) # Resets the plotting region to be for one graphic

plot(temp$Distance, 
     temp$RainVario, 
     pch = 16, 
     col = "black")
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
title("Rainfall Semivariogram", cex = .3)
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

# To save a graphic in Rstudio click "Export" under the plots tab.
# You can save as an image or pdf and choose the dimensions and file location.
# Keep in mind that changing the size and dimensions of the plot pane in RStudio will change the
# size of your saved images for future plots you save.

# To save a graphic in base R, first have the graphic window upfront and then click on File
# on the top menu bar for R and go to "Save as" and there is an option for jpg.


# -------------------------------------------------------------------------------------------------
# Some more spatial analysis using the Wolfcamp Aquifer dataset which is also provided in geoR
# -------------------------------------------------------------------------------------------------

data(wolfcamp)

? wolfcamp  # This will provide some very basic and very brief background for this data

summary(wolfcamp)

plot(variog(wolfcamp, max.dist = 436.2 / 2)) # the max distance was given in the summary values
title("Semivariogram of Pressure")
# variogram estimation and plotting in the same command line, however the computed
# variogram estimate is not saved in any object. Did you see where I got the max dist from?

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


# Trend Analysis and Residual Spatial Variation

pressure <- wolfcamp$data
xcoord <- wolfcamp$coords[, 1]
ycoord <- wolfcamp$coords[, 2]

modelfit3 <- lm(pressure ~ xcoord + ycoord)

summary(modelfit3)

resid.geo <- as.geodata(cbind(wolfcamp$coords, modelfit3$residuals))
names(resid.geo)

plot(resid.geo)
plot(variog(resid.geo, max.dist = 436.2 / 2))

# Do you understand why the max distance of 436.2 did not change.  You can verify that it did not
# by running the command summary(resid.geo)

par(mfrow = c(1, 2))
plot(variog(wolfcamp, max.dist = 436.2 / 2))
title("Semivariogram of Pressure")
plot(variog(resid.geo, max.dist = 436.2 / 2))
title("Residual Semivariogram of Pressure")


# ------------------------------
# Problem Set 1: Geostatistics
#
# Some R Commands for Question 2
# ------------------------------

# load it again if starting from here; if continuing from above, this is redundant
library(geoR)

# Here is one way to read data into R.  The file Q2 Data.csv I provided is a .csv file (comma
# separated) that is easily generated from Excel.  Change the path to match where you have this file
# saved. See 'help(read.csv)' for more details.

# Set the working directory, where the data is to be read in. You choose what goes in between the brackets.
# Put in quotes. If you copy and paste from a file explorer address bar, make sure to check what kind of
# slashes you have (R only reads forward slashes in file locations).

setwd("D:/Documents/TA/Spatial Statistics/R_Projects/Problem_Set_1/")
setwd("C:/Frank/Spatial Analysis III/2018/Problem Sets/Geostatistics/Data and Files to Post")

# To check the working directory in your current R session
getwd()

#note, make sure your file has the same name, when downloading off of courseplus the name may change!
data2 <- read.csv("Q2 Data.csv", header = TRUE)

# You can also read in files with the command below, or you can use the Import Dataset option from
# Rstudio (See below).
data2 <- read.csv(file = file.choose())

# The read.csv command put data into what is called a dataframe in R, which is just a fancy looking
# matrix with column heading and possible row labels.

# In RStudio you can also click "Import Dataset" from the Environment tab. Select "From Text File..."
# then in the dialog box make sure the Name is "data2" and the radio button for "Heading" is "Yes".

# Use the names command to see whats in object data2
names(data2)

# Extract out your variables for easy use.  Note this is optional as we could just as well always
# refer to y1 as data2$y1.

y1 <- data2$y1
y2 <- data2$y2
coordx <- data2$coordx
coordy <- data2$coordy
x1 <- data2$x1
x2 <- data2$x2

# Create geodata objects
y1geo <- as.geodata(cbind(coordx, coordy, y1))
y2geo <- as.geodata(cbind(coordx, coordy, y2))

# Explore these data for large-scale spatial variation
plot(y1geo)
plot(y2geo)

summary(y1geo)
summary(y2geo)
# Output will proivde the maximum inter-point distance for each

# Note the shorthand using the variog function within the plot command
plot(variog(y1geo, max.dist = 1.3 / 2))
title("Semivariogram for Y1")
plot(variog(y2geo, max.dist = 1.3 / 2))
title("Semivariogram for Y2")

# More shorthand using the lm function within the cbind command and
# extracting the residuals and then all that is within the as.geodata
r1geo <- as.geodata(cbind(coordx, coordy, lm(y1 ~ x1)$residuals))
r2geo <- as.geodata(cbind(coordx, coordy, lm(y2 ~ x2)$residuals))

plot(variog(r1geo, max.dist = 1.3 / 2))
title("Residual Semivariogram for Y1")
plot(variog(r2geo, max.dist = 1.3 / 2))
title("Residual Semivariogram for Y2")


# ------------------------------
# Problem Set 1: Geostatistics
#
# Some R Commands for Question 3
# ------------------------------

# If you haven't before installed packages, you may need to install them before
# loading them for this assignment.  Once installed, you won't have to run
# these lines again, and just use the "library" command to reference them
install.packages("gstat")
install.packages("stars")
install.packages("sf")
install.packages("splines")

# Once packages are installed, load them into your current R session

library(geoR)
library(gstat)
library(stars)
library(sf)
library(splines)

# Resetting the plot window
par(mfrow = c(1, 1))

# this will turn off the device, which will also reset the plot window to default
dev.off()

# Again, set the working directory - where your data files are saved - if needed
# if you're running code from above, you don't have to repeat this step
setwd("C:/Frank/Spatial Analysis III/2018/Problem Sets/Geostatistics/Data")
# Make sure the zipped files are unzipped and in your working directory folder

# Import the State data

# In the commands below you are reading in GIS Shapefiles in R and working with data
# contained within these shapefiles.  R has tremendous capabilities now to work with
# these types of file (very convenient).  Remember though a "shapefile" doesn't necessarily
# correspond to a single file, rather a collection of files.  And although in R you are
# reading in only the ".shp" file all the other ones need to be contained in the same
# directory as the .shp file for this to work.

#' In our class assignments we will mostly use the "sf" package to handle and plot spatial objects
#' "sf" stands for Simple Features and is a more recent way of managing spatial data which is more
#' similar the way you would work with a data frame.
#' Previously almost all spatial object management and analysis used the "sp" package
#' which structures spatial data by separating the attribute information from the info
#' about the geometry, projection, etc. in what are known as slots. 
#' Because SF hasn't been fully incorporated into all other packages for spatial statistics,
#' we will sometimes go back and forth between them to adhere to the required data
#' format for certain operations

# Load the shapefile into RStudio as an SF object
States <- st_read("States_reproj_lower48.shp")

class(States)

# the 'geometry' is referenced like a column but actually contains more complex spatial information
head(States$geometry)
# note that this is areal spatial data ('multipolygons')

# Try to plot the States
plot(States)
# with an SF object, R's default plot will plot all the attributes attached to the geometry
# (you should see a warning message that only X out of the total N attributes can be plotted at once)
# To just plot the geometry outline without creating a choropleth for any of the information, specify:
plot(st_geometry(States))


# Read in the Ozone points
Ozone <- st_read("Ozone_Monitors_2007_reproj.shp")

head(Ozone$geometry)
# the ozone monitors are points (geostatistical data)

# we can look at the details of shapefiles' projection information
st_crs(States)
st_crs(Ozone)
# these have the same projection so they will correctly plot onto of one another
# in future work, if you're having trouble getting multiple layers
# to show up in the same plot window, this is a good first place to check

# Map of ozone sampling locations
# Because we are now working with 2 SF objects, we will use the base call 'plot'
# (instead of 'points' command) and specify add=T so the points are on top of the states
plot(st_geometry(States))
plot(st_geometry(Ozone),
       col = "blue",
       cex = 0.5,
       pch = 16, add=T) # Note the size of the points is changed using 'cex = 0.5'
title("EPA AQS 2007 Air Monitoring Locations")

# The title may appear far from the map.  Below is a way how to adjust that.
# See the help form (?title) to learn how to adjust the title to your own preference.

plot(st_geometry(States))
plot(st_geometry(Ozone),
       col = "blue",
       cex = 0.5,
       pch = 16, add=T)
title("EPA AQS 2007 Air Monitoring Locations", line = -3)

# Now we'll use some of the functions provided geoR.
# First we need to create a 'geodata' object.

# Object "temp" is a dataframe object that contains the projected x and y coordinates and
# the outcome which is the daily 8 hour max (D8HOURMAX). Also when creating the geodata
# object I identified the x coordinate (easting) as a covariate "cov.col=1". This will be used later.
names(Ozone) # to see which columns you need to reference

temp <- data.frame(easting = Ozone$x_coord_m,
                   northing = Ozone$y_coord_m,
                   ozone = Ozone$D8HOURMAX)
ozone.geo <- as.geodata(temp, covar.col = 1)

# At this point I want to make a few comments.  The analysis we are presenting here in R is using a
# more sophisticated object class, namely the shapefiles.  As I mentioned in class we could have
# generated the needed attribute table in ArcGIS and export it out as a dbf and then save it as a .csv
# file and read that into R.  But I wanted to give an example of how to work with shapefiles.
# R is able to do much of what you learned in the GIS course.

# Here is an indicator map with ozone split into quartiles.  I used the cex.min and cex.max
# attributes of the function to keep all the dots the same size.  I also included a legend.

# notice the 'points' call can be used again here bc its in 'geodata' form now, not 'sf'
plot(st_geometry(States))
points.geodata(ozone.geo,
               pt.divide = "quartiles",
               cex.min = .8,
               cex.max = .8,
               x.leg = -2529176,
               y.leg = -1200000,
               add.to.plot = TRUE)
title("Average Annual Daily 8hr Max Ozone")

# In the past some students have had some trouble with the legend part of this plot, for example it
# not fitting well and getting cutoff in RStudio.  The x.leg and y.leg attributes above tell R where
# to place the legend and you can change this to hopefully bring up or move over your legend so it
# all fits.  Below is an alternative approach.  I save the output from the points.geodata function
# into an object I just called xx which you can see what is in it by the names(xx) command.  This still
# produces the plot but I took out the x.leg and y.leg attributes so no legend is created.  Now I can use the
# legend command to make one myself but I need the cutoff values which I can calculate by running xx$quantiles
# (since this is how the color scale of the original plot is divided). In this example, I am giving the legend
# the attribute locator(1) which basically lets R know that you will tell it where the legend goes.  So
# move your mouse over to the plot and click once. It usually takes a few times to get the location just right
# so just keep re-plotting till you get it. After running "locator," you cannot move on to other 
# commands until you have clicked the plot (you'll notice the cursor over the plot is a crosshairs and 
# the next '>' in the console won't pop up until you're done)
# Finally you can adjust the font size in the legend by setting the
# cex attribute in the legend command. Here I am making it .9 so 90% of the default size, adjust yours
# accordingly to be smaller or bigger than 1.

plot(st_geometry(States))
xx <- points.geodata(ozone.geo,
                     pt.divide = "quartiles",
                     cex.min = .8,
                     cex.max = .8,
                     add.to.plot = TRUE)
names(xx)
xx$quantiles # just to see what the breaks are for the quantiles
legend(locator(1), 
       legend = c("[24.70,40.63)",
                  "[40.63,45.68)",
                  "[45.68,50.72)",
                  "[50.72,68.94)"),
       pch = 16,
       col = c("blue", "green", "yellow", "red"),
       cex = .9)

# Since ozone.geo is a geodata object, the plot command below knows to call up
# the routine that generates that spatial 2x2 set of plots.
plot(ozone.geo)

# Since this is such a large geographic area let's look at potential trend in coordinates
# on their own plot.  I was also messing around with the function smooth.spline to help
# highlight potential trends.

plot(ozone.geo$coords[, 1],
     ozone.geo$data,
     xlab = "Easting (meters)",
     ylab = "Ozone (ppb)",
     pch = 16,
     pty = "m")
lines(smooth.spline(ozone.geo$coords[, 1], ozone.geo$data, df = 15),
      col = "blue",
      lwd = 2)

plot(ozone.geo$data,
     ozone.geo$coords[, 2],
     xlab = "Ozone (ppb)",
     ylab = "Northing (meters)",
     pch = 16,
     pty = "m")
temp <- smooth.spline(ozone.geo$coords[, 2], ozone.geo$data, df = 4)
lines(temp$y, temp$x, col = "blue", lwd = 2)

# Now let's plot a semivariogram to look at small-scale spatial variation (over different distances)
maxdist <- summary(ozone.geo)$distances.summary[2]

ozone.vario1 <- variog(ozone.geo, max.dist = maxdist / 2)

ozone.vario2 <- variog(ozone.geo)

par(mfrow = c(1, 2))

plot(ozone.vario1,
     pch = 16,
     xlab = "Distance (meters)",
     ylab = "Semivariogram",
     pty = "m")
title("Ozone Semivariogram\n(restricted < 1/2 max distance)") 
#'\n' moves the following text to a new line

plot(ozone.vario2,
     pch = 16,
     xlab = "Distance (meters)",
     ylab = "Semivariogram",
     pty = "m")
title("Ozone Semivariogram\n(all pairwise distances)")

par(mfrow = c(1, 1))
# Eyeball a set of parameter estimates for a spherical semivariogram
psill = 35
range = 1000000
nug = 20
plot(ozone.vario1,
     pch = 16,
     xlab = "Distance (meters)",
     ylab = "Semivariogram",
     pty = "m")
title("Ozone Semivariogram\n(restricted < 1/2 max distance)")
lines.variomodel(cov.model = "spherical",
                 cov.pars = c(psill, range ),
                 nugget = nug,
                 max,
                 max.dist = maxdist/2)

# Code to fit a large-scale trend using natural splines for the easting coordinate
data2use <- data.frame(ozone = Ozone$D8HOURMAX, easting = Ozone$x_coord_m)
model0 <- lm(ozone ~ 1, data = data2use) # An intercept-only model (ozone~1)
summary(model0)
model1 <- lm(ozone ~ ns(easting, 4), data = data2use) # A model regressing a 4-knot spline of
  # the x-coordinate on ozone
summary(model1)
plot(data2use$easting,
     data2use$ozone,
     xlab = "Easting (meters)",
     ylab = "Ozone (ppb)",
     pch = 16)
data2pred <- data.frame(easting = seq(-2329176, 2244937, length = 1000)) # A new dataframe with 1000 x-coordinates
# along the range of observed x-coordinates (prediction point about every 4500 meters)
preds <-  predict.lm(model1, data2pred) # Predicting ozone at each of the 1000 locations from east to west
lines(data2pred$easting, preds, lwd = 3, col = "blue") # Plotting these predicted ozone values

# Now create another geodata object this time for the ozone residuals
ozone.resid.geo <- as.geodata(cbind(ozone.geo$coords,
                                    model1$residuals))

ozone.resid.vario <- variog(ozone.resid.geo, max.dist = maxdist / 2)
plot(ozone.resid.vario,
     pch = 16,
     xlab = "Distance (meters)",
     ylab = "Semivariogram")
title("Ozone Residual Semivariogram")

# Note that another way to plot the residual semivariogram is to include the trend in the Variog function,
# see annotated code below:

# ozone.resid.vario<-variog(ozone.geo,max.dist=maxdist/2, trend=~ns(easting,4))

# Again, eyeball a set of parameter estimates for a spherical semivariogram
# You can either re-assign values for "psill", "range" "maxdist" etc. or you can
# enter the values manually as you approximate the best-fitting line

# ------- #
# Kriging #
# ------- #

# Let's create a model using the parameters we've just estimated using the weighted least squares (WLS) approach.


# Use your eyeballed estimates as the initial parameter estimates.  I chose a spherical function.
# You are welcome to select another function  that appears to fit the shape, just make the appropriate
# change in the command below. The ini.cov.pars, provides initial covariance parameters (partial sill,range),
# As you see the initial "guess" for the nugget is supplied in a separate argument. I have included some
# initial values for all three parameters, replace these with your eyeballed estimates.
ozone.vario.wls <- variofit(ozone.vario1,
                            ini.cov.pars = c(30, 1000000),
                            cov.model = "spherical",
                            nugget = 20,
                            weights = "cressie")

ozone.resid.vario.wls <- variofit(ozone.resid.vario,
                                  ini.cov.pars = c(15, 750000),
                                  cov.model = "spherical",
                                  nugget = 25,
                                  weights = "cressie")

# The following just adds the WLS fitted semivariogram model to the plot
plot(ozone.vario1,
     pch = 16,
     xlab = "Distance (meters)",
     ylab = "Semivariogram",
     pty = "m")
title("Fitted WLS Semivariogram for Ozone Levels")
lines(ozone.vario.wls)

plot(ozone.resid.vario,
     pch = 16,
     xlab = "Distance (meters)",
     ylab = "Semivariogram",
     pty = "m")
title("Fitted WLS Semivariogram for Ozone Residuals")
lines(ozone.resid.vario.wls)

# You can put them both on the same plot as below.  Note the locator(1) deal in the legend
# command is going to prompt you where to place the legend.  Drop it by clicking the left
# mouse button.  It will probably take a few times to get it where you want.
plot(ozone.vario1,
     pch = 16,
     xlab = "Distance (meters)",
     ylab = "Semivariogram",
     pty = "m")
lines(ozone.vario.wls)
points(ozone.resid.vario$u,
       ozone.resid.vario$v,
       pch = 16,
       col = "blue")
lines(ozone.resid.vario.wls, col = "blue")
legend(locator(1),
  legend = c("Ozone", "Ozone Residuals"),
  lwd = c(1, 1),
  col = c("black", "blue"))
title("WLS Fitted Semivariograms for Ozone and Ozone Residuals")

# Create a grid of locations to predict at.  Below I'm just going 10000 meters on either side
# of the min and max easting and northing coordinate and then building a 100 by 100 grid for
# a total of 10,000 prediction locations.
summary(ozone.geo$coords) # the min/max of each coordinate are manually entered to generate a grid
grid <-  expand.grid(easting = seq(-2329176 - 100000, 2244937 + 100000, len = 100),
                     northing = seq(-1452329 - 100000, 1365129 + 100000, len = 100))

plot(st_geometry(States))
points.geodata(ozone.geo,
               pt.divide = "quartiles",
               cex.min = .8,
               cex.max = .8,
               add.to.plot = TRUE)
points(grid, pch = ".")

# To perform IDW I used the library gstat and the code below
# the input information is the Ozone sf object (which contains the observed data values at monitor locations)
# and we will predict to all locations of the grid we defined above after making it an SF object as well
temp.grid <- st_as_sf(grid, coords=c("easting", "northing"), crs=st_crs(Ozone))

# by typing 'gstat::' we are explicitly telling R to use the 'idw' command from the gstat
# package; R can figure this out automatically if there is only one 'idw' function loaded, 
# however, you can run into issues and need to be package-specific in your command call
# if you've loaded other packages with functions of the same name

IDW.pred <- gstat::idw(D8HOURMAX ~ 1, Ozone, temp.grid, idp = 2) # idp is the inverse distance weighting power

# Now, the result of this analysis gives an SF object (same as what we put in)
# If this were stand-alone result, we might go straight to plotting from here;
# However, to compare with the outputs from other spatial interpolation methods,
# read the notes below to see how we'll get all outcomes to plot together


# Trend surface model predictions are just predictions from a linear model
TS.pred <- predict.lm(model1, grid, se.fit = TRUE) # model1 is our 4-knot x-coordinate spline regression on ozone

# The following two commands produce ordinary and universal kriged predictions using the WLS semivariogram
# estimated function. In the ordinary kriging (first command) I supply the semivariogram ozone.vario.wls
# In the universal kriging (second command) I supply the semivariogram ozone.resid.vario.wls and also supply
# the form of the trend.  The function will internally fit the trend and use the residual semivariogram we supplied.
# This may take a second for R to complete. Check 'help(krige.conv)' to see what the function arguments mean.

OK.pred <- krige.conv(ozone.geo,
                      locations = grid,
                      krige = krige.control(obj.model = ozone.vario.wls))

UK.pred <- krige.conv(ozone.geo,
                      locations = grid,
                      krige = krige.control(obj.model = ozone.resid.vario.wls,
                                            trend.d = trend.spatial( ~ ns(easting, 4), ozone.geo),
                                            trend.l = trend.spatial( ~ ns(easting, 4), grid)))

# The code that follows maps the prediction results for the IDW, Trend Surface, and kriging
# approaches.  The library geoR has a function called image.kriging which produces maps
# of kriged predictions.  This function assumes its input is from a kriging analysis and hence
# a specific type of object is required for input to image.kriging.  IDW and Trend Surface though
# are not kriging and hence the output from these analyses are not stored in this specific type
# of kriging object.  So below I use a little trick to allow the image.kriging function to be used
# on results from IDW and Trend Surface.  I copy the ordinary kriging object (OK.pred) into a dummy
# temporary object and then just replace predictions in that object (i.e. the ordinary kriged
# predictions) with the IDW predictions.  Then I do the same for the Trend Surface predictions.  For
# the Trend Surface I also replace the prediction variances with the Trend Surface prediction variance.
# Then just use the image.kriging function on these temporary objects.

# The maps created below for all the predictions do not remove predicted locations that are outside
# the US borders.  As you will see we get a rectangle shape predictions.  I would have to read in
# another shapefile for just the outer border of the US and clip results to that boundary to fit results to frame 
# but I decided not to so you can also see the behavior of predictions and prediction standard errors as we move away from
# the sampled data.  Also, a point previously mentioned is that for publication/presentation purposes
# we would likely import these results into a GIS and make final maps there. While similar quality maps
# can be made in R, they can be more difficult to make and are beyond expectations for this course.

# You'll notice after displaying the results in the form of the kriging objects we can still add the 
# geometries from the original states and ozone objects (and because our work was done from their original
# projected data, the coordinates are in the same units and locations line up correctly)

plot(st_geometry(Ozone), pch=16, col="red") 
points(ozone.geo$coords, pch = 16, cex = .5)
# this confirmed the sf object pts are same locations as geodata object coords

IDW.pred.temp <- OK.pred
IDW.pred.temp$predict <- IDW.pred$var1.pred

plot(st_geometry(States)) # set the plot window in the correct location
image(x = IDW.pred.temp,
      locations = grid,
      axes = FALSE,
      xlab = "", ylab = "",
      breaks = seq(25, 65, 3.1),
      col = heat.colors(12),
      x.leg = c(-2329176, -1000000),
      y.leg = c(-1500000, -1200000))
title("IDW Predictions of Ozone", line = -1.5) # remember "line" move the location of the title; adjust if needed
plot(st_geometry(States), add = TRUE) # re-add the state outlines bc they were drawn over
plot(st_geometry(Ozone), add = TRUE, pch=16, cex=0.5)


## These image calls should come from the "graphics" package which
# loaded as a dependency from other packages in our library
# if it's having trouble plotting, double check you've correctly 
# created the kriging objects by calling "class(IDW.pred.temp)" for example
# and manually respecify the pacakge for the plot commange "graphics::image(...)"
# you can see which objects contain - names, classes, lists, etc. with "str(OK.pred)"

TS.pred.temp <- OK.pred
TS.pred.temp$predict <- TS.pred$fit
TS.pred.temp$krige.var <- TS.pred$se.fit ^ 2
image(x = TS.pred.temp,
      locations = grid,
      axes = FALSE,
      xlab = "", ylab = "",
      breaks = seq(25, 65, 3.1),
      col = heat.colors(12),
      x.leg = c(-2329176, -1000000),
      y.leg = c(-1500000, -1200000))
title("Trend Surface Model Predictions of Ozone", line = -1.5)
plot(st_geometry(States), add = TRUE)
plot(st_geometry(Ozone), add = TRUE, pch=16, cex=0.5)

### Note the different legend scale on this plot for S.E.s
image(TS.pred.temp,
      locations = grid,
      axes = FALSE,
      xlab = "", ylab = "",
      value = sqrt(TS.pred.temp$krige.var))
legend.krige(c(-2329176, -1000000),
             c(-1500000, -1200000),
             sqrt(TS.pred.temp$krige.var),
             offset.leg = .8)
title("Trend Surface Model Prediction Standard Errors", line = -1.5)
plot(st_geometry(States), add = TRUE)
plot(st_geometry(Ozone), add = TRUE, pch=16, cex=0.5)

# Click the zoom button to see the enlarged plot if you want
# for you assignment, you can screenshot each finished product
# or click the export, save as image to generate a separate file


# now the Kriging outputs, already in the kriging format from the command krige.conv()
image(x = OK.pred,
      locations = grid,
      axes = FALSE,
      xlab = "", ylab = "",
      breaks = seq(25, 65, 3.1),
      col = heat.colors(12),
      x.leg = c(-2329176, -1000000),
      y.leg = c(-1500000, -1200000))
title("Ordinary Kriged Predictions of Ozone", line = -1.5)
points(ozone.geo$coords, pch = 16, cex = .5)
plot(st_geometry(States), add = TRUE)

image(x = OK.pred,
      locations = grid,
      axes = FALSE,
      xlab = "", ylab = "",
      value = sqrt(OK.pred$krige.var),
      breaks = seq(3.5, 8, .35),
      col = heat.colors(12),
      x.leg = c(-2329176, -1000000),
      y.leg = c(-1500000, -1200000))
title("Ordinary Kriged Predictions Standard Errors", line = -1.5)
points(ozone.geo$coords, pch = 16, cex = .5)  # remember this version should work too
plot(st_geometry(States), add = TRUE)

image(x = UK.pred,
      locations = grid,
      axes = FALSE,
      xlab = "", ylab = "",
      breaks = seq(25, 65, 3.1),
      col = heat.colors(12),
      x.leg = c(-2329176, -1000000),
      y.leg = c(-1500000, -1200000))
title("Universal Kriged Predictions of Ozone", line = -1.5)
points(ozone.geo$coords, pch = 16, cex = .5)
plot(st_geometry(States), add = TRUE)

image(x = UK.pred,
      locations = grid,
      axes = FALSE,
      xlab = "", ylab = "",
      value = sqrt(UK.pred$krige.var),
      breaks = seq(3.5, 8, .35),
      col = heat.colors(12),
      x.leg = c(-2329176, -1000000),
      y.leg = c(-1500000, -1200000))
title("Universal Kriged Predictions Standard Errors", line = -1.5)
points(ozone.geo$coords, pch = 16, cex = .5)
plot(st_geometry(States), add = TRUE)

### because in this instance we are plotting multiple types of objects at once,
# R plot can get confused overlaying the layers together. If you turn off the 
# device (dev.off()) or close your work and come back before running one of the 
# image commands again, you'll need to re-plot the sf object first to establish
# the correct boundaries of the plotting window before running image() again, as done in line 742



### If you are doing more than taking a first look at your results in R,
#' you may want to save the estimates and SEs of the different models 
#' for additional summary, comparison, analysis into a single file which
#' you can work with again in R or in ArcGIS later. We can do this in two ways:

#' First, store the predictions back into a flat file (like a csv)

Data2Export <- data.frame(x = grid$easting,
                          y = grid$northing,
                          IDWpreds = IDW.pred$var1.pred,
                          TSpreds = TS.pred$fit,
                          TSse = sqrt(TS.pred$se.fit),
                          OKpreds = OK.pred$predict,
                          OKse = sqrt(OK.pred$krige.var),
                          UKpreds = UK.pred$predict,
                          UKse = sqrt(UK.pred$krige.var))
write.csv(Data2Export, "Results_ModelComparison.csv")


#' Second, save the information into a spatial object
#' estimates (and errors) are at the point locations of the prediction grid

# the attribute data is the same as the info we just stored in the flat file
# and now we are just specifying the geometry information to make it an SF object
# give it points of the grid (where we made our predictions) + the coordinate system of the ozone data

sf2Export <- st_as_sf(Data2Export, coords = c("x", "y"), crs = st_crs(Ozone))
# write the spatial data file (shp), saving to default location of working directory
st_write(sf2Export, "Results_ModelComparison.shp")

# Note that this is only points now, not a filled in grid (raster or pixel data)
# In ArcGIS you could now use "out of the box" interpolation tools to map the complete 
# area (and this estimates will be better than using the IDW from the raw ozone data because
# you've already generated more accurate predictions to a fine grid with kriging)



## Depending how you will utilize your results in future work, 
#' you may want to reference the above options for storing and visualizing data
#' However, for this assignment, we have provided all the material you need and
#' the plots from the image commands are sufficient to display your results,  
#' though you might consider resizing your plot dimensions to minimize white space
#' and make sure graphics are easily legible however you choose to export to turn
#' in your answers (exporting plot, taking a screenshot, etc.)

 

#-----------------------------------------
#
#Extra Code: MLE estimation and Kriging
#
#-----------------------------------------


# The following performs maximum likelihood estimation of the semivariogram parameters and regression betas
# and compares to a trend surface model and semivariograms estimated via weighted-least squares.
# This may take a little while to run (1-2 minutes). Note that MLE won't always fit the variogram plot as nicely as
# estimation based on least squares (OLS or WLS).  Make sure you understand why.

### MLE estimation

ozone.vario.mle <- likfit(ozone.geo,
                          ini.cov.pars = c(33.51, 1083000),
                          cov.model = "spherical",
                          nugget = 18.28)

ozone.vario.mle   # Output shows regression estimates of intercept (beta) and covariance parameters

# A large-scale parameter model summary (for example if running analysis for inference via GLS)
cbind(Estimate = ozone.vario.mle$beta,
      S.E. = sqrt(ozone.vario.mle$beta.var),
      t.value = ozone.vario.mle$beta / sqrt(ozone.vario.mle$beta.var),
      p.value = 2 * pt(-abs(ozone.vario.mle$beta / sqrt(ozone.vario.mle$beta.var)),
                       df = nrow(data2use) - ozone.vario.mle$npars))

# Let's compare to the trend surface (OLS) model
summary(model0)
# Notice how the SE has increased over 8x from the trend surface model! The actual estimate only slightly changes

# Estimate for Universal Kriging
ozone.vario.UK.mle <- likfit(ozone.geo,
                             ini.cov.pars = c(15.44, 742000),
                             cov.model = "spherical",
                             nugget = 20.96,
                             trend = trend.spatial( ~ ns(easting, 4), ozone.geo))

ozone.vario.UK.mle #output shows regression estimates of intercept (beta),
# spline parameters (beta1 - beta4) and covariance parameters

# Regression model summary
cbind(Estimate = ozone.vario.UK.mle$beta,
      S.E. = sqrt(diag(ozone.vario.UK.mle$beta.var)),
      t.value = ozone.vario.UK.mle$beta / sqrt(diag(ozone.vario.UK.mle$beta.var)),
      p.value = 2 * pt(-abs(ozone.vario.UK.mle$beta / sqrt(diag(ozone.vario.UK.mle$beta.var))),
                       df = nrow(data2use) - ozone.vario.UK.mle$npars))

# Let's compare to the trend surface (OLS) model again
summary(model1)
# Again, compared to the trend surface model, standard errors for all betas are much higher, and some spline
# estimates change considerably. The 1st and 4th spline are no longer statistically significant at alpha=0.05

# As shown below I had to increase the range of the y-axis to make the MLE semivariogram fit in the plot.
plot(ozone.vario1,
     pch = 16,
     xlab = "Distance (meters)",
     ylab = "Semivariogram",
     ylim = c(0, 65))
title("Fitted WLS and MLE Semivariograms for Ozone levels")
lines(ozone.vario.wls, col = "blue")
lines(ozone.vario.mle, col = "red")
legend(locator(1),
       legend = c("WLS", "MLE"),
       lwd = c(1, 1),
       col = c("blue", "red"))

# Make sure you remember to place your legend if 'legend(locator(1))' is used

plot(ozone.resid.vario,
     pch = 16,
     xlab = "Distance (meters)",
     ylab = "Semivariogram",
     ylim = c(0, 65))
title("Fitted WLS and MLE Semivariograms for Ozone Residuals")
lines(ozone.resid.vario.wls, col = "blue")
lines(ozone.vario.UK.mle, col = "red")
legend(locator(1),
       legend = c("WLS", "MLE"),
       lwd = c(1, 1),
       col = c("blue", "red"))


### MLE Kriging onto grid

OK.pred.mle <- krige.conv(ozone.geo,
                          locations = grid,
                          krige = krige.control(obj.model = ozone.vario.mle))

UK.pred.mle <- krige.conv(ozone.geo,
                          locations = grid,
                          krige = krige.control(obj.model = ozone.vario.UK.mle,
                                                trend.d = trend.spatial( ~ ns(easting, 4), ozone.geo),
                                                trend.l = trend.spatial( ~ ns(easting, 4), grid)))
plot(st_geometry(States))
image(x = OK.pred.mle,
      locations = grid,
      axes = FALSE,
      xlab = "", ylab = "",
      breaks = seq(25, 65, 3.1),
      col = heat.colors(12),
      x.leg = c(-2329176, -1000000),
      y.leg = c(-1500000, -1200000))
title("Ordinary Kriged (MLE) Predictions of Ozone", line = -1.5)
points(ozone.geo$coords, pch = 16, cex = .5)
plot(st_geometry(States), add = TRUE)

image(x = OK.pred.mle,
      locations = grid,
      axes = FALSE,
      xlab = "", ylab = "",
      value = sqrt(OK.pred.mle$krige.var),
      breaks = seq(3.5, 8, .35),
      col = heat.colors(12),
      x.leg = c(-2329176, -1000000),
      y.leg = c(-1500000, -1200000))
title("Ordinary Kriged (MLE) Predictions Standard Errors", line = -1.5)
points(ozone.geo$coords, pch = 16, cex = .5)
plot(st_geometry(States), add = TRUE)

image(x = UK.pred.mle,
      locations = grid,
      axes = FALSE,
      xlab = "", ylab = "",
      breaks = seq(25, 65, 3.1),
      col = heat.colors(12),
      x.leg = c(-2329176, -1000000),
      y.leg = c(-1500000, -1200000))
title("Universal Kriged (MLE) Predictions of Ozone", line = -1.5)
points(ozone.geo$coords, pch = 16, cex = .5)
plot(st_geometry(States), add = TRUE)

image(x = UK.pred.mle,
      locations = grid,
      axes = FALSE,
      xlab = "", ylab = "",
      value = sqrt(UK.pred.mle$krige.var),
      breaks = seq(3.5, 8, .35),
      col = heat.colors(12),
      x.leg = c(-2329176, -1000000),
      y.leg = c(-1500000, -1200000))
title("Universal Kriged (MLE) Predictions Standard Errors", line = -1.5)
points(ozone.geo$coords, pch = 16, cex = .5)
plot(st_geometry(States), add = TRUE)

# Look at your individual files to compare WLS & MLE Kriged maps & discuss


