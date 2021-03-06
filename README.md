# StatEngine
This package provides R tools for methods covered in the course STAT 509 "Statistics for Engineers" at the University of South Carolina.

Prior to installation of this R package, you need to install R and Rstudio first.

# Installation via GitHub

You need to install the R package "devtools" first.

    install.packages("devtools")
    
Then run these two lines to install the StatEngine package.

    library(devtools)
    install_github("Harrindy/StatEngine",force=TRUE) 
    #Press Enter if needed to finish the installation.
    
Run this to check whether the installation is successful:

    library(StatEngine)
   
This command loads the StatEngine Package. If you get 

    Loading required package: car
    Loading required package: carData   

it means you did it!



# Installation from a local file

You need to download the file StatEngine_1.0.0.tar.gz to your own local drive (remember the path).

Then install the R package car first.

    install.packages("car")
    
Now follow the general instruction below to install the StatEngine package from the downloaded file.
    
    install.packages(path_to_file, repos = NULL, type="source")
    
Where path_to_file would represent the full path and file name:

On Windows it will look something like this: 

    "C:\\StatEngine_1.0.0.tar.gz"
    
For example, if you download StatEngine_1.0.0.tar.gz to your C drive (in no folder). Then you should use 

    install.packages("C:\\StatEngine_1.0.0.tar.gz", repos = NULL, type="source")

On UNIX/MAC it will look like this: 

    "/home/blah/StatEngine_1.0.0.tar.gz"
    
For example, I saved StatEngine_1.0.0.tar.gz in the Downloads folder of my Mac, I then use

    install.packages("/Users/harrindy/Downloads/StatEngine_1.0.0.tar.gz", repos = NULL, type="source")
    
After installation, run this to check whether it is successful:

    library(StatEngine)
   
This command loads the StatEngine Package. If you get 

    Loading required package: car
    Loading required package: carData   

it means you did it!

