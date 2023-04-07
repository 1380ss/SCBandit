#install.packages("googledrive")
library("googledrive")
drive_auth()


path <- '~/Simulations/PostDiff/April6'


drive_upload('Rplot4.png',path=path,"1.png")
  
getwd()
