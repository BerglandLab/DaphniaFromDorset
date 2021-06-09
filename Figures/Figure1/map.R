library (leaflet)
library( magrittr)
library(mapview)
#
# In Google Maps, find a location, right click, select 'What's Here', and copy
# the lat/long values. These values can be stored in a .csv file. I've added them
# here to file for simplicity.
#

df <- read.csv(textConnection(
"Name,Lat,Long
Kilwood,50.643015 ,-2.08956
D10,50.709379,-2.206421
"))

p1 <- leaflet(height=750, width=550) %>%
addProviderTiles(providers$OpenStreetMap.HOT) %>%
addMarkers(., lat=df$Lat, lng=df$Long, label=df$Name) %>%
setView(., -2.08956/2 + -2.206421/2 + 0.012, 50.643015/2 + 50.709379/2- .033, zoom = 12) %>%
addScaleBar()

mapshot(p1, file = "~/Rplot.pdf", selfcontained=F)
