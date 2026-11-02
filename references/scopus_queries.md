# Scopus queries: Third Pole, Indus, snow, glaciers, monsoon, ENSO

Paste into Scopus Advanced Search. Not yet run in Scopus: export results as RIS and CSV into this folder, and record the date and count below.

## A. Snow cover across the region
```
TITLE-ABS-KEY ( ( "snow cover" OR "snow extent" OR "snow duration" OR "snow phenology" OR "snow line" OR snowline OR "snow water equivalent" )
AND ( "Third Pole" OR "High Mountain Asia" OR "Hindu Kush" OR Karakoram OR Himalaya* OR Pamir OR "Tian Shan" OR "Tien Shan" OR "Tibetan Plateau" OR Hengduan OR Indus OR Ganges OR Brahmaputra )
AND ( MODIS OR AVHRR OR Landsat OR Sentinel OR satellite OR "remote sensing" OR reanalysis ) )
AND PUBYEAR > 1999
```

## B. Snow and climate drivers
```
TITLE-ABS-KEY ( ( snow* )
AND ( "Third Pole" OR "High Mountain Asia" OR "Hindu Kush" OR Karakoram OR Himalaya* OR Pamir OR "Tian Shan" OR "Tibetan Plateau" OR Indus )
AND ( monsoon OR ENSO OR "El Nino" OR "La Nina" OR "Indian Ocean Dipole" OR "western disturbance*" OR westerlies OR "North Atlantic Oscillation" OR teleconnection* ) )
```

## C. Glaciers and water
```
TITLE-ABS-KEY ( ( glacier* AND ( "mass balance" OR "mass change" OR "Karakoram anomaly" OR runoff OR meltwater ) )
AND ( "Third Pole" OR "High Mountain Asia" OR "Hindu Kush" OR Karakoram OR Himalaya* OR Pamir OR "Tian Shan" OR "Tibetan Plateau" OR Indus ) )
AND PUBYEAR > 2004
```

## D. Testing for gaps in rigour
```
TITLE-ABS-KEY ( ( "snow cover" OR "snow duration" ) AND ( "High Mountain Asia" OR Himalaya* OR Karakoram OR "Hindu Kush" OR "Tibetan Plateau" OR Indus )
AND ( "false discovery" OR "field significance" OR "multiple testing" OR "cloud gap" OR "sensitivity analysis" OR "spatial autocorrelation" ) )
```

## E. Long records, before 2000
```
TITLE-ABS-KEY ( ( "snow cover" OR snowline ) AND ( AVHRR OR "passive microwave" OR Landsat OR "long-term" OR Hexagon ) AND ( "Tibetan Plateau" OR Himalaya* OR Karakoram OR "High Mountain Asia" OR Indus ) )
```

## Mangrove and GEDI
```
TITLE-ABS-KEY ( ( mangrove* OR "Indus Delta" OR "Avicennia marina" OR "Rhizophora mucronata" OR "blue carbon" )
AND ( "canopy height" OR "tree height" OR "forest height" OR CHM OR "aboveground biomass" OR "above-ground biomass" OR AGB OR AGBD OR "carbon stock*" OR allometr* )
AND ( GEDI OR "Global Ecosystem Dynamics Investigation" OR ICESat* OR LiDAR OR "remote sensing" OR satellite OR Sentinel* OR "ALOS PALSAR" OR "TanDEM-X" OR SAR OR "Google Earth Engine" OR "machine learning" OR "random forest" OR "deep learning" ) )
AND PUBYEAR > 2009
```
```
TITLE-ABS-KEY ( ( "Indus Delta" OR "Indus River Delta" OR "Keti Bandar" OR Sindh OR Pakistan ) AND mangrove* )
```
```
TITLE-ABS-KEY ( ( GEDI OR "Global Ecosystem Dynamics Investigation" ) AND ( "canopy height" OR biomass ) AND ( "wall-to-wall" OR map* OR fusion OR "spatial cross-validation" OR validation ) ) AND PUBYEAR > 2018
```
