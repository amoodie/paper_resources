# Data organization

The `toSelectedPlots.mat` file which contains all the data has many fields. This document attempts to describe what _many_ of them are.


## the structure `df`

Each row represents a station, where we stopped to collect data. Each station contains the following information.


### StationID

A nominal station name used to label bottles and refer to it in the field, etc. Includes date information.


### NominalLocation

Name of nearby city.


### Transect

Which transect the station lies on.


### WrtChannel

Where along the transect the station is located. LB = left bank, C = center, RB = right bank.


### StationLocation

The name of the station without date information.


### CollectionDate

Date station visited


### FlowDepthAtCollection

Flow depth at the collection time.


### Samples

Raw-ish information about all the samples collected. See __waterSamplesTable__ below for more refined data. Contains the fields:

* __FullSampleName:__ full name given in field at birth
* __SampleID:__ name refined for data analysis
* __SampleLocation:__ station name of sample
* __WrtChannel:__ where in channel
* __NominalDepth:__ fraction depth from surface
* __NominalZ:__ fraction height above bed
* __IterationNumber:__ count for multiple samples at depth
* __FlowDepthAtCollection:__ flow depth at collection time and location
* __PredictedCollectionDepth:__ product of __NominalDepth__ and __FlowDepthAtCollection__, where sample was intended to be taken
* __ActualCollectionDepth:__ where the sample was actually taken
* __MetersAboveBed:__ meters above the bed of sample location
* __TotalVolume_mL_:__ volume of water and sediment collected
* __TotalMass_g:__ mass of sediment collected
* __RawConcentration_gl:__ concentration of sample
* __conc:__ DOAA (?)
* __concData:__ some more crap relating to the filtering and processing of samples to get concentration data
* __gsSummRawF:__ summary of raw grain size analysis of all sample runs sample
* __gsDistRawF:__ full raw grain size analysis of all sample runs sample
* __gsDistMeanF:__ average distribution of multiple analysis runs
* __gsDistMeanFNum:__ DOAA as a `double` type
* __gsSummMeanF:__ summary of mean distribution
* __gsWLdata:__ calculation data for determining the washload cutoff
* __gsDistMeanNW:__ average distribution with washload subtracted
* __gsDistMeanNWnorm:__ DOAA normalized to sum to 100%
* __concNW:__ concentration without washload
* __gsSummMeanNWnorm:__ summary of DOAA


### Velocity

* __meter:__ data from the mechanical velocimeter
* __adcp:__ data from the aDcp
* __meterTable:__ formatted to a table
* __adcpTable:__ formatted to a table
* __flowDepth:__ depth at the time of measurement
* __slope:__ slope of the water surface
* __Cf:__ estimated friction coefficient
* __taubDSP:__ estimated depth slope product tau_b
* __ustarDSP:__ estimated u_* from depth slope product
* __Ubar_conCf:__ estimated depth averaged velocity
* __Fr:__ estimated Froude number
* __tstar:__ estimated Shields number
* __tstarsk:__ estimated Shields number going to skin friction
* __ustarsk:__ estimated u_* going to skin friction
* __ustarBest:__ subjective selection based on known station information

### waterSamplesTable

Organized table of the suspended sediment samples collected at that station. This table contains fields:

* __sampleDepth:__ fraction depth from surface
* __sampleZ:__ meters above bed
* __sampleZnorm:__ fraction height above bed
* __sampleIter:__ count for multiple samples at depth
* __conc:__ concentration of sediment
* __gsWLpercent:__ percent of sample that is washload
* __concNW:__ concentration with washload subtracted
* __gsClass:__ grain size class bins
* __gsDistNW:__ grain size distribution without washload
* __gsDistNWnorm:__ DOAA normalized to sum to 100%
* __gsSummNWnorm:__ summary of DOAA
* __concNWbyClass:__ per grain size class concentration without washload


### gsDistBed

Grain size distribution of the bed sample.


### gsSummBed

Summary of the bed sample.


### gsDistNearBedF

Grain size distribution of a near-bottom sample, a mean value if available.


### gsDistNearBedNW

DOAA without washload.


### gsDistNearBedNWnorm

DOAA normalized to sum to 100% of a grain size distribution again.


### gsSummNearBedNWnorm

Summary of DOAA.