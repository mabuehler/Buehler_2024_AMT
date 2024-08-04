# Data to Bühler et al. (2024) #
Data to the publication Bühler et al. (2024) - Applicability of the inverse dispersion method to measure emissions from animal housing (https://doi.org/10.5194/amt-2023-258)

## Data availability ##
Provided are the raw data of the instruments, scripts to treat them and reproduce the findings in the publication, R outputs and Figures.


## Scripts ##
In total, there are 10 scripts provided, of which most of them are needed to reproduce the data in the publication.

Below, a brief explanation of the content of the different scripts.

- 01_Datatreatment_01_Weatherstation.r <br /> # This script reads in the weather station data and makes it ready for further use.
- 01_Datatreatment_02_Sonics.r <br /> # This script reads in the 3D ultrasonic data and makes it ready for further use.
- 01_Datatreatment_03_GasFinder.r <br /> # This script reads in the GasFinder data and makes it ready for further use.
- 01_Datatreatment_04_MFC_Pressuresensor.r <br /> # This script reads in the Mass flow controller (MFC) and pressure sensor data and makes it ready for further use.
- 02_Calculation_01_bLS.r <br /> # This script is made to run the bLSmodelR on the number cruncher of the University of Applied Sciences BFH. The code should also work on your computer but you have to adopt the number of cores.
- 02_Calculation_02_Concentration.r <br /> # his script treats the unprocessed concentration data. It removes false concentrations and applies an intercalibration and makes the data ready for further use.
- 02_Calculation_03_Emissions.r <br /> # This script calculates emissions and makes it ready for further use.
- 02_Calculation_04_contourXYZ_Plume.r <br /> # This script calculates the plume contours in the XY and XZ plane. This script is not necessary to reproduce the findings of the publication.
- 03_Apply_filter.r <br /> # This script applies the quality filtering and makes the data ready for further use.
- 04_Plots_Tables.r <br /> # With this scirpt one can create all the plots and values of the tables in the publication, the supplement and the initial submission.

Note, for the geometry, there is no script available. The coordinates of the different sensors and the source is provided as R output.


## Naming of instruments ##
The instruments in the publication have different names than in the scripts. In some scripts the final names are also provided but throughout the evaluation the original device names are used. Only in the script 04_Plots_Tables.r are the final names introduces. Below an overview of what original name corresponds to the final name of the devices:

#### GasFinder instruments called 'OP' in the publication
- OP-UW = GF26
- OP-2.0h = GF17
- OP-5.3h = GF18
- OP-6.8h = GF16
- OP-12h = GF25

#### 3D ultrasonic anemometer instruments called 'UA' in the publication
- UA-UW = SonicC
- UA-2.0h = SonicA
- UA-5.3h = Sonic2
- UA-6.8h = SonicB


#### Source
In some of the scripts, the source might be called 'Schopf' which is a local term for 'shed'.


## Note ## 
This code was written by Marcel Bühler (minor parts might originally be written by Christoph Häni) and is intended to follow the publication 'Applicability of the inverse dispersion method to measure emissions from animal housing' in AMT.
Please feel free to use and modify it (e.g., use it to run different dispersion models), but attribution is appreciated.


## Disclaimer ##
I do not guarantee that everything works. It might be that not all variables were changed to English for a better understanding correctly.
Unfortunately, it is not possible to provide all the catalogs of the bLS run, as the total size is several 100s of GB. In case you run the bLS model on your own, the result will have a minimal difference, as no bLS run produces twice the same result. This should, however, not alter the findings.


## Contact ##
In case you have question, please contact Marcel Bühler (mb@bce.au.dk). In case this does not work, Christoph Häni might also be able to help (christoph.haeni@bfh.ch)
