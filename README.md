# PCB-Reactive-Transport-Model
The Polychlorinated Biphenyl (PCB) Reactive Transport Model produces simulations of congener-specific PCB concentration in the aqueous and gas phase of laboratory-scale bioreactors containing PCB-contaminated sediment.  The script uses experimental passive sampler measurements and a minimization of least squares procedure to fit congener-specific sampling rates for both passive samplers. Finally, the model can simulate the reduction in volatile release of PCBs from sediment to air resulting from the presence of an aerobic PCB-degrading microorganism with a known PCB biotransformation rate. 

This ReadMe file was generated on 2022-12-21 by Andres Martinez and last updated by CMB on 2022-02-08

----------------------
General Information
----------------------

Deposit Title: Polychlorinated Biphenyl (PCB) Reactive Transport Model

Contributor information:

Andres J. Martinez Araneda, PhD
University of Iowa - Department of Civil & Environmental Engineering
Iowa Superfund Research Program (ISRP)
andres-martinez@uiowa.edu
ORCID: 0000-0002-0572-1494

Principal Investigator: Timothy Mattes, PhD
Principal Investigator email: tim-mattes@uiowa.edu

Date of Data Collection:

This work was supported by the National Institutes of Environmental Health Sciences (NIEHS) grant #P42ES013661.  The funding sponsor did not have any role in study design; in collection, analysis, and/or interpretation of data; in creation of the dataset; and/or in the decision to submit this data for publication or deposit it in a repository.

Subject: Polychlorinated Biphenyls; Contaminant fate and transport; Paraburkholderia xenovorans LB400; Kinetic phase passive sampling; Bioremediation; Biodegradation; Biosurfactants; Bioavailability; GC-MS/MS

GeoLocation: The sediment used in this study was taken from a PCB-contaminated emergency overflow lagoon located in Altavista, VA (37°06'52"N, 79°16'21"W). Laboratory and analytical work was done at the University of Iowa in Iowa City, IA, USA.

--------------------------
SHARING/ACCESS/ATTRIBUTION LICENSE INFORMATION
--------------------------

Licenses/restrictions placed on the model documentation, or limitations of reuse: Open Data Commons Attribution License (ODC-By)

Recommended citation for the model documentation:

Bako, Christian M.; Martinez, Andres J. (2022): Polychlorinated Biphenyl (PCB) Reactive Transport Model. University of Iowa. (dataset). https://doi.org/reservedDOI

Citation for and links to publications that cite or use the model:



--------
ABSTRACT & PROJECT DESCRIPTION
--------

This Iowa Research Online (IRO) Data Repository deposit describes documentation of code developed in "R" software for a "Polychlorinated Biphenyl (PCB) reactive transport model". The PCB Reactive Transport Model produces simulations of congener-specific PCB concentration in the aqueous and gas phase of laboratory-scale bioreactors containing PCB-contaminated sediment slurry.  The script uses experimental passive sampler measurements and minimization of least squares procedures to fit congener-specific adsorption and desorption coefficients to the model.  Finally, the model can simulate the reduction in volatile release of PCBs from sediment slurry to air resulting from the presence of an aerobic PCB-degrading microorganism with a known PCB biotransformation rate. As written, the model script calls data resulting from bioreactor experiments described in the peer-reviewed journal article ("Aerobic Bioaugmentation to Decrease Polychlorinated Biphenyl (PCB) Emissions from Contaminated Sediments to Air") linked to this IRO deposit. In this deposit, the full "R" script is available for PCB 52 to provide documentation for the underlying code and to provide a framework for future re-use. The code is documented and described in detail within an "R Markdown" (.RMD) file which can be viewed as .html, .pdf, or .docx. The purpose of this IRO deposit is to promote reuse of the underlying experimental dataset from the linked publication and the PCB reactive transport model.  We hope other researches will use this model to conduct inform other analyses, investigations, and / or interpretive insights.  It is recommended that a user first reads the linked publication and one of the R Markdown output files to interpret and contexualize the model script.

--------------------
IOWA RESEARCH ONLINE DEPOSIT & FILE OVERVIEW
--------------------
File Name: 00_PCB_ReactiveTransportModel_ReadMe
File Size: 7.35kb
File Format: .txt
File description: Contains a brief description of the project, a description of files also associated with the IRO deposit, prerequisite software, installation instructions for software, and user contribution information.

File Name: 01_PCB_ReactiveTransportModel_R
File Size: kb
File Format: .R
File description: Contains the raw code for the PCB Reactive Transport model in "R" which uses PCB 52 data from xxx.csv to compare simulated results to experimental observations.

--------
PREREQUISITES & DEPENDENCIES
--------
This section of the ReadMe file lists the necessary software and packages required to run the PCB reactive transport model in "R" and to view documentation provided in the "R Markdown" file, as well as in the various output files (e.g., .html, .pdf, and .docx).

For more information on "R Markdown" documentation please see: https://rmarkdown.rstudio.com/index.html

Software:
- Any web browser (e.g., Google Chrome, Microsoft Edge, Mozilla Firefox, etc.)
- R-studio for easily viewing, editing, and executing "R" code as a regular "R script" file and as an "R markdown" file: https://www.rstudio.com/products/rstudio/download/
- MiKTeX for outputting "R Markdown" file as a PDF: https://miktex.org/download
- Adobe Acrobat Reader for viewing PDF output of "R Markdown" file: https://get.adobe.com/reader/otherversions/
- Word for viewing .docx output of "R Markdown" file (requires software license)

R Packages:
- dplyr: Package used for reading and formatting data from files - https://dplyr.tidyverse.org/
- ggplot2: Package used for plotting data - https://ggplot2.tidyverse.org/
- reshape2: Package used for reshaping data (tall-narrow <-> short-wide) - https://cran.r-project.org/web/packages/reshape2/index.html
- deSolve: Package used for solving differential equations - http://desolve.r-forge.r-project.org/
- minpack.lm: Package for least squares fitting procedure with the levenberg-marquart algorithm - https://rdrr.io/cran/minpack.lm/

Other:
- Environmental Science & Technology citation style file (.csl) for including references in "R Markdown" file: https://www.zotero.org/styles

--------
SOFTWARE INSTALLATION
--------

This section of the ReadMe file provides short instructions on how to download and install "R Studio".  "R Studio" is an open source (no product license required) integrated development environment (IDE) for "R" and completely free to use.  To install "R Studio" follow the instructions below:

1. Visit the following web address: https://www.rstudio.com/products/rstudio/download/
2. Click the "download" button beneath RStudio Desktop
3. Click the button beneath "Download RStudio Desktop".  This will download the correct installation file based on the operating system detected.
4. Run the installation file and follow on-screen instructions.  

--------
USER CONTRIBUTION & ISSUE TRACKING
--------
https://github.com/valdiman/PCB-Reactive-Transport-Model

--------
METADATA
--------
See the "XXX.json" file for all available CodeMeta metadata.

