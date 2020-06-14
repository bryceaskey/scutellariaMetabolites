# scutellariaMetabolites
The scutellariaMetabolites project is an investigation of the remarkable flavonoid diversity present across the hundreds of species of plants in the *Scutellaria* genus. This project is a collaborative effect between researchers from the Horticultural Sciences Department at the University of Florida, the Department of Biology at Sungshin Women's Univeristy, and the Gulf Coast Research and Education Center at the Univeristy of Florida.

A brief summary of the goals of this project is as follows:
1. Chemically assay a subset of *Scutellaria* species to identify those accumulating high concentrations of medicinally relevant flavonoids. 
2. Combine chemical and phylogenetic data to provide insight into the evolutionary history of flavonoid diversity, and identify "medicinal hotspots" in the genus.
3. Perform organ-specific chemical profiling to determine differences in site of flavonoid accumulation across multiple species.
4. Apply transcriptomic profiling to identify differentially expressed genes potentially responsible for organ-specific differences in flavonoid accumulation.

Overall, a better understanding of flavonoid diversity and the mechanisms controlling it will facilitate the development of drugs and treatments from *Scutellaria* plants.

The purpose of this repository is to serve as a bookkeeping tool for all of code, data, and figures generated for each part of this project.

# Bioinformatics analysis
The documentation in this section details the bioinformatics approaches and methods applied in this project. Ideally, the documentation here should allow someone with little to no bioinformatics experience to understand and recreate the analysis method that was used. 

Topics covered:
1. [Setting up a hipergator account](#setting-up-hipergator)
2. [Submitting jobs in hipergator](#submitting-jobs) \
   2.1 [SLURM scripts](#slurm-scripts) \
   2.2 [Bash basics](#bash-basics) \
   2.3 [Output files](#output-files)
3. [Downloading datasets](#downloading-datasets)

Note that *Scutellaria* plants are "non-model" organisms, which means that there are a couple of extra processing and error checking steps necessary as compared to an analysis workflow that would be applied to data from a "model" organism, such as *Arabidopsis thaliana*.

## 1. Setting up a hipergator account
<a name="setting-up-hipergator"></a>
The ufrc wiki provides a [pretty good guide](https://help.rc.ufl.edu/doc/Getting_Started) on how to set up a hipergator account. The first step is to request an account from your supervisor. Assuming that your supervisor has already set up a group with allocated resources in hipergator, just submit a [request account form](https://www.rc.ufl.edu/access/request-account/). Once your supervisor approves the account request, you should receive an email from ufrc stating that your account has been created. 

Once your account has been created, the next step is to connect to hipergator using a secure shell, or SSH. An SSH is what you will use to communicate with hipergator (e.g. submit jobs, check the status of jobs). This process will differ depending on whether your computer is running Windows, Linux, or MacOS. 

**Windows users -** Since Windows doesn't come with an SSH preinstalled, you will need to download one. [MobaXTerm](https://mobaxterm.mobatek.net/) has a pretty friendly user interface, so I recommend it for newer users. [This video](https://mediasite.video.ufl.edu/Mediasite/Play/2bf4c860f19b48a593fb581018b813a11d) from ufrc provides a good tutorial for first time setup of MobaXTerm. In short, after opening MobaXTerm and clicking the "New session" button, select the "SSH" option at the top of the popup window. Then enter the name of the remote host (hpg.rc.ufl.edu), check the "Specify username" box, enter your gatorlink username, and click the "OK" button.

![Image of MobaXTerm SSH setup](https://github.com/bryceaskey/scutellariaMetabolites/tree/master/figures/docImages/MobaXTerm.png)

*An image of what MobaXTerm SSH setup would look like for my account (with Gatorlink username braskey).*

**Linux and MacOS users -** Both Linux and MacOS systems come with an SSH by default, which can be accessed from 

## 2. Submitting jobs in hipergator
<a name="submitting-jobs"></a>

### 2.1 SLURM scipts
<a name="slurm-scripts"></a>

### 2.2 Bash basics
<a name="bash-basics"></a>

### 2.3 Output files
<a name="output-files"></a>

## 3. Downloading datasets
<a name="downloading-datasets"></a>