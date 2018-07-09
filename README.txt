DosageConvertor is a C++ tool to convert dosage files (in VCF format) from Minimac3/4 to ther formats such as MaCH or PLINK.

<<< SEE http://genome.sph.umich.edu/wiki/DosageConvertor FOR DOCUMENTATION >>>

 To install, type the following in the main folder: bash install.sh
 
Users should follow the following steps to compile DosageConvertor 

## Prerequisites

Automatic installation of DosageConvertor requires [cget](http://cget.readthedocs.io/en/latest/src/intro.html#installing-cget) and cmake v3.2. These prerequisites can be installed as follows:

Ubuntu 16.04
```
sudo apt-get install cmake python-pip python-dev
pip install cget
```
Ubuntu 14.04
```
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:george-edison55/cmake-3.x
sudo apt-get update
sudo apt-get install cmake python-pip python-dev
pip install cget
```
MacOS
```
brew install cmake
sudo easy-install pip
pip install --user cget --ignore-installed six
```

## Installation
The easiest way to install DosageConvertor and its dependencies is to use the install.sh file provided.

```
cd DosageConvertor
bash install.sh
```

Alternatively, you can setup a dev environment cmake directly.
```bash
cd DosageConvertor
cget install -f ./requirements.txt                      # Install dependencies locally.
mkdir build && cd build                                 # Create out of source build directory.
cmake -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake .. # Configure project with dependency paths.
make                                                    # Build.
```



## Usage


 Usage: ./DosageConvertor  --vcfDose      TestDataImputedVCF.dose.vcf.gz
                           --info         TestDataImputedVCF.info
                           --prefix       OutputFilePrefix
                           --type         plink OR mach   // depending on output format
                           --format       DS or GP        // based on if you want to output
                                                          // dosage (DS) or genotype prob (GP)
                           --buffer       10000           // Number of Markers to import and
                                                          // print at a time (valid only for
                                                          // MaCH format)
                           --idDelimiter  _               // Delimiter to Split VCF Sample ID into
                                                          // FID and IID for PLINK format

 
<<< SEE http://genome.sph.umich.edu/wiki/DosageConvertor FOR DOCUMENTATION >>>
