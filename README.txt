DosageConvertor is a C++ tool to convert dosage files (in VCF format) from Minimac3 to ther formats such as MaCH or PLINK.

<<< SEE http://genome.sph.umich.edu/wiki/DosageConvertor FOR DOCUMENTATION >>>

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
