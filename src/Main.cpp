

#include <iostream>
#include <ctime>
#include <stdio.h>
#include "Parameters.h"
#include "StringBasics.h"
#include "HaplotypeSet.h"

using namespace std;
void DosageConvertorVersion();
void helpFile();

int main(int argc, char ** argv)
{
	// Parameter Options

	bool   help = false, params = false;
    HaplotypeSet inputDose;
    bool nobgzip=false;

	ParameterList inputParameters;
	PhoneHome::allThinning = 50;

	BEGIN_LONG_PARAMETERS(longParameterList)
		LONG_PARAMETER_GROUP("Input Dosage")
		LONG_STRINGPARAMETER("vcfDose", &inputDose.VcfDose)
		LONG_STRINGPARAMETER("info", &inputDose.Info)
		LONG_STRINGPARAMETER("tag", &inputDose.tag)
		LONG_PARAMETER_GROUP("Output Parameters")
		LONG_STRINGPARAMETER("prefix", &inputDose.outFile)
		LONG_STRINGPARAMETER("format", &inputDose.Format)
		LONG_STRINGPARAMETER("type", &inputDose.Type)
		LONG_PARAMETER("nobgzip", &nobgzip)
		LONG_PARAMETER_GROUP("Other Parameters")
		LONG_INTPARAMETER("buffer", &inputDose.BufferSize)
		LONG_PARAMETER("allDiploid", &inputDose.samePloidy)
		LONG_PARAMETER("trimNames", &inputDose.TrimAlleles)
		LONG_INTPARAMETER("trimLength", &inputDose.TrimLength)
		LONG_STRINGPARAMETER("SexFile", &inputDose.SexFile)
		LONG_STRINGPARAMETER("idDelimiter", &inputDose.IdDelimiter)
		LONG_PARAMETER("help", &help)
		LONG_PARAMETER("params", &params)
		LONG_PHONEHOME(VERSION)
		BEGIN_LEGACY_PARAMETERS()
		LONG_STRINGPARAMETER("MyChromosome", &inputDose.MyChromosome)
//		LONG_STRINGPARAMETER("golden", &golden)
//		LONG_INTPARAMETER("sample", &max_indiv)
//		LONG_INTPARAMETER("marker", &max_marker)
//		LONG_STRINGPARAMETER("remove", &removeSam)
		END_LONG_PARAMETERS();

    int start_time = time(0);

	inputParameters.Add(new LongParameters(" Command Line Options: ",longParameterList));
	DosageConvertorVersion();
    String compStatus;
	inputParameters.Read(argc, &(argv[0]));
	if (help)
	{
		helpFile();
		return(-1);
	}

	inputParameters.Status();
    inputDose.gzip=!nobgzip;


    if(inputDose.checkValidity()==-1)
    {
        cout<< " Try --help or Contact author [sayantan@umich.edu] if you still need help ... "<<endl<<endl;
        cout<< " Program Exiting ...\n\n";
        compStatus="Command.Line.Error";
        PhoneHome::completionStatus(compStatus.c_str());
        return(-1);
    }

    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                             INPUT VCF DOSAGE FILE                             "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

	if (!inputDose.FastLoadHaplotypes())
	{
		cout << "\n Program Exiting ... \n\n";
		compStatus="Input.VCF.Dose.Load.Error";
        PhoneHome::completionStatus(compStatus.c_str());
        return(-1);

	}

    cout<<"\n ------------------------------------------------------------------------------"<<endl;
    cout<<"                                END OF PROGRAM                                 "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

	int time_tot = time(0) - start_time;

    cout << "\n Program Successfully Implemented... \n ";


	printf("\n Total Run completed in %d hours, %d mins, %d seconds.\n",
		time_tot / 3600, (time_tot % 3600) / 60, time_tot % 60);

    cout<<"\n Thank You for using DosageConvertor !!! "<<endl<<endl;


    compStatus="Success";
    PhoneHome::completionStatus(compStatus.c_str());

	return 0;

}




void DosageConvertorVersion()
{
	printf("\n\n -------------------------------------------------------------------------------- \n");
	printf("          DosageConvertor - Converting Dosage from Minimac3/4 to other Formats     \n");
	printf(" --------------------------------------------------------------------------------\n");
    printf(" (c) 2015 - Sayantan Das \n");
	cout<<"\n Version: " << VERSION<< ";\n Built: " << DATE << " by " << USER << std::endl;

  printf("\n URL = http://genome.sph.umich.edu/wiki/DosageConvertor\n");
  printf(" GIT = https://github.com/Santy-8128/DosageConvertor\n");
}

void helpFile()
{



    printf("\n DosageConvertor converts dosage files from minimac3/4 to other formats.\n");


	printf("\n Usage: ./DosageConvertor  --vcfDose      TestDataImputedVCF.dose.vcf.gz");
	printf("\n                           --info         TestDataImputedVCF.info (NOT mandatory)");
	printf("\n                           --prefix       OutputFilePrefix");
	printf("\n                           --type         plink OR mach   // depending on output format");
	printf("\n                           --format       1 or 2          // based on if you want to output");
	printf("\n                                                          // dosage (1) or genotype prob (2)");
    printf("\n                           --buffer       10000           // Number of Markers to import and ");
	printf("\n                                                          // print at a time (valid only for ");
	printf("\n                                                          // MaCH format)");
    printf("\n                           --idDelimiter  _               // Delimiter to Split VCF Sample ID into");
	printf("\n                                                          // FID and IID for PLINK format ");


	printf("\n                           --allDiploid                   // Handle to assume all samples diploid");
	printf("\n                           --TrimAlleles                  // Handle to trim length of alleles");
	printf("\n                           --SexFile      TestDataSexFile // First column sample names");
	printf("\n                                                          // Second columns M or F");

    printf("\n\n URL = http://genome.sph.umich.edu/wiki/DosageConvertor\n");
    printf(" GIT = https://github.com/Santy-8128/DosageConvertor\n");
  printf("\n Visit website for more details ...\n");

cout<<endl<<endl;
	return;
}
