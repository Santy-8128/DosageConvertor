

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

	String vcfDose = "";
	String info = "", snps = "",removeSam="";
	String outfile = "Converted.Dosage";
	String format = "DS";
	String type = "mach", errFile = "",chr="",golden="";
	String idDelimiter = "";
//	int max_indiv = 0, max_marker = 0;
    vector<bool> formatVector(3,false);

    int buffer=10000;

	bool  gzip = true, nobgzip = false;
	bool   help = false, params = false;

	ParameterList inputParameters;
	PhoneHome::allThinning = 50;

	BEGIN_LONG_PARAMETERS(longParameterList)
		LONG_PARAMETER_GROUP("Input Dosage")
		LONG_STRINGPARAMETER("vcfDose", &vcfDose)
		LONG_STRINGPARAMETER("info", &info)
		//LONG_PARAMETER("rsid", &rsid)
		LONG_PARAMETER_GROUP("Output Parameters")
		LONG_STRINGPARAMETER("prefix", &outfile)
		LONG_STRINGPARAMETER("format", &format)
		LONG_STRINGPARAMETER("type", &type)
		LONG_PARAMETER("nobgzip", &nobgzip)
		LONG_PARAMETER_GROUP("Other Parameters")
		LONG_INTPARAMETER("buffer", &buffer)
		LONG_STRINGPARAMETER("idDelimiter", &idDelimiter)
		LONG_PARAMETER("help", &help)
		LONG_PARAMETER("params", &params)
		LONG_PHONEHOME(VERSION)
		BEGIN_LEGACY_PARAMETERS()
//		LONG_PARAMETER("onlyRefMarkers", &onlyRefMarkers)
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
    if(nobgzip)
        gzip=false;

    string tempString(type);

    char temptype[tempString.length() + 1];
    std::strcpy(temptype,tempString.c_str());
    for (char *iter = temptype; *iter != '\0'; ++iter)
    {
       *iter = std::tolower(*iter);
    }

	if (vcfDose == "")
	{
		cout<< " Missing \"--vcfDose\", a required parameter.\n\n";
		cout<< " Try \"--help\" for usage ...\n\n";
		cout<< " Program Exiting ...\n\n";
		compStatus="Command.Line.Error";
		PhoneHome::completionStatus(compStatus.c_str());
		return(-1);
	}

	if (info == "")
	{
		cout<< " Missing \"--info\", a required parameter.\n\n";
		cout<< " Try \"--help\" for usage ...\n\n";
		cout<< " Program Exiting ...\n\n";
		compStatus="Command.Line.Error";
		PhoneHome::completionStatus(compStatus.c_str());
		return(-1);
	}

	if(format!="GP" && format!="DS")
    {
		cout<< " Invalid input for \"--format\" parameter : "<<format<<endl;
		cout<< " Parameter must be equal to \"GP\" or \"DS\" ..."<<endl<<endl;
		cout<< " Try \"--help\" for usage ...\n\n";
		cout<< " Program Exiting ...\n\n";
		compStatus="Command.Line.Error";
		PhoneHome::completionStatus(compStatus.c_str());
		return(-1);
	}

	if((string)temptype!="mach" && (string)temptype!="plink")
    {
		cout<< " Invalid input for \"--type\" parameter : "<<type<<endl;
		cout<< " Parameter must be equal to \"mach\" or \"plink\" ..."<<endl<<endl;
		cout<< " Try \"--help\" for usage ...\n\n";
		cout<< " Program Exiting ...\n\n";
		compStatus="Command.Line.Error";
		PhoneHome::completionStatus(compStatus.c_str());
		return(-1);
	}

    if(buffer<=0)
    {
        cout<< " Invalid input for parameter \"--buffer\" : "<<buffer<<endl;
        cout<< " \"--buffer\"  can only take Positive Integers ...\n";
        cout<< " Try \"--help\" for usage ...\n\n";
        cout<< " Program Exiting ...\n\n";
        compStatus="Command.Line.Error";
        PhoneHome::completionStatus(compStatus.c_str());
        return(-1);
    }
    else
    {
        if((string)temptype=="plink")
        {
            cout<< " WARNING: Parameter \"--buffer\" will be ignored when be used with \"--type\" = PLINK !!!"<< endl<<endl;
        }



    }




    type=(String)temptype;


    HaplotypeSet inputDose(vcfDose, info, outfile,gzip,format,type,idDelimiter,buffer);

    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                             INPUT VCF DOSAGE FILE                             "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

//
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
	printf("          DosageConvertor - Converting Dosage from Minimac3 to other Formats     \n");
	printf(" --------------------------------------------------------------------------------\n");
    printf(" (c) 2015 - Sayantan Das \n");
//	printf(" Version	: Undocumented Release\n");
//	printf(" Built		: sayantan\n\n");
	cout<<"\n Version: " << VERSION<< ";\n Built: " << DATE << " by " << USER << std::endl;
}

void helpFile()
{

  printf("\n URL = http://genome.sph.umich.edu/wiki/DosageConvertor\n");
  printf(" GIT = https://github.com/Santy-8128/DosageConvertor\n");


    printf("\n DosageConvertor converts dosage files from minimac3 to other formats.\n");


	printf("\n Usage: ./DosageConvertor  --vcfDose      TestDataImputedVCF.dose.vcf.gz");
	printf("\n                           --info         TestDataImputedVCF.info");
	printf("\n                           --prefix       OutputFilePrefix");
	printf("\n                           --type         plink OR mach   // depending on output format");
	printf("\n                           --format       DS or GP        // based on if you want to output");
	printf("\n                                                          // dosage (DS) or genotype prob (GP)");
    printf("\n                           --buffer       10000           // Number of Markers to import and ");
	printf("\n                                                          // print at a time (valid only for ");
	printf("\n                                                          // MaCH format)");
    printf("\n                           --idDelimiter  _               // Delimiter to Split VCF Sample ID into");
	printf("\n                                                          // FID and IID for PLINK format ");

    printf("\n\n URL = http://genome.sph.umich.edu/wiki/DosageConvertor\n");
    printf(" GIT = https://github.com/Santy-8128/DosageConvertor\n");
  printf("\n Visit website for more details ...\n");

cout<<endl<<endl;
	return;
}


