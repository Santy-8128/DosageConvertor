#ifndef HAPLOTYPESET_H_INCLUDED
#define HAPLOTYPESET_H_INCLUDED

#include "StringBasics.h"
#include "VcfFileReader.h"
#include <fstream>
#include <sstream>


using namespace std;

class variant
{
public:

    string name;
    int bp;
    string chr;
    char refAllele,altAllele;
    string refAlleleString,altAlleleString;
    string MajAlleleString,MinAlleleString;

    string Af, Maf, Rsq, Ersq;
    bool genotyped,genotyped_only;
    string Tag;



    bool swapped;
    variant(){};
    variant(string &id,string &CHR,int &BP)
    {
        name=id;
        bp=BP;
        chr=CHR;
        genotyped=false;
        genotyped_only=false;
        Af="-";
        Rsq="-";
        Ersq="-";
        Tag="Imputed";

    };
    void assignValues(string &id,string &CHR,int &BP)
    {
        name=id;
        bp=BP;
        chr=CHR;
    };
    void assignSwap(bool Swap)
    {
        swapped=Swap;
    };
     void assignRefAlt(string &refe,string &alt)
    {
        refAlleleString=refe;
        altAlleleString=alt;
    };
    void assignMajMin(string &refe,string &alt)
    {
        MajAlleleString=refe;
        MinAlleleString=alt;
    };


};


class HaplotypeSet
{

	public:
		int         numHaplotypes,numSamples;
		int         numMarkers;
		bool        IsInfo;
		string MyChromosome;

        vector<variant> VariantList;
		vector<vector<float> > dosage;
		vector<vector<float> > GP1;
		vector<vector<float> > GP2;
		bool GT,DS,GP;
		int BufferSize;
        String VcfDose;
        String Info,IdDelimiter;
        String outFile;
        bool gzip;
        String Format;
        String tag;
        String Type;
        String SexFile;
        bool samePloidy;
		vector<string> individualName,familyName;
		vector<string> Sex;
		vector<string> SampleNoHaplotypes;
		String finChrom;
        int doseIndex,gpIndex,gtIndex;
        int TrimLength;
        IFILE machPartial,plinkMain,plinkMap;

        bool TrimAlleles;

//
//
//		vector<vector<char> > ImputedAlleles;
//		int PrintStartIndex,PrintEndIndex;
//        bool onlyCompress,filter;
//        bool EstimateNcompress;
//        bool PseudoAutosomal;
//        bool AllMaleTarget;
//        String removeSam;
//        vector<double>      Recom,Error;
//        vector<string> markerName;
//		vector<int> SampleNoHaplotypes;
//		vector<char> refAlleleList,major, minor;
//		vector<bool> missing, MarkerIndices;
//		string finChromosome;
//		bool allowMissing, vcfType,m3vcfxType,machType;


        void        InitializeVariables         ();
        bool        GetSampleInformation        ();
        void        PrintPlinkMessage           ();
        bool        CreateInfoFile              ();
        void        MergeMaCHOutputFile         (int &part);
        void        WriteMaCHThisChunk          (int &markerCount,int &totmarkerCount, int &part);
        void        WritePlinkThisVariant       (int &markerCount,int &totmarkerCount);
        bool        PrintErr                    (string text);
        bool        ConvertDosageData           ();
        bool        LoadInfoFromVCFfile         ();
        bool        CheckValidChrom             (string chr);
        bool        FastLoadHaplotypes          ();
        string      DetectReferenceFileType     (String filename);



        HaplotypeSet()
        {
            VcfDose = "";
            Info = "";
            tag = "DS";
            outFile = "Converted.Dosage";
            Format = "1";
            Type = "plink";
            IdDelimiter = "";
            MyChromosome="";
            SexFile="";
            samePloidy=false;
            TrimAlleles=false;
            BufferSize =10000;
            gzip = true;
            TrimLength=100;
        }

        HaplotypeSet(String vcfDose, String info, String Outfile,bool Gzip,
                     String format,String type,String idDelimiter,int bufferSize,String myChromosome, String Tag)
        {
            VcfDose=vcfDose;
            Info=info;
            outFile=Outfile;
            gzip=Gzip;
            Format=format;
            Type=type;
            IdDelimiter=idDelimiter;
            BufferSize=bufferSize;
            tag=Tag;
            MyChromosome=myChromosome.c_str();
        }

        int checkValidity()
        {
            if (VcfDose == "")
            {
                cout<< " ERROR !!! Missing --vcfDose, a required parameter.\n\n";
                return(-1);
            }

            string tempString(Type);

            char temptype[tempString.length() + 1];
            std::strcpy(temptype,tempString.c_str());
            for (char *iter = temptype; *iter != '\0'; ++iter)
            {
               *iter = std::tolower(*iter);
            }

            if((string)temptype!="mach" && (string)temptype!="plink")
            {
                cout<< " ERROR !!! Invalid input for --type parameter : "<<Type<<endl;
                cout<< " Parameter must be equal to \"mach\" or \"plink\" ..."<<endl<<endl;
                return(-1);
            }

            if(tag!="DS" && tag!="GP" && tag!="GT")
            {
                cout<< " ERROR !!! Invalid input for --tag parameter : "<<Type<<endl;
                cout<< " Parameter must be equal to \"DS\" or \"GP\" or \"GT\" ..."<<endl<<endl;
                return(-1);
            }

            if (Info == "")
            {
                if((string)temptype=="mach")
                {
                    cout<< " WARNING !!! Missing --info parameter.\n";
                    cout<< " Final MaCH info file will have some missing columns.\n";
                    cout<< " Please use INFO file from Minimac3 for complete MaCH info file ... \n\n";
                }
            }
            if(tag=="DS")
            {
                if(Format!="1")
                {
                    cout<< " ERROR !!! Cannot generate genotype likelihoods from --tag \"DS\" "<<endl;
                    cout<< " Please use --tag \"GP\" or \"GT\" OR use --format \"1\" "<<Format<<endl<<endl;
                    return(-1);
                }
            }

            if(Format!="1" && Format!="2" && Format!="3")
            {
                cout<< " ERROR !!! Invalid input for --format parameter : "<<Format<<endl;
                cout<< " Parameter must be equal to \"1\" or \"2\" or \"3\" ..."<<endl<<endl;
                return(-1);
            }

             if((string)temptype=="mach" && Format=="3")
            {
                cout<< " ERROR !!! Invalid input for --format parameter : "<<Format<<endl;
                cout<< " If --type is \"mach\" then --format must be equal to \"1\" or \"2\" ..."<<endl<<endl;
                return(-1);
            }

            if(BufferSize<=0)
            {
                cout<< " ERROR !!! Invalid input for parameter --buffer : "<<BufferSize<<endl;
                cout<< " --buffer  can only take Positive Integers ...\n";
                return(-1);
            }

            if(TrimLength<20 || TrimLength>16000)
            {
                cout<< " ERROR !!! Invalid input for parameter --trimLength : "<<TrimLength<<endl;
                cout<< " --trimLength  can only take values between 20 and 16,000 ...\n";
                return(-1);
            }
            else
            {
                if((string)temptype=="plink")
                {
                    cout<< " WARNING !!! Parameter --buffer will be ignored when be used with --type \"plink\" !"<< endl<<endl;
                }
            }

            Type=(String)temptype;
            return 1;
        }





};









#endif // HAPLOTYPESET_H_INCLUDED
