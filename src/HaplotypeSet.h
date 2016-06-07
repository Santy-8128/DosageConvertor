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
        String Type;
		vector<string> individualName,familyName;
        int doseIndex,gpIndex,gtIndex;



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


        HaplotypeSet(String vcfDose, String info, String Outfile,bool Gzip,
                     String format,String type,String idDelimiter,int bufferSize,String myChromosome)
        {
            VcfDose=vcfDose;
            Info=info;
            outFile=Outfile;
            gzip=Gzip;
            Format=format;
            Type=type;
            IdDelimiter=idDelimiter;
            BufferSize=bufferSize;
            MyChromosome=myChromosome.c_str();
        }



        bool        ConvertDosageData           ();
        bool        LoadInfoFromVCFfile         ();
        bool        CheckValidChrom             (string chr);
        bool        FastLoadHaplotypes          ();
        string      DetectReferenceFileType     (String filename);

};









#endif // HAPLOTYPESET_H_INCLUDED
