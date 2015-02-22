#ifndef HAPLOTYPESET_H_INCLUDED
#define HAPLOTYPESET_H_INCLUDED

#include "StringBasics.h"
#include "VcfFileReader.h"
#include <fstream>
#include <sstream>


using namespace std;
class ReducedHaplotypeInfo
{
    public:
        int startIndex,endIndex;
        //vector<int> uniqueRep;  // in reduced state space (which is finally enumerated by this vector) ... indices of original haplotypes which are representatives
        vector<int> uniqueCardinality; // has number of representatives for each unique representative
        vector<int> uniqueIndexMap; // maps to corresponding item in the uniqueRep... required to map Constants and unfold fold probs


        vector<vector<char> > uniqueHaps;


        char returnHapAtPosition(int i,int position)
        {
            //assert((position-startIndex)>=0);
            return uniqueHaps[i][position-startIndex];
        }



    ReducedHaplotypeInfo()
    {

        startIndex=0;
        endIndex=0;
        //uniqueRep.clear();
        uniqueCardinality.clear();
        uniqueIndexMap.clear();
        uniqueHaps.clear();
    }

};


class variant
{
public:

    string name;
    int bp;
    string chr;
    char refAllele,altAllele;
    string refAlleleString,altAlleleString;
    string MajAlleleString,MinAlleleString;
    bool swapped;
    variant(){};
    variant(string &id,string &CHR,int &BP)
    {
        name=id;
        bp=BP;
        chr=CHR;
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

//        bool DS,GP;
        vector<int>        optEndPoints;
		vector<int>        ScaffoldIndex;
		vector<int>        UnScaffoldIndex;
		vector<ReducedHaplotypeInfo> ReducedStructureInfo;
		//vector<vector<char> >     haplotypes;
		vector<vector<char> >     haplotypesUnscaffolded;
		vector<vector<double> > alleleFreq;
		vector<vector<double> > Dosage;

		vector<vector<float> > dosage;
		vector<vector<float> > GP1;
		vector<vector<float> > GP2;
		vector<vector<char> > ImputedAlleles;
		int PrintStartIndex,PrintEndIndex;
        bool GT,DS,GP;
        int doseIndex,gpIndex,gtIndex;
        bool onlyCompress,filter;
        bool EstimateNcompress;
        bool PseudoAutosomal;
        bool AllMaleTarget;
        String removeSam;
        vector<double>      Recom,Error;
        vector<string> markerName;
		vector<string> individualName;
		vector<int> SampleNoHaplotypes;
		vector<char> refAlleleList,major, minor;
		vector<bool> missing, MarkerIndices;
		vector<variant> VariantList;
        string finChromosome;
		bool allowMissing, vcfType,m3vcfxType,machType;

        String VcfDose;
        String Info;
        String outFile;
        bool gzip;
        String Format;
        String Type;

        bool WriteMachFile();
//        numHaplotypes = 0;
//        numMarkers = 0;
//        ReducedStructureInfo.clear();
//        alleleFreq.clear();
//        //haplotypes.clear();
//        individualName.clear();
//        SampleNoHaplotypes.clear();
//
//
//        markerName.clear();
//        refAlleleList.clear();
//        major.clear();
//        minor.clear();
//        missing.clear();
//        MarkerIndices.clear();
//        allowMissing = true;
//        vcfType = false;
//        m3vcfxType=false;
//        machType=false;
//        PrintStartIndex=0;
//        PrintEndIndex=0;
//        Recom.clear();
//        Error.clear();


        HaplotypeSet(String vcfDose, String info, String Outfile,bool Gzip,String format,String type)
        {
            VcfDose=vcfDose;
            Info=info;
            outFile=Outfile;
            gzip=Gzip;
            Format=format;
            Type=type;

        }


        char    getScaffoldedHaplotype          (int sample,int marker);
        void    Create                          (vector<char> &tempHaplotype);
        void    calculateFreq                   ();
		bool    LoadInfoFile                     (String filename);
		bool    LoadMachHaplotypes              (String filename, String targetSnpfile, vector<string> &refSnpList);
		bool    LoadMachHaplotypes              (String filename, String targetSnpfile);
        char    convertAlleles                  (string markerId, string indivId, const char *alleles, string refAlleleString, string altAlleleString);
		bool    LoadHaplotypes                  (String filename, String snpNames, int maxIndiv, int maxMarker,bool optmx,String CNO,int START,int END);
		bool    FastLoadHaplotypes              ();
		bool    LoadTargetHaplotypes            (String filename, String targetSnpfile, vector<string> &refSnpList, HaplotypeSet &rHap);
		bool    LoadVcfTargetHaplotypes         (String filename, String snpNames, vector<string> &refSnpList,HaplotypeSet &rHap);
        void    convertToReducedStructure       (vector<int> &optStructure);
        void    writem3vcfFile                  (String &filename,bool &gzip);
        bool    readm3vcfFile                   (String m3vcfFile,String CHR,int START,int END,int WINDOW);
        void    reconstructHaplotype            (vector<char> &reHaplotypes,int &index);
        void    SaveDosageForVcfOutput          (int hapID,vector<double> dose,vector<char> impAlleles);
        void    SaveDosageForVcfOutputSampleWise(int SamID,string &SampleName, vector<double> &dose1,vector<double> &dose2,vector<char> &impAlleles1,vector<char> &impAlleles2);
        void    InitializeDosageForVcfOutput    (int NHaps,int NMarkers);
        void    InitializePartialDosageForVcfOutput    (int NHaps,int NMarkers, vector<bool> &Format);
        void    InitializePartialDosageForVcfOutputMaleSamples    (int NHaps,int NMarkers, vector<bool> &Format);
        void    PrintDosageForVcfOutputForID    (IFILE vcfdose, int MarkerIndex,bool majorIsReference,char refAllele);
        void    PrintPartialDosageForVcfOutputForID    (IFILE vcfdose, int MarkerIndex,bool majorIsReference,char refAllele);
        bool    BasicCheckForTargetHaplotypes   (String filename);
        string    DetectTargetFileType   (String filename);
        string    DetectReferenceFileType   (String filename);
        void    SaveDosageForVcfOutputSampleWiseChrX(int SamID,string &SampleName, vector<double> &dose1,
                                                    vector<char> &impAlleles1);

void PrintDosageForVcfOutputForIDMaleSamples(IFILE vcfdose, int MarkerIndex,bool majorIsReference,char refAllele);


};









#endif // HAPLOTYPESET_H_INCLUDED
