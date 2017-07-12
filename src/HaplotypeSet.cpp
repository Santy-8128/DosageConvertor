#include "HaplotypeSet.h"


bool HaplotypeSet::ConvertDosageData()
{

    String filename=VcfDose;

    int numSamplesRead = numSamples,markerCount=0,totmarkerCount=0;



    machPartial=NULL,plinkMain=NULL,plinkMap=NULL;
    int part=1;

    if(Type=="mach")
        machPartial = ifopen(outFile + ".output.part.0" + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
    else
    {
        plinkMain = ifopen(outFile + ".plink.fam", "wb",InputFile::UNCOMPRESSED);
        for (int i = 0; i < numSamplesRead; i++)
        {
             ifprintf(plinkMain,"%s\t%s\t0\t0\t%s\t-9\n",familyName[i].c_str(),individualName[i].c_str(), Sex[i].c_str());
        }
        ifclose(plinkMain);
        plinkMap = ifopen(outFile + ".plink.map", "wb",InputFile::UNCOMPRESSED);
        plinkMain = ifopen(outFile + ".plink.dosage" + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);

        ifprintf(plinkMain,"SNP\tA1\tA2");
        for (int i = 0; i < numSamplesRead; i++)
        {
             ifprintf(plinkMain,"\t%s\t%s",familyName[i].c_str(),individualName[i].c_str());
        }
        ifprintf(plinkMain,"\n");
    }


    if(Type=="mach")
    for (int i = 0; i < numSamplesRead; i++)
    {
         ifprintf(machPartial,"%s->%s\tDOSE\n",familyName[i].c_str(),individualName[i].c_str());
    }

    if(Type=="mach")
        ifclose(machPartial);

    InitializeVariables();

    int factor=10000;
    string line;

    std::cout << " Reading and Importing Data from VCF File ..."<<endl<<endl ;
	VcfFileReader inFile;
	VcfHeader header;
	VcfRecord record;
    if (!inFile.open(VcfDose, header))
	{
		cout << "\n Program could NOT open file : " << VcfDose << endl<<endl;
		return false;
	}

    while (totmarkerCount<numMarkers)
    {
        markerCount=0;
        if(Type=="plink")
        {
            if(totmarkerCount%factor==0)
            {
                printf("    Writing to PLINK output file ...");
                cout<<endl;
            }
        }
        while(markerCount<BufferSize && totmarkerCount<numMarkers)
        {
            inFile.readRecord(record);
            if(tag=="DS")
            {

                for (int i = 0; i<(numSamplesRead); i++)
                {
                    const string *s=(record.getGenotypeInfo().getString("DS",i));
                    if(s!= NULL)
                        dosage[markerCount][i]=(SampleNoHaplotypes[i]=="1"?2:1)*atof(s->c_str());
                    else
                        return PrintErr(" Tag DS does NOT exist for "+(string)record.getIDStr());

                    if(Format!="1")
                        return PrintErr(" Cannot generate genotype likelihoods from --tag DS \n Please use --tag GP or GT OR use --format 1 ");
                }
            }
            else if(tag=="GP")
            {
                for (int i = 0; i<(numSamplesRead); i++)
                {
                    const string *s=(record.getGenotypeInfo().getString("GP",i));
                    if(s!= NULL)
                    {
                        char *pch, *end_str;
                        pch = strtok_r ((char*)s->c_str(),",",&end_str);
                        GP1[markerCount][i]=atof(pch);
                        pch = strtok_r (NULL,",", &end_str);
                        GP2[markerCount][i]=(SampleNoHaplotypes[i]=="1"?0:atof(pch));
                        dosage[markerCount][i]=2*(1-GP1[markerCount][i]-GP2[markerCount][i]) + GP2[markerCount][i];
                    }
                    else
                        return PrintErr(" Tag GP does NOT exist for "+(string)record.getIDStr());
                }
            }
            else if(tag=="GT")
            {
                for (int i = 0; i<(numSamplesRead); i++)
                {
                    int a1,a2;
                    const string *s=(record.getGenotypeInfo().getString("GT",i));
                    if(s!= NULL)
                    {
                        char *pch, *end_str;
                        pch = strtok_r ((char*)s->c_str(),"|/",&end_str);
                        a1=atoi(pch);
                        if(SampleNoHaplotypes[i]=="2")
                        {
                            pch = strtok_r (NULL,"|/", &end_str);
                            if(pch==NULL) return PrintErr(" Incorrect Ploidy for Sample ID : " +individualName[i] + " at Marker ID : " + (string)record.getIDStr());
                            a2=atoi(pch);
                        }
                        else
                        {
                            a2=a1;
                        }

                        GP1[markerCount][i]=0;
                        GP2[markerCount][i]=0;
                        if((1-a1)*(1-a2) == 1)
                        {
                            GP1[markerCount][i]=1;
                        }
                        else if (a1+a2==1)
                        {
                            GP2[markerCount][i]=1;
                        }

                        dosage[markerCount][i]=2*(1-GP1[markerCount][i]-GP2[markerCount][i]) + GP2[markerCount][i];
//                        dosage[markerCount][i]=2*GP1[markerCount][i]+GP2[markerCount][i];
                    }
                    else
                        return PrintErr(" Tag GT does NOT exist for "+(string)record.getIDStr());
                }
            }

            if(Type=="plink")
                WritePlinkThisVariant(markerCount,totmarkerCount);

            markerCount++;
            totmarkerCount++;
        }
        if(Type=="plink")
        {
            if(totmarkerCount%factor==0)
            {
                printf("    Finished Writing %d markers to PLINK output file ...", totmarkerCount);
                cout<<endl;
            }
        }
        else if(Type=="mach")
        {
             WriteMaCHThisChunk(markerCount,totmarkerCount,part);
        }

    }
    if(Type=="plink")
    {
        printf("    Finished Writing ALL markers to PLINK output file ...");
        cout<<endl;
    }

    std::cout << "\n\n Number of Samples read from VCF Dosage File              : " << numSamplesRead ;
	std::cout << "\n Number of Markers read from VCF Dosage File              : " << totmarkerCount ;

    if(Type=="mach")
    {
        MergeMaCHOutputFile(part);
        if(!CreateInfoFile())
            return false;
    }
    else
    {
        PrintPlinkMessage();
    }

    cout<<endl;

    return true;

}


void HaplotypeSet::InitializeVariables()
{
    BufferSize=min(BufferSize,numMarkers);
    if(Type=="plink")
        BufferSize=1;

    dosage.clear();
    GP1.clear();
    GP2.clear();

    dosage.resize(BufferSize);
    for (int i = 0; i < BufferSize; i++)
    {
         dosage[i].resize(numSamples,-9.0);
    }

    GP1.resize(BufferSize);
    GP2.resize(BufferSize);
    for (int i = 0; i < BufferSize; i++)
    {
        GP1[i].resize(numSamples,-9.0);
        GP2[i].resize(numSamples,-9.0);
    }
}







void HaplotypeSet:: PrintPlinkMessage()
{

    cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
    cout<<"                                PLINK OUTPUT FILE                            "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl<<endl;


    printf("\n Dosage Information written to PLINK Dosage File  : %s ",(outFile + ".plink.dosage" + (gzip ? ".gz" : "")).c_str());
    printf("\n Map Information written to PLINK Map File        : %s ",(outFile + ".plink.map").c_str());
    printf("\n Family Information written to PLINK Fam File     : %s ",(outFile + ".plink.fam").c_str());

    ifclose(plinkMain);
    ifclose(plinkMap);

}



void HaplotypeSet:: WritePlinkThisVariant(int &markerCount,int &totmarkerCount)
{
    ifprintf(plinkMap,"%s\t%s\t0\t%d\n",VariantList[totmarkerCount].chr.c_str(),
            VariantList[totmarkerCount].name.c_str(),
            VariantList[totmarkerCount].bp);
    ifprintf(plinkMain,"%s\t%s\t%s",VariantList[totmarkerCount].name.c_str(),
            VariantList[totmarkerCount].refAlleleString.c_str(),
            VariantList[totmarkerCount].altAlleleString.c_str());

    for (int i = 0; i < numSamples; i++)
    {
        if(Format=="1")
            ifprintf(plinkMain,"\t%.3f",2-dosage[markerCount][i]);
        else if(Format=="2")
            ifprintf(plinkMain,"\t%.3f\t%.3f",GP1[markerCount][i],GP2[markerCount][i]);
        else if(Format=="3")
            ifprintf(plinkMain,"\t%.3f\t%.3f\t%.3f",GP1[markerCount][i],GP2[markerCount][i],1-GP1[markerCount][i]-GP2[markerCount][i]);
    }
    ifprintf(plinkMain,"\n");
}



void HaplotypeSet:: WriteMaCHThisChunk(int &markerCount,int &totmarkerCount, int &part)
{
    printf("    Finished Reading Chunk %d of %d markers ... \n", part, markerCount);
    machPartial = ifopen(outFile + ".output.part."+part + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
    for (int i = 0; i < numSamples; i++)
    {
        for (int j = 0; j < markerCount; j++)
        {
            if(Format=="1")
                ifprintf(machPartial,"\t%.3f",dosage[j][i]);
            else if(Format=="2")
                ifprintf(machPartial,"\t%.3f\t%.3f",GP1[j][i],GP2[j][i]);
        }
        ifprintf(machPartial,"\n");
    }
    ifclose(machPartial);
    printf("    Finished Writing Chunk %d of %d markers to MaCH Output File ... \n", part, markerCount);
    part++;

}


void HaplotypeSet:: MergeMaCHOutputFile(int &part)
{
    cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
    cout<<"                             MERGE MACH OUTPUT FILE                            "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl<<endl;
    IFILE machFinal=NULL;

    if(Format=="1")
    {
        printf("\n Merging Split Files to MaCH Dosage File : %s \n\n",(outFile + ".mach.dose" + (gzip ? ".gz" : "")).c_str());
        machFinal = ifopen(outFile + ".mach.dose" + (gzip ? ".gz" : ""),  "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
    }
    else if(Format=="2")
    {
        printf("\n Merging Split Files to MaCH GP File : %s \n\n",(outFile + ".mach.gprob" + (gzip ? ".gz" : "")).c_str());
        machFinal = ifopen(outFile + ".mach.gprob" + (gzip ? ".gz" : ""),  "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
    }

    vector<IFILE> machFinalList(part);
    string line;
    for(int i=0;i<part;i++)
    {
        string tempFileIndex(outFile);
        stringstream strs;
        strs<<(i);

        tempFileIndex+=(".output.part." +
                        (string)(strs.str())+(gzip ? ".gz" : ""));
        machFinalList[i] = ifopen(tempFileIndex.c_str(), "r");
    }

    for(int j=0;j<numSamples;j++)
    {
        if(j%500==0)
        {
              printf("     Merging Sample %d of %d to MaCH File ...", j + 1, numSamples);
              cout<<endl;
        }
        for(int i=0;i<part;i++)
        {
            line.clear();
            machFinalList[i]->readLine(line);
            ifprintf(machFinal,"%s",line.c_str());
        }
        ifprintf(machFinal,"\n");
    }
    ifclose(machFinal);
    printf("     Merging Finished ... \n");

    for(int i=0;i<part;i++)
    {
        ifclose(machFinalList[i]);
        string tempFileIndex(outFile);
        stringstream strs;
        strs<<(i);

        tempFileIndex+=(".output.part." +
                        (string)(strs.str())+(gzip ? ".gz" : ""));
        remove(tempFileIndex.c_str());

    }

}


bool HaplotypeSet:: CreateInfoFile()
{

    printf("\n Creating Info file : %s \n\n",(outFile + ".mach.info").c_str());
    IFILE machInfo = ifopen(outFile + ".mach.info",  "wb", InputFile::UNCOMPRESSED);
    string line;

    if(IsInfo)
    {
        printf("\n Copying from input Info file : %s \n\n",Info.c_str());
        line.clear();
        IFILE ifs = ifopen(Info, "r");
        int Indic=ifs->readLine(line);
        int Infocount=-1;

        while(Indic==-1 ||  Indic==0)
        {
            if(Indic==-1 && line.length()==0)
            {
                break;
            }

            ifprintf(machInfo,"%s\n",line.c_str());
            Infocount++;
            line.clear();
            Indic=ifs->readLine(line);
        }
        if(Infocount!=numMarkers)
        {
            cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
            cout<<"                                    ERROR !!!                                 "<<endl;
            cout<<" ------------------------------------------------------------------------------"<<endl;
          std::cout << "\n\n             "<<Infocount<<" markers found in Input info file <"<< Info <<"> " ;
            std::cout << "\n             However, "<<numMarkers <<" markers were read from dosage VCF file <"<< VcfDose<< "> " ;
            std::cout << "\n             Please use info file from same imputation run ! " << endl;
            remove((outFile + ".mach." + (Format=="1" ? "dose" : "gprob")  + (gzip ? ".gz" : "")).c_str());
            return false;
        }
    }
    else
    {
        cout<<" ------------------------------------------------------------------------------"<<endl;
        cout<<"                                   WARNING !!!                                 "<<endl;
        cout<<" ------------------------------------------------------------------------------"<<endl<<endl;

        cout<< " No Info file provided by user. \n";
        cout<< " Final MaCH info file <"<<outFile + ".mach.info" << "> will have some missing columns.\n";
        cout<< " Please use INFO file from Minimac3 for complete MaCH info file ... \n\n";


        ifprintf(machInfo, "SNP\tREF(0)\tALT(1)\tALT_Frq\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose0\tDose1\n");

        for (int i =0; i < numMarkers; i++)
        {
            variant &tempVariant=VariantList[i];

        ifprintf(machInfo, "%s\t%s\t%s\t%s\t%s\t-\t%s\t%s\t-\t-\t%s\t-\t-\n",
            tempVariant.name.c_str(),
            tempVariant.refAlleleString.c_str(),
            tempVariant.altAlleleString.c_str(),
            tempVariant.Af.c_str(),
            tempVariant.Maf.c_str(),
            tempVariant.Rsq.c_str(),
            tempVariant.Tag.c_str(),
            tempVariant.Ersq.c_str());
        }

    }

    ifclose(machInfo);
    cout<<" \n [NOTE: In output MaCH files, dosage values are twice the probability of";
    printf("\n        ALTERNATE ALLELE (in the info file) and NOT twice the probability");
    cout<< "\n        of MAJOR ALLELE (as is in usual Minimac output)"<<endl;
    cout<<   "        Similarly, Genotype probabilities are given for REF|REF followed"<<endl;
    cout<<   "        REF|ALT (allele given in info file)."<<endl<<endl;


    if(Format=="1")
    {
        printf("\n Dosage Information written to MaCH Dosage File  : %s ",(outFile + ".mach.dose" + (gzip ? ".gz" : "")).c_str());
        printf("\n Information written to MaCH Info File           : %s ",(outFile + ".mach.info").c_str());

    }
    else if(Format=="2")
    {
        printf("\n Genotype Probability written to MaCH Dosage File  : %s ",(outFile + ".mach.gprob" + (gzip ? ".gz" : "")).c_str());
        printf("\n Information written to MaCH Info File             : %s ",(outFile + ".mach.info").c_str());
    }

    return true;
}



bool HaplotypeSet::LoadInfoFromVCFfile()
{

	VcfFileReader inFile;
	VcfHeader header;
	VcfRecord record;
    int TypedOnlyCount=0,GenotypedCount=0,ImputeCount=0;

	if(!GetSampleInformation())
        return false;

    if (!inFile.open(VcfDose, header))
	{
		cout << "\n Program could NOT open file : " << VcfDose << endl<<endl;
		return false;
	}

    VariantList.clear();
    string name,refAlleleString,altAlleleString;
    int bp;
    string chr;
    std::cout << "\n Calculating number Markers in VCF File ..."<<endl;

    while (inFile.readRecord(record))
	{
		chr=record.getChromStr();
        bp=record.get1BasedPosition();
        name=record.getIDStr();

        if(!CheckValidChrom(chr))
        {
            return PrintErr(" Reference VCF File contains chromosome : "+chr+"\n VCF File can only contain chromosomes 1-22 ! ");
        }


        refAlleleString = record.getRefStr();
		altAlleleString = record.getAltStr();


		if(TrimAlleles)
        {
            name = name.substr(0, TrimLength);
//            refAlleleString = refAlleleString.substr(0, 100);
//            altAlleleString = altAlleleString.substr(0, 100);
        }


        variant tempVariant(name,chr,bp);
        tempVariant.assignRefAlt(refAlleleString,altAlleleString);

        VcfRecordFilter& ThisFilter=record.getFilter();
        for(int i=0;i<ThisFilter.getNumFilters();i++)
        {
            string Filter=ThisFilter.getString(i);

            if(Filter=="GENOTYPED_ONLY")
            {
                tempVariant.genotyped_only=true;
                tempVariant.Tag="Typed_Only";
            }
            else if(Filter=="GENOTYPED")
            {
                tempVariant.genotyped=true;
                tempVariant.Tag="Genotyped";
            }

        }

        VcfRecordInfo& ThisInfo=record.getInfo();

        if(ThisInfo.getString("TYPED")!=NULL)
        {
            tempVariant.genotyped=true;
            tempVariant.Tag="Genotyped";
        }
        if(ThisInfo.getString("TYPED_ONLY")!=NULL)
        {
            tempVariant.genotyped_only=true;
            tempVariant.Tag="Typed_Only";
        }

        if(ThisInfo.getString("AF")!=NULL)
            tempVariant.Af=*ThisInfo.getString("AF");
        if(ThisInfo.getString("MAF")!=NULL)
            tempVariant.Maf=*ThisInfo.getString("MAF");

        if(ThisInfo.getString("R2")!=NULL)
            tempVariant.Rsq=*ThisInfo.getString("R2");
        else
        {
            if(!tempVariant.genotyped_only)
            {
                return PrintErr(" Variant ["+name+"] does NOT have R2 in the INFO column, although it is NOT GENOTYPED_ONLY ! ");
            }
        }

        if(ThisInfo.getString("ER2")!=NULL)
        {
            if(!tempVariant.genotyped)
            {
                return PrintErr(" Variant ["+name+"] was GENOTYPED, but does NOT have ER2 in the INFO column ! ");
            }
            tempVariant.Ersq=*ThisInfo.getString("ER2");
        }


        if(tempVariant.Tag=="Genotyped")
            GenotypedCount++;
        else if(tempVariant.Tag=="Typed_Only")
            TypedOnlyCount++;
        else
            ImputeCount++;


        VariantList.push_back(tempVariant);
	}

    numMarkers=VariantList.size();
	inFile.close();



  std::cout << "\n\n Number of Samples read from VCF Dosage File              : " << numSamples ;
    std::cout << "\n Number of Markers read from VCF Dosage File              : " << numMarkers <<endl<<endl;


  std::cout << "\n\n Number of Imputed Variants                               : " << ImputeCount;
	std::cout << "\n Number of Genotyped and Imputed Variants                 : " << GenotypedCount;
	std::cout << "\n Number of Genotyped Only Variants                        : " << TypedOnlyCount<<endl<<endl;


    return true;


}



bool HaplotypeSet::GetSampleInformation()
{
	VcfFileReader inFile;
	VcfHeader header;
	VcfRecord record;

    individualName.clear();
    familyName.clear();
    String filename=VcfDose;

    if (!inFile.open(filename, header))
	{
		cout << "\n Program could NOT open file : " << filename << endl<<endl;
		return false;
	}

    std::cout << "\n Calculating number of Samples in VCF File ..." ;

	numSamples = header.getNumSamples();

    if(numSamples==0)
    {
        cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
        cout<<"                                    ERROR !!!                                  "<<endl;
        cout<<" ------------------------------------------------------------------------------"<<endl;
        return PrintErr("  NO samples found in Input VCF Dosage File ! ");
    }

    vector<string> OrigNames;

    for (int i = 0; i < numSamples; i++)
	{
		string tempName(header.getSampleName(i));
		OrigNames.push_back(tempName);
		if(IdDelimiter=="")
        {
            individualName.push_back(tempName);
            familyName.push_back(tempName);
        }
		else
        {

            size_t pos = 0;
            std::string delimiter(IdDelimiter) ;
            std::string token;
            int Count=0;
            string temptempName=tempName;

            while ((pos = tempName.find(delimiter)) != std::string::npos)
            {
                token = tempName.substr(0, pos);

                if(Count==0)
                    familyName.push_back(token);

                tempName.erase(0, pos + delimiter.length());
                if(Count==1)
                {
                    string temp=" Program could NOT parse the following sample name with ID delimiter ["
                                + (string)IdDelimiter + "] : " + temptempName + "\n More than TWO tokens found : ["
                                + familyName.back()  + "] [" + token + "] [" + tempName + "] \n ";
                    return PrintErr(temp);
                }

                Count++;


            }
            if(Count==0)
            {
                string temp=" Program could NOT parse the following sample name with ID delimiter ["
                                + (string)IdDelimiter + "] : " + temptempName
                                + "\n Delimiter NOT FOUND in the sample name ! \n ";
                return PrintErr(temp);
            }

            individualName.push_back(tempName);
        }
    }


    inFile.readRecord(record);

    finChrom=record.getChromStr();
    if(!CheckValidChrom((string)finChrom))
    {
        return PrintErr(" Reference VCF File contains chromosome : "+(string)finChrom+"\n VCF File can only contain chromosomes 1-22, X ! ");
    }

    SampleNoHaplotypes.resize(numSamples,"2");
    Sex.resize(numSamples,"0");

    if(finChrom=="X")
    {
        cout<<"\n\n Chromosome X detected !"<<endl;
        if(samePloidy==true)
        {
            cout<<" Assuming all samples are diploid (--allDiploid) ! "  <<endl;
            SampleNoHaplotypes.resize(numSamples,"2");
        }

        if(SexFile!="")
        {
            vector<string> tempNames;
            vector<string> tempPloidy;

            cout<<" Reading sex information from : "<< SexFile  <<endl;

            IFILE sexInfo = ifopen(SexFile, "r");
            string line;
            sexInfo->readLine(line);
            while(line!="")
            {
                char *pch, *end_str;
                pch = strtok_r ((char*)line.c_str()," \t",&end_str);
                if(pch==NULL) return PrintErr(" Invalid input for --sexFile !");
                tempNames.push_back(pch);

                pch = strtok_r (NULL," \t", &end_str);
                if(pch==NULL) return PrintErr(" Invalid input for --sexFile !");
                if(strcmp(pch,"M")==0)
                    tempPloidy.push_back("1");
                else if(strcmp(pch,"F")==0)
                    tempPloidy.push_back("2");
                else
                    tempPloidy.push_back("3");
                line.clear();
                sexInfo->readLine(line);
            }

            for(int i=0;i<numSamples;i++)
            {
                bool found=false;

                for(int j=0;j<(int)tempNames.size();j++)
                {
                    if(tempNames[j]==OrigNames[i])
                    {
                        found=true;
                        if(strcmp(tempPloidy[j].c_str(),"3")==0)
                            return PrintErr(" Incorrect Ploidy for Sample ID : "+ OrigNames[i] +" \n");
                        SampleNoHaplotypes[i] = samePloidy ? SampleNoHaplotypes[i] : tempPloidy[j];
                        Sex[i] = tempPloidy[j];


                        break;
                    }
                }
                if(found==false)
                {
                    return PrintErr(" Sample ID : "+ OrigNames[i] +" not found in SexFile [" + (string)SexFile + "] \n");
                }
            }

            ifclose(sexInfo);
        }


        if(SexFile=="" && !samePloidy)
        {
            cout<<" Determining ploidy from GT tag in VCF file ! "  <<endl;
            bool AllTwoPloidy=true;

            for (int i = 0; i<(numSamples); i++)
            {
                if(record.getNumGTs(i)==0)
                {
                    return PrintErr(" Chromosome X found and GT tags are NOT available ! \n Cannot determine ploidy for the samples ! \n Either turn ON handle --allDiploid, or else provide a --sexFile \n The 'sexFile' should have the sample names in the first column \n and 'M' or 'F' in the second column\n");
                }
                else
                {
                    stringstream ss;
                    ss<<record.getNumGTs(i);
                    SampleNoHaplotypes[i]=ss.str();
                    Sex[i]=SampleNoHaplotypes[i];
//                    tempHapCount+=SampleNoHaplotypes[i];
                }
//                cout<<" WHY "<<i<<"\t"<<Sex[i]<<"\t"<<Sex[i]<<endl;

                if(Sex[i]=="1")
                {
                    AllTwoPloidy=false;
                }

            }

            if(AllTwoPloidy)
            {
                Sex.clear();
                Sex.resize(numSamples,"0");
            }

        }

        int males = 0, females = 0;
        if(Sex[0]!="0")
        {
            for(int i=0;i<numSamples;i++)
            {
                if(strcmp(Sex[i].c_str(),"1")==0)
                    males++;
                else
                    females++;
            }
            cout<<"\n Number of males   = "<<males<<endl;
            cout<<" Number of females = "<<females<<endl;
        }
        else
        {
            cout<<" Could NOT determine sex since all samples have same ploidy 2 !"<<endl;
        }
    }



    return true;


}





bool HaplotypeSet::CheckValidChrom(string chr)
{
    bool result=false;

    if(MyChromosome!="" && chr==MyChromosome)
        return true;

    string temp[]={"1","2","3","4","5","6","7","8","9","10","11"
                    ,"12","13","14","15","16","17","18","19","20","21","22","X"};
    std::vector<string> ValidChromList (temp, temp + sizeof(temp) / sizeof(string) );

    for(int counter=0;counter<(int)ValidChromList.size();counter++)
        if(chr==ValidChromList[counter])
            result=true;

    return result;

}





bool HaplotypeSet::FastLoadHaplotypes()
{
    String filename=VcfDose;
    if(Info=="")
        IsInfo=false;
    else
        IsInfo=true;



    cout<<"\n Detecting Dosage File Type ... "<<endl;

    string FileType=DetectReferenceFileType(filename);

    if(FileType.compare("NA")==0)
    {
        cout<<"\n Following File File Not Available : "<<filename<<endl;
        return false;
    }
    else if(FileType.compare("Invalid")==0)
    {

        cout << "\n Dosage File provided by --vcfDose must be a VCF file ! \n";
        cout << " Please check the following file : "<<filename<<endl;
        return false;
    }

    cout<<"\n Format = VCF (Variant Call Format) "<<endl;



    if(!LoadInfoFromVCFfile())
        return false;



    std::cout << " Loading Dosage Data from VCF File       : " << filename << endl<<endl;

    return ConvertDosageData();

}



string HaplotypeSet::DetectReferenceFileType(String filename)
{
    IFILE fileStream = ifopen(filename, "r");
    string line;
    if(fileStream)
    {
        fileStream->readLine(line);
        if(line.length()<1)
            return "Invalid";
         string tempString;
        tempString=(line.substr(0,17));

        char temp[tempString.length() + 1];
        std::strcpy(temp,tempString.c_str());
        for (char *iter = temp; *iter != '\0'; ++iter)
        {
           *iter = std::tolower(*iter);
        }

        if(((string)temp).compare("##fileformat=vcfv")==0)
        {
            return "vcf";
        }
        else
            return "Invalid";

    }
    else
    {
        return "NA";
    }
    ifclose(fileStream);

    return "NA";
}


bool HaplotypeSet::PrintErr(string text)
{

cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
    cout<<"                                    ERROR !!!                                  "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl<<endl;
    cout << text <<endl;
    cout << " Please verify the input file(s) ! "<<endl;
    cout << " Type --help for usage details \n Contact author [sayantan@umich.edu] if you still need help ... "<<endl;
    cout << " Program Aborting ... "<<endl;
    return false;
}

