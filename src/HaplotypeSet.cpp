#include "HaplotypeSet.h"


bool HaplotypeSet::ConvertDosageData()
{

    String filename=VcfDose;

    int numSamplesRead = numSamples,markerCount=0,totmarkerCount=0,NumGP=0,NumDS=0,NumGT=0;

    char * pch_split,* pch_split3;
    char * pch;
    char *end_str1,*end_str3;


    IFILE machPartial=NULL,plinkMain=NULL,plinkMap=NULL;
    int part=1;

    if(Type=="mach")
        machPartial = ifopen(outFile + ".output.part.0" + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
    else
    {
        plinkMain = ifopen(outFile + ".plink.fam", "wb",InputFile::UNCOMPRESSED);
        for (int i = 0; i < numSamplesRead; i++)
        {
             ifprintf(plinkMain,"%s\t%s\t0\t0\t0\t-9\n",familyName[i].c_str(),individualName[i].c_str());
        }
        ifclose(plinkMain);
        plinkMap = ifopen(outFile + ".plink.map", "wb",InputFile::UNCOMPRESSED);
        plinkMain = ifopen(outFile + ".plink."+(Format=="DS" ? "dose" : "gprob") + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);

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

//    VariantList.clear();
    string name,refAlleleString,altAlleleString;
    int bp;
    string chr;

    int bufferSize=BufferSize;
    if(Type=="plink")
        bufferSize=1;

    dosage.clear();
    GP1.clear();
    GP2.clear();

    if(Format=="DS")
    {
        dosage.resize(bufferSize);
        for (int i = 0; i < bufferSize; i++)
        {
             dosage[i].resize(numSamplesRead,-9.0);
        }
    }


    if(Format=="GP")
    {
        GP1.resize(bufferSize);
        GP2.resize(bufferSize);
        for (int i = 0; i < bufferSize; i++)
        {
            GP1[i].resize(numSamplesRead,-9.0);
            GP2[i].resize(numSamplesRead,-9.0);
        }
    }
    int factor=10000;

    IFILE DoseRead = ifopen(filename, "r");
    string line;

    std::cout << " Reading and Importing Data from VCF File ..."<<endl<<endl ;

    if(DoseRead)
    {
        bool Header=true;

        while(Header)
        {
            line.clear();
            DoseRead->readLine(line);
            if(line.substr(0,1).compare("#")==0)
                Header=true;
            else
                Header=false;
        }

        while(line.compare("")!=0)
        {
            markerCount=0;

            while(markerCount<bufferSize && line.compare("")!=0)
            {

                char *end_str9;
                pch = strtok_r ((char*)line.c_str(),"\t", &end_str1);
                if(pch==NULL)
                {
                    cout<<"\n Input VCF file has invalid format at Marker No : "<<totmarkerCount+1<<endl;
                    cout<<" Please verify the following file : "<<filename<<endl;
                    return false;
                }
                else
                    chr=pch;
                pch = strtok_r (NULL, "\t", &end_str1);

                if(pch==NULL)
                {
                    cout<<"\n Input VCF file has invalid format at Marker No : "<<totmarkerCount+1<<endl;
                    cout<<" Please verify the following file : "<<filename<<endl;
                    return false;
                }
                else
                    bp=atoi(pch);

                pch = strtok_r (NULL, "\t", &end_str1);
                if(pch==NULL)
                {
                    cout<<"\n Input VCF file has invalid format at Marker No : "<<totmarkerCount+1<<endl;
                    cout<<" Please verify the following file : "<<filename<<endl;
                    return false;
                }
                else
                    name=pch;


                variant tempVariant(name,chr,bp);

                pch = strtok_r (NULL, "\t", &end_str1);
                pch = strtok_r (NULL, "\t", &end_str1);

                VariantList[totmarkerCount].assignValues(name,chr,bp);

                pch = strtok_r (NULL, "\t", &end_str1);
                pch = strtok_r (NULL, "\t", &end_str1);
                pch = strtok_r (NULL, "\t", &end_str1);
                pch = strtok_r (NULL, "\t", &end_str1);

                pch_split = strtok_r (pch,":", &end_str9);

                DS=false;
                GP=false;
                GT=false;

                int index=0;

                while (pch_split != NULL)
                {
                    if(strcmp(pch_split,"DS")==0)
                    {
                        DS=true;
                        doseIndex=index;
                    }
                    if(strcmp(pch_split,"GP")==0)
                    {
                        GP=true;
                        gpIndex=index;
                    }
                    if(strcmp(pch_split,"GT")==0)
                    {
                        GT=true;
                        gtIndex=index;
                    }
                    index++;
                    pch_split = strtok_r (NULL, ":", &end_str9);
                }



                if(!DS && !GP && !GT)
                {
                    cout<<"\n Input VCF file has NONE of the following fields : GT, DS or GP "<<endl;
                    cout<<" VCF file must have at least one of these fields to work !!! "<<endl;
                    cout<<" Please check the following file : "<<filename<<endl;
                    return false;
                }

                int indCount=0;

                pch = strtok_r (NULL, "\t", &end_str1);
                while (pch != NULL)
                {

                    char *end_str2;
                    pch_split = strtok_r (pch,":", &end_str2);
                    index=0;

                    while (pch_split != NULL)
                    {
                        if(GP)
                        {
                            if(Format=="DS")
                                if(gpIndex==index)
                                {
                                    NumGP++;

                                    pch_split3 = strtok_r (pch_split,",",
                                                        &end_str3);



                                    double GP11=atof(pch_split3);
                                    pch_split3 = strtok_r (NULL,",", &end_str3);
                                    double GP22=atof(pch_split3);
                                    dosage[markerCount][indCount]=GP22+(2*(1-GP11-GP22));
                                }
                            if(Format=="GP")
                                if(gpIndex==index)
                                {
                                    NumGP++;

                                    pch_split3 = strtok_r (pch_split,",",
                                                        &end_str3);
                                    GP1[markerCount][indCount]=atof(pch_split3);
                                    pch_split3 = strtok_r (NULL,",", &end_str3);
                                    GP2[markerCount][indCount]=atof(pch_split3);
                                }

                        }
                        else if(DS)
                        {
                            if(Format=="DS")
                                if(doseIndex==index)
                                {
                                    dosage[markerCount][indCount]=atof(pch_split);
                                    NumDS++;
                                }
                            if(Format=="GP")
                            {
                                if(!GT)
                                {
                                    cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
                                    cout<<"                                    ERROR !!!                                  "<<endl;
                                    cout<<" ------------------------------------------------------------------------------"<<endl;

                                    std::cout << "\n           No Genotype or Genotype Probability found for variant "<<VariantList[totmarkerCount].name ;
                                    std::cout << "\n           Only dosage (DS) format found in dosage VCF file ...";
                                    std::cout << "\n           GP cannot be calculated without GT or GP in VCF file." << endl;
                                    return false;
                                }
                                else
                                {
                                    if(gtIndex==index)
                                    {
                                         NumGT++;
                                         end_str3=NULL;

                                         if((string)pch_split=="0|0")
                                         {
                                             GP1[markerCount][indCount]=1.0;
                                             GP2[markerCount][indCount]=0.0;
                                         }
                                         else if((string)pch_split=="0|1")
                                         {
                                             GP1[markerCount][indCount]=0.0;
                                             GP2[markerCount][indCount]=1.0;
                                         }
                                         else if((string)pch_split=="1|0")
                                         {
                                             GP1[markerCount][indCount]=0.0;
                                             GP2[markerCount][indCount]=1.0;
                                         }
                                         else if((string)pch_split=="1|1")
                                        {
                                             GP1[markerCount][indCount]=0.0;
                                             GP2[markerCount][indCount]=0.0;
                                         }
                                         else
                                         {
                                            cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
                                            cout<<"                                    ERROR !!!                                  "<<endl;
                                            cout<<" ------------------------------------------------------------------------------"<<endl;
                                            std::cout << "\n           Invalid value of GT ["<<(string)pch_split<<"] found for variant "<<VariantList[totmarkerCount].name ;
                                            std::cout << "\n           Only 0|0, 0|1 and 1|1 values are supported ...";
                                            std::cout << "\n           Please check input VCF file ..."<<endl;
                                            return false;
                                         }

                                    }
                                }
                            }
                        }
                        else if(GT)
                        {
                            if(Format=="DS")
                                if(gtIndex==index)
                                {
                                    NumGT++;
                                     end_str3=NULL;
                                     pch_split3 = strtok_r (pch_split,"|",
                                                            &end_str3);

                                    dosage[markerCount][indCount]=atoi(pch_split3);
                                    pch_split3 = strtok_r (NULL,"|", &end_str3);
                                    dosage[markerCount][indCount]+=atoi(pch_split3);
                                }
                            if(Format=="GP")
                            {
                                if(gtIndex==index)
                                {
                                     NumGT++;
                                     end_str3=NULL;
                                     if((string)pch_split=="0|0")
                                     {
                                         GP1[markerCount][indCount]=1.0;
                                         GP2[markerCount][indCount]=0.0;
                                     }
                                     else if((string)pch_split=="0|1")
                                     {
                                         GP1[markerCount][indCount]=0.0;
                                         GP2[markerCount][indCount]=1.0;
                                     }
                                     else if((string)pch_split=="1|0")
                                     {
                                         GP1[markerCount][indCount]=0.0;
                                         GP2[markerCount][indCount]=1.0;
                                     }
                                     else if((string)pch_split=="1|1")
                                    {
                                         GP1[markerCount][indCount]=0.0;
                                         GP2[markerCount][indCount]=0.0;
                                     }
                                     else
                                    {
                                            cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
                                            cout<<"                                    ERROR !!!                                  "<<endl;
                                            cout<<" ------------------------------------------------------------------------------"<<endl;
                                            std::cout << "\n           Invalid value of GT ["<<(string)pch_split<<"] found for variant "<<VariantList[totmarkerCount].name ;
                                            std::cout << "\n           Only 0|0, 0|1 and 1|1 values are supported ...";
                                            std::cout << "\n           Please check input VCF file ..."<<endl;
                                            return false;
                                    }

                                }
                            }
                        }


                        pch_split = strtok_r (NULL,":", &end_str2);
                        index++;
                    }

                    pch = strtok_r (NULL, "\t", &end_str1);
                    indCount++;
                }


                if(Type=="plink")
                {

                    ifprintf(plinkMap,"%s\t%s\t0\t%d\n",VariantList[totmarkerCount].chr.c_str(),
                             VariantList[totmarkerCount].name.c_str(),
                             VariantList[totmarkerCount].bp);
                    ifprintf(plinkMain,"%s\t%s\t%s",VariantList[totmarkerCount].name.c_str(),
                             VariantList[totmarkerCount].altAlleleString.c_str(),
                             VariantList[totmarkerCount].refAlleleString.c_str());

                    for (int i = 0; i < numSamplesRead; i++)
                    {
                        if(Format=="DS")
                            ifprintf(plinkMain,"\t%.3f",dosage[markerCount][i]);
                        else
                            ifprintf(plinkMain,"\t%.3f\t%.3f",GP1[markerCount][i],GP2[markerCount][i]);
                    }
                    ifprintf(plinkMain,"\n");
                }

                markerCount++;
                totmarkerCount++;
                line.clear();
                DoseRead->readLine(line);
            }

            if(Type=="plink")
            {
                if(totmarkerCount%factor==0)
                    printf("    Finished Writing %d markers to PLINK %s Format \n", totmarkerCount,Format=="DS"?"Dosage":"GP");
            }
            else
                printf("    Finished Reading Chunk %d of %d markers ... \n", part, bufferSize);

            if(Type=="mach")
            {

                machPartial = ifopen(outFile + ".output.part."+part + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
                for (int i = 0; i < numSamplesRead; i++)
                {
                    for (int j = 0; j < markerCount; j++)
                    {
                        if(Format=="DS")
                            ifprintf(machPartial,"\t%.3f",dosage[j][i]);
                        else
                            ifprintf(machPartial,"\t%.3f\t%.3f",GP1[j][i],GP2[j][i]);
                    }

                    ifprintf(machPartial,"\n");
                }
                ifclose(machPartial);

                if(Format=="DS")
                    printf("    Finished Writing Chunk %d of %d markers to MaCH Dosage Format ... \n", part, bufferSize);
                else
                    printf("    Finished Writing Chunk %d of %d markers to MaCH GP Format ... \n", part, bufferSize);

            }

            part++;



        }


    }
    else
    {
        cout<<"\n Following File File Not Available : "<<filename<<endl;
        return false;
    }

    ifclose(DoseRead);


    std::cout << "\n\n Number of Samples read from VCF Dosage File              : " << numSamplesRead ;
	std::cout << "\n Number of Markers read from VCF Dosage File              : " << totmarkerCount ;


	std::cout << "\n Number of Variants imported from GP (Genotype Prob)      : " << NumGP/numSamplesRead;
	std::cout << "\n Number of Variants imported from DS (Dosage)             : " << NumDS/numSamplesRead;
	std::cout << "\n Number of Variants imported from GT (Genotype)           : " << NumGT/numSamplesRead<<endl<<endl;


    if((NumGP+NumDS)==0)
    {
        if(Format=="DS")
        {
            cout<<" ------------------------------------------------------------------------------"<<endl;
            cout<<"                                   WARNING !!!                                 "<<endl;
            cout<<" ------------------------------------------------------------------------------"<<endl;
            std::cout << "\n WARNING !!! All Dosage values imported from GT values " ;
            std::cout << "\n             No GP or DS values were found in VCF file ! " ;
            std::cout << "\n             GT values (although may be available,  " ;
            std::cout << "\n             might not have been meant to be used ";
                cout<<   "\n             for DS calculation ! " ;

        }
        else
        {
            cout<<" ------------------------------------------------------------------------------"<<endl;
            cout<<"                                   WARNING !!!                                 "<<endl;
            cout<<" ------------------------------------------------------------------------------"<<endl;
            std::cout << "\n WARNING !!! All Genotype Probabilities calculated from GT values " ;
            std::cout << "\n             No GP values were found in VCF file ! " ;
            std::cout << "\n             GT values (although may be available,  " ;
            std::cout << "\n             might not have been meant to be used ";
                cout<<   "\n             for GP calculation ! " ;
        }
        std::cout << "\n             Please Verify Dosage File ! " << endl;
    }



    IFILE machFinal;
    if(Type=="mach")
    {
        cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
        cout<<"                             MERGE MACH OUTPUT FILE                            "<<endl;
        cout<<" ------------------------------------------------------------------------------"<<endl<<endl;

        if(Format=="DS")
        {
            printf("\n Merging Split Files to MaCH Dosage File : %s \n\n",(outFile + ".mach.dose" + (gzip ? ".gz" : "")).c_str());
            machFinal = ifopen(outFile + ".mach.dose" + (gzip ? ".gz" : ""),  "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
        }
        else
        {
            printf("\n Merging Split Files to MaCH GP File : %s \n\n",(outFile + ".mach.gprob" + (gzip ? ".gz" : "")).c_str());
            machFinal = ifopen(outFile + ".mach.gprob" + (gzip ? ".gz" : ""),  "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
        }

        vector<IFILE> machFinalList(part);

        for(int i=0;i<part;i++)
        {
            string tempFileIndex(outFile);
            stringstream strs;
            strs<<(i);

            tempFileIndex+=(".output.part." +
                            (string)(strs.str())+(gzip ? ".gz" : ""));
            machFinalList[i] = ifopen(tempFileIndex.c_str(), "r");
        }

        for(int j=0;j<numSamplesRead;j++)
        {
            if(j%500==0)
            {
                  printf("     Merging Sample %d of %d to MaCH File ...", j + 1, numSamplesRead);
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

        printf("\n Creating Info file : %s \n\n",(outFile + ".mach.info").c_str());


        machFinal = ifopen(outFile + ".mach.info",  "wb", InputFile::UNCOMPRESSED);

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

                ifprintf(machFinal,"%s\n",line.c_str());
                Infocount++;
                line.clear();
                Indic=ifs->readLine(line);
            }
            if(Infocount!=numMarkers)
            {
                cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
                cout<<"                                    ERROR !!!                                  "<<endl;
                cout<<" ------------------------------------------------------------------------------"<<endl;
              std::cout << "\n\n             "<<Infocount<<" markers found in Input info file <"<< Info <<"> " ;
                std::cout << "\n             However, "<<numMarkers <<" markers were read from dosage VCF file <"<< filename<< "> " ;
                std::cout << "\n             Please use info file from same imputation run !!! " << endl;
                remove((outFile + ".mach." + (Format=="DS" ? "dose" : "gprob")  + (gzip ? ".gz" : "")).c_str());
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


            ifprintf(machFinal, "SNP\tREF(0)\tALT(1)\tALT_Frq\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose0\tDose1\n");

            for (int i =0; i < totmarkerCount; i++)
            {
                variant &tempVariant=VariantList[i];



            ifprintf(machFinal, "%s\t%s\t%s\t%s\t%s\t-\t%s\t%s\t-\t-\t%s\t-\t-\n",
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

        ifclose(machFinal);



        cout<<" \n [NOTE: In output MaCH files, dosage values are twice the probability of";
        printf("\n        ALTERNATE ALLELE (in the info file) and NOT twice the probability");
        cout<< "\n        of MAJOR ALLELE (as is in usual Minimac output)"<<endl;
        cout<<   "        Similarly, Genotype probabilities are given for REF|REF followed"<<endl;
        cout<<   "        REF|ALT (allele given in info file)."<<endl<<endl;


        if(Format=="DS")
        {
            printf("\n Dosage Information written to MaCH Dosage File  : %s ",(outFile + ".mach.dose" + (gzip ? ".gz" : "")).c_str());
            printf("\n Information written to MaCH Info File           : %s ",(outFile + ".mach.info").c_str());

        }
        else
        {
            printf("\n Genotype Probability written to MaCH Dosage File  : %s ",(outFile + ".mach.gprob" + (gzip ? ".gz" : "")).c_str());
            printf("\n Information written to MaCH Info File             : %s ",(outFile + ".mach.info").c_str());
        }

    }
    else
    {

        cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
        cout<<"                                PLINK OUTPUT FILE                            "<<endl;
        cout<<" ------------------------------------------------------------------------------"<<endl<<endl;

        cout<<" \n [NOTE: In output PLINK files, dosage values are twice the probability of";
        printf("\n        A1 ALLELE (in the dosage/gprob file) and NOT twice the probability");
        cout<< "\n        of MAJOR ALLELE (as is in usual Minimac output)"<<endl;
        cout<<   "        Similarly, Genotype probabilities are given for A1|A1 followed"<<endl;
        cout<<   "        A1|A2 (alleles given in dosage/gprob file)."<<endl<<endl;



        if(Format=="DS")
        {
            printf("\n Dosage Information written to PLINK Dosage File  : %s ",(outFile + ".plink.dose" + (gzip ? ".gz" : "")).c_str());
            printf("\n Map Information written to PLINK Map File        : %s ",(outFile + ".plink.map").c_str());
            printf("\n Family Information written to PLINK Fam File     : %s ",(outFile + ".plink.fam").c_str());
        }
        else
        {
            printf("\n Genotype Probability written to PLINK Dosage File  : %s ",(outFile + ".plink.gprob" + (gzip ? ".gz" : "")).c_str());
            printf("\n Map Information written to PLINK Map File          : %s ",(outFile + ".plink.map").c_str());
            printf("\n Family Information written to PLINK Fam File       : %s ",(outFile + ".plink.fam").c_str());
        }
        ifclose(plinkMain);
        ifclose(plinkMap);
    }

    cout<<endl;

    return true;

}




bool HaplotypeSet::LoadInfoFromVCFfile()
{

	VcfFileReader inFile;
	VcfHeader header;
	VcfRecord record;
    individualName.clear();
    familyName.clear();
    String filename=VcfDose;
    int TypedOnlyCount=0,GenotypedCount=0,ImputeCount=0;
    int numSamplesRead = 0;
    if (!inFile.open(filename, header))
	{
		cout << "\n Program could NOT open file : " << filename << endl<<endl;
		return false;
	}

    std::cout << "\n Calculating number of Samples and Markers in VCF File ..." ;

    inFile.setSiteOnly(true);
	numSamplesRead = header.getNumSamples();
    numSamples=numSamplesRead;

    if(numSamplesRead==0)
    {
        cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
        cout<<"                                    ERROR !!!                                  "<<endl;
        cout<<" ------------------------------------------------------------------------------"<<endl;
        std::cout << "\n Number of Samples read from VCF File    : " << numSamplesRead << endl;
        std::cout << "\n ERROR !!! "<<endl;
        cout << "\n NO samples found in VCF File !! \n Please Check Input File !!!  "<< endl;
        return false;
    }

    for (int i = 0; i < numSamplesRead; i++)
	{
		string tempName(header.getSampleName(i));

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
                    cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
                    cout<<"                                    ERROR !!!                                  "<<endl;
                    cout<<" ------------------------------------------------------------------------------"<<endl;
                    cout << "\n Program could NOT parse the following sample name with ID delimiter ["<< IdDelimiter<<"] : " << temptempName << endl;
                    cout<<    " More than TWO tokens found : ["<<familyName.back()  <<"] ["
                    <<token << "] ["<<tempName  <<"] " <<endl;
                    cout << "\n Please verify ID Delimitier : [" << IdDelimiter <<"]"<< endl;
                    return false;
                }

                Count++;


            }
            if(Count==0)
            {
                cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
                cout<<"                                    ERROR !!!                                  "<<endl;
                cout<<" ------------------------------------------------------------------------------"<<endl;
                cout << "\n Program could NOT parse the following sample name with ID delimiter ["<< IdDelimiter<<"] : " << temptempName << endl;
                cout<<    " Delimiter NOT FOUND in Sample Namesss !!! " <<endl;
                cout << "\n Please verify ID Delimitier : [" << IdDelimiter <<"]"<< endl;
                return false;
            }

            individualName.push_back(tempName);
        }
    }

    VariantList.clear();
    string name,refAlleleString,altAlleleString;
    int bp;
    string chr;

    while (inFile.readRecord(record))
	{
		chr=record.getChromStr();
        bp=record.get1BasedPosition();
        name=record.getIDStr();

        if(!CheckValidChrom(chr))
        {
            cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
            cout<<"                                    ERROR !!!                                  "<<endl;
            cout<<" ------------------------------------------------------------------------------"<<endl;
          cout << "\n Reference VCF File contains chromosome : "<<chr<<endl;
            cout << " VCF File can only contain chromosomes 1-22 !!! "<<endl;
            cout << " Please contact author sayantan@umich.edu for help if you need to use other chromosomes ... "<<endl;
            cout << " Program Aborting ... "<<endl;
            return false;
        }




        refAlleleString = record.getRefStr();
		altAlleleString = record.getAltStr();

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



        if(ThisInfo.getString("AF")!=NULL)
            tempVariant.Af=*ThisInfo.getString("AF");

        if(ThisInfo.getString("AF")!=NULL)
            tempVariant.Af=*ThisInfo.getString("AF");

        if(ThisInfo.getString("MAF")!=NULL)
            tempVariant.Maf=*ThisInfo.getString("MAF");

        if(ThisInfo.getString("R2")!=NULL)
            tempVariant.Rsq=*ThisInfo.getString("R2");
        else
        {
            tempVariant.genotyped_only=true;
            tempVariant.Tag="Typed_Only";
        }

        if(ThisInfo.getString("ER2")!=NULL)
        {
            tempVariant.Ersq=*ThisInfo.getString("ER2");
            tempVariant.genotyped=true;
            tempVariant.Tag="Genotyped";
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



  std::cout << "\n\n Number of Samples read from VCF Dosage File              : " << numSamplesRead ;
    std::cout << "\n Number of Markers read from VCF Dosage File              : " << numMarkers <<endl<<endl;


  std::cout << "\n\n Number of Imputed Variants                               : " << ImputeCount;
	std::cout << "\n Number of Genotyped and Imputed Variants                 : " << GenotypedCount;
	std::cout << "\n Number of Genotyped Only Variants                        : " << TypedOnlyCount<<endl<<endl;


    return true;


}


bool HaplotypeSet::CheckValidChrom(string chr)
{
    bool result=false;

    if(MyChromosome!="" && chr==MyChromosome)
        return true;

    string temp[]={"1","2","3","4","5","6","7","8","9","10","11"
                    ,"12","13","14","15","16","17","18","19","20","21","22"};
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

        cout << "\n Dosage File provided by \"--vcfDose\" must be a VCF file !!! \n";
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

//
//bool HaplotypeSet::ConvertDosageDataOld()
//{
//
//	VcfFileReader inFile;
//	VcfHeader header;
//	VcfRecord record;
//    individualName.clear();
//    familyName.clear();
//    String filename=VcfDose;
//    int numSamplesRead = 0,markerCount=0,totmarkerCount=0,NumGP=0,NumDS=0,NumGT=0;
//    if (!inFile.open(filename, header))
//	{
//		cout << "\n Program could NOT open file : " << filename << endl<<endl;
//		return false;
//	}
//
//    inFile.setSiteOnly(false);
//	numSamplesRead = header.getNumSamples();
//    numSamples=numSamplesRead;
//    char * pch_split,* pch_split3;
//    char * pch;
//    char *end_str1,*end_str3;
//
//
//    if(numSamplesRead==0)
//    {
//        cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
//        cout<<"                                    ERROR !!!                                  "<<endl;
//        cout<<" ------------------------------------------------------------------------------"<<endl;
//
//        std::cout << "\n Number of Samples read from VCF File    : " << numSamplesRead << endl;
//        cout << "\n NO samples found in VCF File !! \n Please Check Input File !!!  "<< endl;
//        return false;
//    }
//
//    for (int i = 0; i < numSamplesRead; i++)
//	{
//		string tempName(header.getSampleName(i));
//
//		if(IdDelimiter=="")
//        {
//            individualName.push_back(tempName);
//            familyName.push_back(tempName);
//        }
//		else
//        {
//
//            size_t pos = 0;
//            std::string delimiter(IdDelimiter) ;
//            std::string token;
//            int Count=0;
//            string temptempName=tempName;
//
//            while ((pos = tempName.find(delimiter)) != std::string::npos)
//            {
//                token = tempName.substr(0, pos);
//
//                if(Count==0)
//                    familyName.push_back(token);
//
//                tempName.erase(0, pos + delimiter.length());
//                if(Count==1)
//                {
//                    cout<<" ------------------------------------------------------------------------------"<<endl;
//                    cout<<"                                    ERROR !!!                                  "<<endl;
//                    cout<<" ------------------------------------------------------------------------------"<<endl;
//
//                    cout << "\n Program could NOT parse the following sample name with ID delimiter ["<< IdDelimiter<<"] : " << temptempName << endl;
//                    cout<<    " More than TWO tokens found : ["<<familyName.back()  <<"] ["
//                    <<token << "] ["<<tempName  <<"] " <<endl;
//                    cout << "\n Please verify ID Delimitier : [" << IdDelimiter <<"]"<< endl;
//                    return false;
//                }
//
//                Count++;
//
//
//            }
//            if(Count==0)
//            {
//                cout<<" ------------------------------------------------------------------------------"<<endl;
//                cout<<"                                    ERROR !!!                                  "<<endl;
//                cout<<" ------------------------------------------------------------------------------"<<endl;
//                cout << "\n Program could NOT parse the following sample name with ID delimiter ["<< IdDelimiter<<"] : " << temptempName << endl;
//                cout<<    " Delimiter NOT FOUND in Sample Namessss !!! " <<endl;
//                cout << "\n Please verify ID Delimitier : [" << IdDelimiter <<"]"<< endl;
//                return false;
//            }
//
//            individualName.push_back(tempName);
//        }
//    }
//
//    IFILE machPartial=NULL,plinkMain=NULL,plinkMap=NULL;
//    int part=1;
//
//     std::cout << "\n Number of Samples read from VCF File    : " << numSamplesRead << endl;
//
//
//    if(Type=="mach")
//        machPartial = ifopen(outFile + ".output.part.0" + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
//    else
//    {
//        plinkMain = ifopen(outFile + ".plink.fam", "wb",InputFile::UNCOMPRESSED);
//        for (int i = 0; i < numSamplesRead; i++)
//        {
//             ifprintf(plinkMain,"%s\t%s\t0\t0\t0\t-9\n",familyName[i].c_str(),individualName[i].c_str());
//        }
//        ifclose(plinkMain);
//        plinkMap = ifopen(outFile + ".plink.map", "wb",InputFile::UNCOMPRESSED);
//        plinkMain = ifopen(outFile + ".plink."+(Format=="DS" ? "dose" : "gprob") + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
//
//        ifprintf(plinkMain,"SNP\tA1\tA2");
//        for (int i = 0; i < numSamplesRead; i++)
//        {
//             ifprintf(plinkMain,"\t%s\t%s",familyName[i].c_str(),individualName[i].c_str());
//        }
//        ifprintf(plinkMain,"\n");
//    }
//
//
//    if(Type=="mach")
//        for (int i = 0; i < numSamplesRead; i++)
//        {
//             ifprintf(machPartial,"%s->%s\tDOSE\n",familyName[i].c_str(),individualName[i].c_str());
//        }
//
//    if(Type=="mach")
//        ifclose(machPartial);
//
////    VariantList.clear();
//    string name,refAlleleString,altAlleleString;
//    int bp;
//    string chr;
//
//    int bufferSize=BufferSize;
//    if(Type=="plink")
//        bufferSize=1;
//
//    dosage.clear();
//    GP1.clear();
//    GP2.clear();
//
//    if(Format=="DS")
//    {
//        dosage.resize(bufferSize);
//        for (int i = 0; i < bufferSize; i++)
//        {
//             dosage[i].resize(numSamplesRead,-9.0);
//        }
//    }
//
//
//    if(Format=="GP")
//    {
//        GP1.resize(bufferSize);
//        GP2.resize(bufferSize);
//        for (int i = 0; i < bufferSize; i++)
//        {
//            GP1[i].resize(numSamplesRead,-9.0);
//            GP2[i].resize(numSamplesRead,-9.0);
//        }
//    }
//    int factor=10000;
//
//    IFILE DoseRead = ifopen(filename, "r");
//    string line;
//
//    std::cout << "\n Reading and Importing Data from VCF File ..."<<endl<<endl ;
//
//    if(DoseRead)
//    {
//        bool Header=true;
//
//        while(Header)
//        {
//            line.clear();
//            DoseRead->readLine(line);
//            if(line.substr(0,1).compare("#")==0)
//                Header=true;
//            else
//                Header=false;
//        }
//
//        while(line.compare("")!=0)
//        {
//            markerCount=0;
//
//            while(markerCount<bufferSize && line.compare("")!=0)
//            {
//                if(totmarkerCount>=numMarkers)
//                {
//                    cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
//                    cout<<"                                    ERROR !!!                                  "<<endl;
//                    cout<<" ------------------------------------------------------------------------------"<<endl;
//
//                    std::cout << "\n           Number of markers found in info file = "<<numMarkers ;
//                    std::cout << "\n           More than "<<numMarkers <<" markers are found in dosage VCF file ... " ;
//                    std::cout << "\n           Please use info file from same imputation run !!! " << endl;
//                    return false;
//                }
//
//
//                char *end_str9;
//                pch = strtok_r ((char*)line.c_str(),"\t", &end_str1);
//                if(pch==NULL)
//                {
//                    cout<<"\n Input VCF file has invalid format at Marker No : "<<totmarkerCount+1<<endl;
//                    cout<<" Please verify the following file : "<<filename<<endl;
//                    return false;
//                }
//                else
//                    chr=pch;
//                pch = strtok_r (NULL, "\t", &end_str1);
//
//                if(pch==NULL)
//                {
//                    cout<<"\n Input VCF file has invalid format at Marker No : "<<totmarkerCount+1<<endl;
//                    cout<<" Please verify the following file : "<<filename<<endl;
//                    return false;
//                }
//                else
//                    bp=atoi(pch);
//
//                pch = strtok_r (NULL, "\t", &end_str1);
//                if(pch==NULL)
//                {
//                    cout<<"\n Input VCF file has invalid format at Marker No : "<<totmarkerCount+1<<endl;
//                    cout<<" Please verify the following file : "<<filename<<endl;
//                    return false;
//                }
//                else
//                    name=pch;
//
//                if(VariantList[totmarkerCount].name!=name)
//                {
//                    std::cout << "\n ERROR !!! Mismatching variant name between VCF and INFO file !!!" ;
//                    std::cout << "\n           Marker #"<<totmarkerCount+1 <<" ["<<name<<"] in VCF file does NOT match with"<<
//                    " Marker #"<<totmarkerCount+1 <<" ["<<VariantList[totmarkerCount].name <<"] in INFO file" ;
//                    std::cout << "\n           Please use info file from same imputation run !!! " << endl;
//                    return false;
//                }
//
//
//                variant tempVariant(name,chr,bp);
//
//                pch = strtok_r (NULL, "\t", &end_str1);
//                pch = strtok_r (NULL, "\t", &end_str1);
//
//                VariantList[totmarkerCount].assignValues(name,chr,bp);
//
//                pch = strtok_r (NULL, "\t", &end_str1);
//                pch = strtok_r (NULL, "\t", &end_str1);
//                pch = strtok_r (NULL, "\t", &end_str1);
//                pch = strtok_r (NULL, "\t", &end_str1);
//
//                pch_split = strtok_r (pch,":", &end_str9);
//
//                DS=false;
//                GP=false;
//                GT=false;
//
//                int index=0;
//
//                while (pch_split != NULL)
//                {
//                    if(strcmp(pch_split,"DS")==0)
//                    {
//                        DS=true;
//                        doseIndex=index;
//                    }
//                    if(strcmp(pch_split,"GP")==0)
//                    {
//                        GP=true;
//                        gpIndex=index;
//                    }
//                    if(strcmp(pch_split,"GT")==0)
//                    {
//                        GT=true;
//                        gtIndex=index;
//                    }
//                    index++;
//                    pch_split = strtok_r (NULL, ":", &end_str9);
//                }
//
//
//
//                if(!DS && !GP && !GT)
//                {
//                    cout<<"\n Input VCF file has NONE of the following fields : GT, DS or GP "<<endl;
//                    cout<<" VCF file must have at least one of these fields to work !!! "<<endl;
//                    cout<<" Please check the following file : "<<filename<<endl;
//                    return false;
//                }
//
//                int indCount=0;
//
//                pch = strtok_r (NULL, "\t", &end_str1);
//                while (pch != NULL)
//                {
//
//                    char *end_str2;
//                    pch_split = strtok_r (pch,":", &end_str2);
//                    index=0;
//
//                    while (pch_split != NULL)
//                    {
//                        if(GP)
//                        {
//                            if(Format=="DS")
//                                if(gpIndex==index)
//                                {
//                                    NumGP++;
//
//                                    pch_split3 = strtok_r (pch_split,",",
//                                                        &end_str3);
//
//
//
//                                    double GP11=atof(pch_split3);
//                                    pch_split3 = strtok_r (NULL,",", &end_str3);
//                                    double GP22=atof(pch_split3);
//                                    dosage[markerCount][indCount]=GP22+(2*(1-GP11-GP22));
//                                }
//                            if(Format=="GP")
//                                if(gpIndex==index)
//                                {
//                                    NumGP++;
//
//                                    pch_split3 = strtok_r (pch_split,",",
//                                                        &end_str3);
//                                    GP1[markerCount][indCount]=atof(pch_split3);
//                                    pch_split3 = strtok_r (NULL,",", &end_str3);
//                                    GP2[markerCount][indCount]=atof(pch_split3);
//                                }
//
//                        }
//                        else if(DS)
//                        {
//                            if(Format=="DS")
//                                if(doseIndex==index)
//                                {
//                                    dosage[markerCount][indCount]=atof(pch_split);
//                                    NumDS++;
//                                }
//                            if(Format=="GP")
//                            {
//                                if(!GT)
//                                {
//                                    std::cout << "\n ERROR !!! No Genotype or Genotype Probability found for variant "<<VariantList[totmarkerCount].name ;
//                                    std::cout << "\n           Only dosage (DS) format found in dosage VCF file ...";
//                                    std::cout << "\n           GP cannot be calculated without GT or GP in VCF file." << endl;
//                                    return false;
//                                }
//                                else
//                                {
//                                    if(gtIndex==index)
//                                    {
//                                         NumGT++;
//                                         end_str3=NULL;
//
//                                         if((string)pch_split=="0|0")
//                                         {
//                                             GP1[markerCount][indCount]=1.0;
//                                             GP2[markerCount][indCount]=0.0;
//                                         }
//                                         else if((string)pch_split=="0|1")
//                                         {
//                                             GP1[markerCount][indCount]=0.0;
//                                             GP2[markerCount][indCount]=1.0;
//                                         }
//                                         else if((string)pch_split=="1|0")
//                                         {
//                                             GP1[markerCount][indCount]=0.0;
//                                             GP2[markerCount][indCount]=1.0;
//                                         }
//                                         else if((string)pch_split=="1|1")
//                                        {
//                                             GP1[markerCount][indCount]=0.0;
//                                             GP2[markerCount][indCount]=0.0;
//                                         }
//                                         else
//                                         {
//                                            std::cout << "\n ERROR !!! Invalid value of GT ["<<(string)pch_split<<"] found for variant "<<VariantList[totmarkerCount].name ;
//                                            std::cout << "\n           Only 0|0, 0|1 and 1|1 values are supported ...";
//                                            std::cout << "\n           Please check input VCF file ..."<<endl;
//                                            return false;
//                                         }
//
//                                    }
//                                }
//                            }
//                        }
//                        else if(GT)
//                        {
//                            if(Format=="DS")
//                                if(gtIndex==index)
//                                {
//                                    NumGT++;
//                                     end_str3=NULL;
//                                     pch_split3 = strtok_r (pch_split,"|",
//                                                            &end_str3);
//
//                                    dosage[markerCount][indCount]=atoi(pch_split3);
//                                    pch_split3 = strtok_r (NULL,"|", &end_str3);
//                                    dosage[markerCount][indCount]+=atoi(pch_split3);
//                                }
//                            if(Format=="GP")
//                            {
//                                if(gtIndex==index)
//                                {
//                                     NumGT++;
//                                     end_str3=NULL;
//                                     if((string)pch_split=="0|0")
//                                     {
//                                         GP1[markerCount][indCount]=1.0;
//                                         GP2[markerCount][indCount]=0.0;
//                                     }
//                                     else if((string)pch_split=="0|1")
//                                     {
//                                         GP1[markerCount][indCount]=0.0;
//                                         GP2[markerCount][indCount]=1.0;
//                                     }
//                                     else if((string)pch_split=="1|0")
//                                     {
//                                         GP1[markerCount][indCount]=0.0;
//                                         GP2[markerCount][indCount]=1.0;
//                                     }
//                                     else if((string)pch_split=="1|1")
//                                    {
//                                         GP1[markerCount][indCount]=0.0;
//                                         GP2[markerCount][indCount]=0.0;
//                                     }
//                                     else
//                                    {
//                                            std::cout << "\n ERROR !!! Invalid value of GT ["<<(string)pch_split<<"] found for variant "<<VariantList[totmarkerCount].name ;
//                                            std::cout << "\n           Only 0|0, 0|1 and 1|1 values are supported ...";
//                                            std::cout << "\n           Please check input VCF file ..."<<endl;
//                                            return false;
//                                    }
//
//                                }
//                            }
//                        }
//
//
//                        pch_split = strtok_r (NULL,":", &end_str2);
//                        index++;
//                    }
//
//                    pch = strtok_r (NULL, "\t", &end_str1);
//                    indCount++;
//                }
//
//
//                if(Type=="plink")
//                {
//
//                    ifprintf(plinkMap,"%s\t%s\t0\t%d\n",VariantList[totmarkerCount].chr.c_str(),
//                             VariantList[totmarkerCount].name.c_str(),
//                             VariantList[totmarkerCount].bp);
//                    ifprintf(plinkMain,"%s\t%s\t%s",VariantList[totmarkerCount].name.c_str(),
//                             VariantList[totmarkerCount].altAlleleString.c_str(),
//                             VariantList[totmarkerCount].refAlleleString.c_str());
//
//                    for (int i = 0; i < numSamplesRead; i++)
//                    {
//                        if(Format=="DS")
//                            ifprintf(plinkMain,"\t%.3f",dosage[markerCount][i]);
//                        else
//                            ifprintf(plinkMain,"\t%.3f\t%.3f",GP1[markerCount][i],GP2[markerCount][i]);
//                    }
//                    ifprintf(plinkMain,"\n");
//                }
//
//                markerCount++;
//                totmarkerCount++;
//                line.clear();
//                DoseRead->readLine(line);
//            }
//
//            if(Type=="plink")
//            {
//                if(totmarkerCount%factor==0)
//                    printf("    Finished Writing %d markers to PLINK %s Format \n", totmarkerCount,Format=="DS"?"Dosage":"GP");
//            }
//            else
//                printf("    Finished Reading Chunk %d of %d markers ... \n", part, bufferSize);
//
//            if(Type=="mach")
//            {
//
//                machPartial = ifopen(outFile + ".output.part."+part + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
//                for (int i = 0; i < numSamplesRead; i++)
//                {
//                    for (int j = 0; j < markerCount; j++)
//                    {
//                        if(Format=="DS")
//                            ifprintf(machPartial,"\t%.3f",dosage[j][i]);
//                        else
//                            ifprintf(machPartial,"\t%.3f\t%.3f",GP1[j][i],GP2[j][i]);
//                    }
//
//                    ifprintf(machPartial,"\n");
//                }
//                ifclose(machPartial);
//
//                if(Format=="DS")
//                    printf("    Finished Writing Chunk %d of %d markers to MaCH Dosage Format ... \n", part, bufferSize);
//                else
//                    printf("    Finished Writing Chunk %d of %d markers to MaCH GP Format ... \n", part, bufferSize);
//
//            }
//
//            part++;
//
//
//
//        }
//
//
//    }
//    else
//    {
//        cout<<"\n Following File File Not Available : "<<filename<<endl;
//        return false;
//    }
//
//
//
//
//
//
//    std::cout << "\n\n Number of Samples read from VCF Dosage File              : " << numSamplesRead ;
//	std::cout << "\n Number of Markers read from VCF Dosage File              : " << totmarkerCount ;
//
//    if(totmarkerCount<numMarkers)
//    {
//        std::cout << "\n\n ERROR !!! Number of markers found in info file = "<<numMarkers ;
//        std::cout << "\n           Less than "<<numMarkers <<" markers were read from dosage VCF file ... " ;
//        std::cout << "\n           Please use info file from same imputation run !!! " << endl;
//        return false;
//    }
//	std::cout << "\n Number of Variants imported from GP (Genotype Prob)      : " << NumGP/numSamplesRead;
//	std::cout << "\n Number of Variants imported from DS (Dosage)             : " << NumDS/numSamplesRead;
//	std::cout << "\n Number of Variants imported from GT (Genotype)           : " << NumGT/numSamplesRead<<endl<<endl;
//
//
//    if((NumGP+NumDS)==0)
//    {
//        if(Format=="DS")
//        {
//            std::cout << "\n WARNING !!! All Dosage values imported from GT values " ;
//            std::cout << "\n             No GP or DS values were found in VCF file ! " ;
//        }
//        else
//        {
//            std::cout << "\n WARNING !!! All Genotype Probabilities calculated from GT values " ;
//            std::cout << "\n             No GP values were found in VCF file ! " ;
//            std::cout << "\n             GT values (although may be available,  " ;
//            std::cout << "\n             cannot be used for GP calculation ! " ;
//        }
//        std::cout << "\n             Please Verify Dosage File ! " << endl;
//    }
//
//
//
//    IFILE machFinal;
//    if(Type=="mach")
//    {
//        cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
//        cout<<"                             MERGE MACH OUTPUT FILE                            "<<endl;
//        cout<<" ------------------------------------------------------------------------------"<<endl<<endl;
//
//        if(Format=="DS")
//        {
//            printf("\n Merging Split Files to MaCH Dosage File : %s \n\n",(outFile + ".mach.dose" + (gzip ? ".gz" : "")).c_str());
//            machFinal = ifopen(outFile + ".mach.dose" + (gzip ? ".gz" : ""),  "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
//        }
//        else
//        {
//            printf("\n Merging Split Files to MaCH GP File : %s \n\n",(outFile + ".mach.gprob" + (gzip ? ".gz" : "")).c_str());
//            machFinal = ifopen(outFile + ".mach.gprob" + (gzip ? ".gz" : ""),  "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
//        }
//
//        vector<IFILE> machFinalList(part);
//
//        for(int i=0;i<part;i++)
//        {
//            string tempFileIndex(outFile);
//            stringstream strs;
//            strs<<(i);
//
//            tempFileIndex+=(".output.part." +
//                            (string)(strs.str())+(gzip ? ".gz" : ""));
//            machFinalList[i] = ifopen(tempFileIndex.c_str(), "r");
//        }
//
//        for(int j=0;j<numSamplesRead;j++)
//        {
//            if(j%500==0)
//            {
//                  printf("     Merging Sample %d of %d to MaCH File ...", j + 1, numSamplesRead);
//                cout<<endl;
//            }
//
//            for(int i=0;i<part;i++)
//            {
//                line.clear();
//                machFinalList[i]->readLine(line);
//                ifprintf(machFinal,"%s",line.c_str());
//            }
//            ifprintf(machFinal,"\n");
//        }
//        ifclose(machFinal);
//
//        for(int i=0;i<part;i++)
//        {
//            ifclose(machFinalList[i]);
//            string tempFileIndex(outFile);
//            stringstream strs;
//            strs<<(i);
//
//            tempFileIndex+=(".output.part." +
//                            (string)(strs.str())+(gzip ? ".gz" : ""));
//            remove(tempFileIndex.c_str());
//
//        }
//
//
//        cout<<" \n [NOTE: In output MaCH files, dosage values are twice the probability of";
//        printf("\n        ALTERNATE ALLELE (in the info file) and NOT twice the probability");
//        cout<< "\n        of MAJOR ALLELE (as is in usual Minimac output)"<<endl;
//        cout<<   "        Similarly, Genotype probabilities are given for REF|REF followed"<<endl;
//        cout<<   "        REF|ALT (allele given in info file)."<<endl<<endl;
//
//
//
//        IFILE ifs = ifopen(Info, "r");
//        machFinal = ifopen(outFile + ".mach.info",  "wb", InputFile::UNCOMPRESSED);
//
//
//        line.clear();
//        int Indic=ifs->readLine(line);
//
//        while(Indic==-1 ||  Indic==0)
//        {
//            if(Indic==-1 && line.length()==0)
//            {
//                break;
//            }
//
//            ifprintf(machFinal,"%s\n",line.c_str());
//            line.clear();
//            Indic=ifs->readLine(line);
//
//        }
//        ifclose(machFinal);
//
//
//        if(Format=="DS")
//        {
//            printf("\n Dosage Information written to MaCH Dosage File  : %s ",(outFile + ".mach.dose" + (gzip ? ".gz" : "")).c_str());
//            printf("\n Information written to MaCH Info File           : %s ",(outFile + ".mach.info").c_str());
//
//        }
//        else
//        {
//            printf("\n Genotype Probability written to MaCH Dosage File  : %s ",(outFile + ".mach.gprob" + (gzip ? ".gz" : "")).c_str());
//            printf("\n Information written to MaCH Info File             : %s ",(outFile + ".mach.info").c_str());
//        }
//
//    }
//    else
//    {
//
//        cout<<"\n\n ------------------------------------------------------------------------------"<<endl;
//        cout<<"                                PLINK OUTPUT FILE                            "<<endl;
//        cout<<" ------------------------------------------------------------------------------"<<endl<<endl;
//
//        cout<<" \n [NOTE: In output PLINK files, dosage values are twice the probability of";
//        printf("\n        A1 ALLELE (in the dosage/gprob file) and NOT twice the probability");
//        cout<< "\n        of MAJOR ALLELE (as is in usual Minimac output)"<<endl;
//        cout<<   "        Similarly, Genotype probabilities are given for A1|A1 followed"<<endl;
//        cout<<   "        A1|A2 (alleles given in dosage/gprob file)."<<endl<<endl;
//
//
//
//        if(Format=="DS")
//        {
//            printf("\n Dosage Information written to PLINK Dosage File  : %s ",(outFile + ".plink.dose" + (gzip ? ".gz" : "")).c_str());
//            printf("\n Map Information written to PLINK Map File        : %s ",(outFile + ".plink.map").c_str());
//            printf("\n Family Information written to PLINK Fam File     : %s ",(outFile + ".plink.fam").c_str());
//        }
//        else
//        {
//            printf("\n Genotype Probability written to PLINK Dosage File  : %s ",(outFile + ".plink.gprob" + (gzip ? ".gz" : "")).c_str());
//            printf("\n Map Information written to PLINK Map File          : %s ",(outFile + ".plink.map").c_str());
//            printf("\n Family Information written to PLINK Fam File       : %s ",(outFile + ".plink.fam").c_str());
//        }
//        ifclose(plinkMain);
//        ifclose(plinkMap);
//    }
//
//    cout<<endl;
//
//    return true;
//
//}
//
//bool HaplotypeSet::LoadInfoFile(String filename)
//{
//
//    IFILE ifs = ifopen(filename, "r");
//
//    VariantList.clear();
//    string name,refAlleleString,altAlleleString;
////    string MajAlleleString,MinAlleleString;
//    int bp=100;
//    string chr="Z";
//    char *pch,*end_str1;
//
//
//
//    int RowNo=1;
//
//    string line;
//    if(ifs)
//    {
//        ifs->readLine(line);
//        line.clear();
//        int Indic=ifs->readLine(line);
//
//        while(Indic==-1 ||  Indic==0)
//        {
//
//            if(Indic==-1 && line.length()==0)
//            {
//                break;
//            }
//
//
//            RowNo++;
//            pch = strtok_r ((char*)line.c_str(),"\t", &end_str1);
//
//            if(pch==NULL)
//            {
//                cout<<"\n Info file does NOT have 14 columns at Row : "<<RowNo<<endl;
//                cout<<" Please verify the following file : "<<filename<<endl;
//                return false;
//            }
//            else
//                name=pch;
//
//            variant tempVariant(name,chr,bp);
//
//
//            pch = strtok_r (NULL, "\t", &end_str1);
//            if(pch==NULL)
//            {
//                cout<<"\n Info file does NOT have 14 columns at Row : "<<RowNo<<endl;
//                cout<<" Please verify the following file : "<<filename<<endl;
//                return false;
//            }
//            else
//                refAlleleString=pch;
//
//            pch = strtok_r (NULL, "\t", &end_str1);
//            if(pch==NULL)
//            {
//                cout<<"\n Info file does NOT have 14 columns at Row : "<<RowNo<<endl;
//                cout<<" Please verify the following file : "<<filename<<endl;
//                return false;
//            }
//            else
//                altAlleleString=pch;
//
////            pch = strtok_r (NULL, "\t", &end_str1);
////            if(pch==NULL)
////            {
////                cout<<"\n Info file does NOT have 14 columns at Row : "<<RowNo<<endl;
////                cout<<" Please verify the following file : "<<filename<<endl;
////                return false;
////            }
////            else
////                MajAlleleString=pch;
////
////            pch = strtok_r (NULL, "\t", &end_str1);
////            if(pch==NULL)
////     ifs       {
////                cout<<"\n Info file does NOT have 14 columns at Row : "<<RowNo<<endl;
////                cout<<" Please verify the following file : "<<filename<<endl;
////                return false;
////            }
////            else
////                MinAlleleString=pch;
//
//            tempVariant.assignRefAlt(refAlleleString,altAlleleString);
////            tempVariant.assignMajMin(MajAlleleString,MinAlleleString);
//
////            if(refAlleleString.compare(altAlleleString)==0)
////                tempVariant.assignSwap(true);
////            else
////                tempVariant.assignSwap(false);
//
//            VariantList.push_back(tempVariant);
//            line.clear();
//
//            Indic=ifs->readLine(line);
//
//        }
//    }
//    else
//    {
//        cout<<"\n Following File File Not Available : "<<filename<<endl;
//        return false;
//    }
//
//    numMarkers=VariantList.size();
//
//
//	std::cout << "\n Number of Markers read from Info File   : " << numMarkers << endl;
//	ifclose(ifs);
//	return true;
//
//}
//
