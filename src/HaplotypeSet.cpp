#include "HaplotypeSet.h"


bool HaplotypeSet::WriteMachFile()
{

	VcfFileReader inFile;
	VcfHeader header;
	VcfRecord record;
    individualName.clear();

    String filename=VcfDose;
    int numSamplesRead = 0,markerCount=0,totmarkerCount=0,NumGP=0,NumDS=0,NumGT=0;
    if (!inFile.open(filename, header))
	{
		cout << "\n Program could NOT open file : " << filename << endl<<endl;
		return false;
	}

    inFile.setSiteOnly(false);
	numSamplesRead = header.getNumSamples();
    numSamples=numSamplesRead;
    for (int i = 0; i < numSamplesRead; i++)
	{
		string tempName(header.getSampleName(i));
		individualName.push_back(tempName);
    }

    IFILE machPartial=NULL,plinkMain=NULL,plinkMap=NULL;
    int part=1;

    if(Type=="mach")
        machPartial = ifopen(outFile + ".output.part.0" + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
    else
    {
        plinkMain = ifopen(outFile + ".plink.fam", "wb",InputFile::UNCOMPRESSED);
        for (int i = 0; i < numSamplesRead; i++)
        {
             ifprintf(plinkMain,"%s\t%s\t0\t0\t0\t-9\n",individualName[i].c_str(),individualName[i].c_str());
        }
        ifclose(plinkMain);
        plinkMap = ifopen(outFile + ".plink.map", "wb",InputFile::UNCOMPRESSED);
        plinkMain = ifopen(outFile + ".plink."+(Format=="DS" ? "dose" : "gprob") + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);

        ifprintf(plinkMain,"SNP\tA1\tA2");
        for (int i = 0; i < numSamplesRead; i++)
        {
             ifprintf(plinkMain,"\t%s\t%s",individualName[i].c_str(),individualName[i].c_str());
        }
        ifprintf(plinkMain,"\n");
    }

    char * pch_split,* pch_split3;
    char * pch;
    char *end_str1,*end_str3;

    if(Type=="mach")
        for (int i = 0; i < numSamplesRead; i++)
        {
             ifprintf(machPartial,"%s\tDOSE\n",individualName[i].c_str());
        }

    if(Type=="mach")
        ifclose(machPartial);

//    VariantList.clear();
    string name,refAlleleString,altAlleleString;
    int bp;
    string chr;

    int bufferSize=100000;

    dosage.resize(bufferSize);
    for (int i = 0; i < bufferSize; i++)
    {
         dosage[i].resize(numSamplesRead,-9.0);
    }


    GP1.resize(bufferSize);
    GP2.resize(bufferSize);
    for (int i = 0; i < bufferSize; i++)
    {
        GP1[i].resize(numSamplesRead,-9.0);
        GP2[i].resize(numSamplesRead,-9.0);
    }



    IFILE DoseRead = ifopen(filename, "r");
    string line;
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

            printf("    Reading Chunk %d of %d markers ... \n", part, bufferSize);

            while(markerCount<bufferSize && line.compare("")!=0)
            {
                if(totmarkerCount>=numMarkers)
                {
                    std::cout << "\n ERROR !!! Number of markers found in info file = "<<numMarkers ;
                    std::cout << "\n           More than "<<numMarkers <<" markers are found in dosage VCF file ... " ;
                    std::cout << "\n           Please use info file from same imputation run !!! " << endl;
                    return false;
                }


                char *end_str9;
                pch = strtok_r ((char*)line.c_str(),"\t", &end_str1);
                chr=pch;
                pch = strtok_r (NULL, "\t", &end_str1);

                bp=atoi(pch);
                pch = strtok_r (NULL, "\t", &end_str1);
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

                                    GP1[markerCount][indCount]=atof(pch_split3);
                                    pch_split3 = strtok_r (NULL,",", &end_str3);
                                    GP2[markerCount][indCount]=atof(pch_split3);
                                    dosage[markerCount][indCount]=GP2[markerCount][indCount]+(2*(1-GP1[markerCount][indCount]-GP2[markerCount][indCount]));
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
                                    std::cout << "\n ERROR !!! No Genotype or Genotype Probability found for variant "<<VariantList[totmarkerCount].name ;
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
                                            std::cout << "\n ERROR !!! Invalid value of GT ["<<(string)pch_split<<"] found for variant "<<VariantList[totmarkerCount].name ;
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
                                            std::cout << "\n ERROR !!! Invalid value of GT ["<<(string)pch_split<<"] found for variant "<<VariantList[totmarkerCount].name ;
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

                    if(markerCount==bufferSize-1)
                    {
                        if(Format=="DS")
                            printf("    Writing Chunk %d of %d markers to PLINK Dosage Format ... \n", part, bufferSize);
                        else
                            printf("    Writing Chunk %d of %d markers to PLINK GP Format ... \n", part, bufferSize);
                    }

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

            if(Type=="mach")
            {

                    if(Format=="DS")
                        printf("    Writing Chunk %d of %d markers to MaCH Dosage Format ... \n", part, bufferSize);
                    else
                        printf("    Writing Chunk %d of %d markers to MaCH GP Format ... \n", part, bufferSize);
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
                }

            part++;



        }


    }
    else
    {
        cout<<"\n Following File File Not Available : "<<filename<<endl;
        return false;
    }






    std::cout << "\n\n Number of Samples read from VCF Dosage File              : " << numSamplesRead ;
	std::cout << "\n Number of Markers read from VCF Dosage File              : " << totmarkerCount ;

    if(totmarkerCount<numMarkers)
    {
        std::cout << "\n\n ERROR !!! Number of markers found in info file = "<<numMarkers ;
        std::cout << "\n           Less than "<<numMarkers <<" markers were read from dosage VCF file ... " ;
        std::cout << "\n           Please use info file from same imputation run !!! " << endl;
        return false;
    }
	std::cout << "\n Number of Variants imported from GP (Genotype Prob)      : " << NumGP/numSamplesRead;
	std::cout << "\n Number of Variants imported from DS (Dosage)             : " << NumDS/numSamplesRead;
	std::cout << "\n Number of Variants imported from GT (Genotype)           : " << NumGT/numSamplesRead<<endl<<endl;


    if((NumGP+NumDS)==0)
    {
        if(Format=="DS")
        {
            std::cout << "\n WARNING !!! All Dosage values imported from GT values " ;
            std::cout << "\n             No GP or DS values were found in VCF file ! " ;
        }
        else
        {
            std::cout << "\n WARNING !!! All Genotype Probabilities calculated from GT values " ;
            std::cout << "\n             No GP values were found in VCF file ! " ;
            std::cout << "\n             GT values (although may be available,  " ;
            std::cout << "\n             cannot be used for GP calculation ! " ;
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
            if(j%1000==0)
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


        cout<<" \n [NOTE: In output MaCH files, dosage values are twice the probability of";
        printf("\n        ALTERNATE ALLELE (in the info file) and NOT twice the probability");
        cout<< "\n        of MAJOR ALLELE (as is in usual Minimac output)"<<endl;
        cout<<   "        Similarly, Genotype probabilities are given for REF|REF followed"<<endl;
        cout<<   "        REF|ALT (allele given in info file)."<<endl<<endl;



        IFILE ifs = ifopen(Info, "r");
        machFinal = ifopen(outFile + ".mach.info",  "wb", InputFile::UNCOMPRESSED);

        line.clear();
        while ((ifs->readLine(line))!=-1)
        {
            ifprintf(machFinal,"%s\n",line.c_str());
            line.clear();
        }
        ifclose(machFinal);


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


bool HaplotypeSet::LoadInfoFile(String filename)
{

    IFILE ifs = ifopen(filename, "r");

    VariantList.clear();
    string name,refAlleleString,altAlleleString;
    string MajAlleleString,MinAlleleString;
    int bp=100;
    string chr="Z";
    char *pch,*end_str1;





    string line;
    if(ifs)
    {
        ifs->readLine(line);
        line.clear();
        while ((ifs->readLine(line))!=-1)
        {

            pch = strtok_r ((char*)line.c_str(),"\t", &end_str1);
            name=pch;

            variant tempVariant(name,chr,bp);


            pch = strtok_r (NULL, "\t", &end_str1);
            refAlleleString=pch;
            pch = strtok_r (NULL, "\t", &end_str1);
            altAlleleString=pch;
            pch = strtok_r (NULL, "\t", &end_str1);
            MajAlleleString=pch;
            pch = strtok_r (NULL, "\t", &end_str1);
            MinAlleleString=pch;

            tempVariant.assignRefAlt(refAlleleString,altAlleleString);
            tempVariant.assignMajMin(MajAlleleString,MinAlleleString);

            if(refAlleleString.compare(altAlleleString)==0)
                tempVariant.assignSwap(true);
            else
                tempVariant.assignSwap(false);

            VariantList.push_back(tempVariant);
            line.clear();

        }
    }
    else
    {
        cout<<"\n Following File File Not Available : "<<filename<<endl;
        return false;
    }

    numMarkers=VariantList.size();


	std::cout << "\n Number of Markers read from Info File   : " << numMarkers << endl;
	ifclose(ifs);
	return true;

}





bool HaplotypeSet::FastLoadHaplotypes()
{
    String filename=VcfDose;



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


    cout<<"\n Reading Info File                       : "<<Info<<endl;


    if(!LoadInfoFile(Info))
        return false;



    std::cout << "\n Loading Dosage Data from VCF File       : " << filename << endl<<endl;

    return WriteMachFile();

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



