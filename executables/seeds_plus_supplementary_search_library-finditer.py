#!/usr/bin/env python
# coding: utf-8

# In[50]:


import pandas as pd
import regex as re

#mRNA regions dict


dict_mRNARegions = {'5UTR':'./mart_export_5UTRs_2021_oneline.txt',                   '3UTR':'./mart_export_3UTRs_2021_oneline.txt',                   'CDS':'./mart_export_CDS_2021_oneline.txt'}


# In[51]:


def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U':'A'}
    return ''.join([complement[base] for base in dna[::-1]])


# In[52]:


def generarDiccionarioMicroRNAs():
       
   r = pd.read_csv('./miR_Family_Info_targetScan_2016.txt', sep='\t', index_col=3)
   miRseq = pd.DataFrame(r)

   miRseqHuman = miRseq[(miRseq['Species ID'] == 9606)]
   
   dict_regexp = {} 


   for index,row in miRseqHuman.iterrows():
       
       rev = reverse_complement(row['Mature sequence'])
       #seed
       seed6 = row['Mature sequence'][1:7]
       seed7 = row['Mature sequence'][1:8]
       #seed complement
       seed6compl = reverse_complement(seed6)
       seed7compl = reverse_complement(seed7)
       #next nucleotide
       no7 = reverse_complement(row['Mature sequence'][7])
       no8 = reverse_complement(row['Mature sequence'][8])

       #print(row['Mature sequence'][8])

       extra4_13_16 = rev[-16:-12]


       site_6mer = '[^'+ str(no7) + ']' + seed6compl + '[^A]'
       site_7merA1 = '[^'+ str(no7) + ']' + seed6compl + 'A'
       site_7merm8 = '[^'+ str(no8) + ']' + seed7compl + '[^A]'
       site_8mer = '[^'+ str(no8) + ']' + seed7compl + 'A'


       site_6mer_suppl = '(' + str(extra4_13_16) +')' +'((?:(?!'+ str(extra4_13_16) + ').){0,14}?)' +'([^'+ str(no7) + ']' + seed6compl + '[^A]' + ')'
       site_7merA1_suppl = '(' + str(extra4_13_16) +')' +'((?:(?!'+ str(extra4_13_16) + ').){0,14}?)' +'([^'+ str(no7) + ']' + seed6compl + 'A' + ')'
       site_7merm8_suppl = '(' + str(extra4_13_16) +')' +'((?:(?!'+ str(extra4_13_16) + ').){0,14}?)' +'([^'+ str(no8) + ']' + seed7compl + '[^A]' + ')'
       site_8mer_suppl = '(' + str(extra4_13_16) +')' +'((?:(?!'+ str(extra4_13_16) + ').){0,14}?)' +'([^'+ str(no8) + ']' + seed7compl + 'A' + ')'

       dict_regexp[index] = [site_6mer, site_7merA1, site_7merm8, site_8mer, site_6mer_suppl, site_7merA1_suppl, site_7merm8_suppl, site_8mer_suppl]
   
   return dict_regexp
       


# In[53]:


def grabarfile(linea1, linea2, miR, motivo, start, end, wholematch , seqmedia, seqmediaLargo, arch):
                    
    linea1 = linea1.strip()
    linea2 = linea2.strip()
    miR = miR.strip()
    motivo = motivo.strip()
    
    arch.write(linea1)
    arch.write(',')
    arch.write(linea2)
    arch.write(',')
    arch.write(miR)
    arch.write(',')
    arch.write(motivo)
    arch.write(',')
    arch.write(str(start))
    arch.write(',')
    arch.write(str(end))
    arch.write(',')
    arch.write(wholematch)
    arch.write(',')
    arch.write(seqmedia)
    arch.write(',')
    arch.write(str(seqmediaLargo))
    arch.write('\n')   


# In[54]:


def usoREGEXP(microRNAname, mRNAregion):
    #output file name
    nombreDelARchivo = str(microRNAname) + '_finditer_microRNA_6mer_and_7mer_plus_supplementaryNucleotides13_16' + '_mRNA_region_' + str(mRNAregion) + '.csv'
    #open output file
    archivo1 = open(nombreDelARchivo, 'w')
    archivo1.write('GeneName,TranscriptID,TranscriptType,mRNAregion,Interaction,microRNA,Type_of_Interaction,start,end,wholematch,BridgeSequence,BridgeLength\n')
    
    
    #compile regexp
    dict_regexp = generarDiccionarioMicroRNAs() 
    
    prog6mer = re.compile(dict_regexp[microRNAname][0])
    prog7merA1 = re.compile(dict_regexp[microRNAname][1])
    prog7merm8 = re.compile(dict_regexp[microRNAname][2])
    prog8mer = re.compile(dict_regexp[microRNAname][3])
    
    prog6mer_SUPPL = re.compile(dict_regexp[microRNAname][4])
    prog7merA1_SUPPL = re.compile(dict_regexp[microRNAname][5])
    prog7merm8_SUPPL = re.compile(dict_regexp[microRNAname][6])
    prog8mer_SUPPL = re.compile(dict_regexp[microRNAname][7])
    
    
    ciento = open(dict_mRNARegions[mRNAregion], 'r') 
    
    
    line1 = ciento.readline()
    line2 = ciento.readline()
    
    
    
    while line1:
        
        if line1.startswith('>'):
        
            
            line1 = line1.strip()[1:] + '|' + str(mRNAregion)
            line1 = line1.replace('|', ',')                
                

            #-------------------------
            #seeds
            #-------------------------

            for interaction in re.finditer(prog6mer, line2):
                
                start = interaction.start() + 1
                end = interaction.end() - 1
                wholematch= line2[start:end]

                
                grabarfile(line1, 'seed', microRNAname, '6mer_site', start, end, wholematch, '-', 0, archivo1)
                
                

            for interaction in re.finditer(prog7merA1, line2):

                start = interaction.start() + 1
                end = interaction.end()
                wholematch= line2[start:end]

                grabarfile(line1, 'seed' , microRNAname, '7merA1_site', start, end, wholematch,'-', 0, archivo1)
                
                      

                
            for interaction in re.finditer(prog7merm8, line2):

                start = interaction.start() + 1
                end = interaction.end() -1
                wholematch= line2[start:end]

                grabarfile(line1, 'seed' , microRNAname, '7merm8_site', start, end, wholematch, '-', 0, archivo1)




            for interaction in re.finditer(prog8mer, line2):

                start = interaction.start() + 1
                end = interaction.end()
                wholematch= line2[start:end]

                grabarfile(line1, 'seed' , microRNAname, '8mer_site', start, end, wholematch,  '-', 0, archivo1)
                    
                    
                    
                    
            #------------------------
            #Suplementary
            #------------------------

                
            for interaction in re.finditer(prog6mer_SUPPL, line2):              

                start = interaction.start()
                end = interaction.end() - 1
                wholematch= line2[start:end]

                grabarfile(line1, 'suppl' , microRNAname, '6mer_site_SUPPL', start, end, wholematch,  wholematch[4:-6], len(wholematch[4:-6]), archivo1)




            for interaction in re.finditer(prog7merA1_SUPPL, line2):
                
                start = interaction.start()
                end = interaction.end()
                wholematch= line2[start:end]                

                grabarfile(line1, 'suppl' , microRNAname, '7merA1_site_SUPPL', start, end, wholematch, wholematch[4:-7], len(wholematch[4:-7]), archivo1)




            for interaction in re.finditer(prog7merm8_SUPPL, line2):
                start = interaction.start()
                end = interaction.end() - 1
                wholematch= line2[start:end]                

                grabarfile(line1, 'suppl' , microRNAname, '7merm8_site_SUPPL', start, end, wholematch, wholematch[4:-7], len(wholematch[4:-7]), archivo1)




            for interaction in re.finditer(prog8mer_SUPPL, line2):
                
                start = interaction.start()
                end = interaction.end()
                wholematch= line2[start:end]                

                grabarfile(line1, 'suppl' , microRNAname, '8mer_site_SUPPL', start, end, wholematch, wholematch[4:-8], len(wholematch[4:-8]), archivo1)
                    

                

    
        line1 = ciento.readline()
        line2 = ciento.readline()
        
        
        

    archivo1.close()      


# In[55]:


#usoREGEXP('hsa-miR-301b-3p', 'CDS')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[56]:


def diccionariomRNARegions():
    
    diccionario_mRNA_regions_def2 = {}
    
    for region in dict_mRNARegions:
        
            ciento = open(dict_mRNARegions[region], 'r')
    
    
            line1 = ciento.readline()
            line2 = ciento.readline()
            
            
            while line1:

                if line1.startswith('>') and line2.strip() != 'Sequence unavailable':
                    
                    nombre = line1.strip()[1:] + '|' + region
                    
                    diccionario_mRNA_regions_def2[nombre] = line2.strip()
                    
                line1 = ciento.readline()
                line2 = ciento.readline()
                
        
    return diccionario_mRNA_regions_def2


# In[61]:


def buscarDesdeElmRNABlancos(gen_de_interes):
    
    #output file name
    nombreDelARchivo = str(gen_de_interes).replace('|','_') + '_finditer_microRNA_6mer_and_7mer_plus_supplementary_13_16nt' + '.csv'
    #open output file
    archivo1 = open(nombreDelARchivo, 'w')
    archivo1.write('GeneName,TranscriptID,TranscriptType,mRNAregion,Interaction,microRNA,Type_of_Interaction,start,end,wholematch,BridgeSequence,BridgeLength\n')
    
    
    diccionario_mRNA_regions_def2 = diccionariomRNARegions()
    dict_regexp = generarDiccionarioMicroRNAs()
    
    
    mRNAsequence = diccionario_mRNA_regions_def2[gen_de_interes]
    
    for microRNAname in dict_regexp:
        #compile regexp

        prog6mer = re.compile(dict_regexp[microRNAname][0])
        prog7merA1 = re.compile(dict_regexp[microRNAname][1])
        prog7merm8 = re.compile(dict_regexp[microRNAname][2])
        prog8mer = re.compile(dict_regexp[microRNAname][3])
        
        prog6mer_SUPPL = re.compile(dict_regexp[microRNAname][4])
        prog7merA1_SUPPL = re.compile(dict_regexp[microRNAname][5])
        prog7merm8_SUPPL = re.compile(dict_regexp[microRNAname][6])
        prog8mer_SUPPL = re.compile(dict_regexp[microRNAname][7])
        

        
        line1 = gen_de_interes
        line1 = line1.replace('|', ',')


        #-------------------------
        #seeds
        #-------------------------
        for interaction in re.finditer(prog6mer, mRNAsequence):

            start = interaction.start() + 1
            end = interaction.end() - 1
            wholematch= mRNAsequence[start:end]


            grabarfile(line1, 'seed', microRNAname, '6mer_site', start, end, wholematch, '-', 0, archivo1)



        for interaction in re.finditer(prog7merA1, mRNAsequence):

            start = interaction.start() + 1
            end = interaction.end()
            wholematch= mRNAsequence[start:end]

            grabarfile(line1, 'seed' , microRNAname, '7merA1_site', start, end, wholematch,'-', 0, archivo1)




        for interaction in re.finditer(prog7merm8, mRNAsequence):

            start = interaction.start() + 1
            end = interaction.end() -1
            wholematch= mRNAsequence[start:end]

            grabarfile(line1, 'seed' , microRNAname, '7merm8_site', start, end, wholematch, '-', 0, archivo1)




        for interaction in re.finditer(prog8mer, mRNAsequence):

            start = interaction.start() + 1
            end = interaction.end()
            wholematch= mRNAsequence[start:end]

            grabarfile(line1, 'seed' , microRNAname, '8mer_site', start, end, wholematch,  '-', 0, archivo1)




        #------------------------
        #Suplementary
        #------------------------


        for interaction in re.finditer(prog6mer_SUPPL, mRNAsequence):              

            start = interaction.start()
            end = interaction.end() - 1
            wholematch= mRNAsequence[start:end]

            grabarfile(line1, 'suppl' , microRNAname, '6mer_site_SUPPL', start, end, wholematch,  wholematch[4:-6], len(wholematch[4:-6]), archivo1)




        for interaction in re.finditer(prog7merA1_SUPPL, mRNAsequence):

            start = interaction.start()
            end = interaction.end()
            wholematch= mRNAsequence[start:end]                

            grabarfile(line1, 'suppl' , microRNAname, '7merA1_site_SUPPL', start, end, wholematch, wholematch[4:-7], len(wholematch[4:-7]), archivo1)




        for interaction in re.finditer(prog7merm8_SUPPL, mRNAsequence):
            start = interaction.start()
            end = interaction.end() - 1
            wholematch= mRNAsequence[start:end]                

            grabarfile(line1, 'suppl' , microRNAname, '7merm8_site_SUPPL', start, end, wholematch, wholematch[4:-7], len(wholematch[4:-7]), archivo1)




        for interaction in re.finditer(prog8mer_SUPPL, mRNAsequence):

            start = interaction.start()
            end = interaction.end()
            wholematch= mRNAsequence[start:end]                

            grabarfile(line1, 'suppl' , microRNAname, '8mer_site_SUPPL', start, end, wholematch, wholematch[4:-8], len(wholematch[4:-8]), archivo1)
     
        
    archivo1.close()
    
    


# In[62]:


#buscarDesdeElmRNABlancos('AGO2|ENST00000220592.10|protein_coding|3UTR')


# ## def fastatooneline(fastafile):
#     
#     file = open(fastafile, 'r')
#     
#     outputname= fastafile + '_oneline.txt'
#     
#     print(outputname)
# 
#     output= open(outputname, 'w')
#     
#     
#     for i in file.readlines():
#         if i.startswith('>'):
#             output.write('\n')
#             output.write(i)
# 
#         else:
#             output.write(i.strip())
#             
#     
#     

# In[ ]:




