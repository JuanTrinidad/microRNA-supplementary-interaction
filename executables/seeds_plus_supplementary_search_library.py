#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import regex as re

#Genero un diccionario con los nombres de las regiones y los files asociados
###### ACORDARSE DE CAMBIAR LOS ARCHIVOS #########################################
dict_mRNARegions = {'5UTR':'./mart_export_5UTRs_2021_oneline.txt',                   '3UTR':'./mart_export_3UTRs_2021_oneline.txt',                   'CDS':'./mart_export_CDS_2021_oneline.txt'}


# ### Se importa y filtra el archivo para generar el diccionario de expresiones regulares

# In[2]:


def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U':'A'}
    return ''.join([complement[base] for base in dna[::-1]])


# In[3]:


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
       


# Esta def la hice para grabar mi archivo de salida

# In[4]:


def grabarfile(linea1, linea2, miR, motivo, NumMotivo , seqmedia, seqmediaLargo, arch):
                    
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
    arch.write(str(NumMotivo))
    arch.write(',')
    arch.write(seqmedia)
    arch.write(',')
    arch.write(str(seqmediaLargo))
    arch.write('\n')   


# In[15]:


def usoREGEXP(microRNAname, mRNAregion):
    #nombre del archivo de salida
    nombreDelARchivo = str(microRNAname) + '_' + 'microRNA_6mer_and_7mer_plus_supplementaryNucleotides13_16' + '_mRNA_region_' + str(mRNAregion) + '.csv'
    #abro el archivo de salida
    archivo1 = open(nombreDelARchivo, 'w')
    archivo1.write('GeneName,TranscriptID,TranscriptType,mRNAregion,Interaction,microRNA,Type_of_Interaction,Numer_of_Interactions,BridgeSequence,BridgeLength\n')
    
    
    #compilo la expresion regular del microRNA
    dict_regexp = generarDiccionarioMicroRNAs() 
    
    prog6mer = re.compile(dict_regexp[microRNAname][0])
    prog7merA1 = re.compile(dict_regexp[microRNAname][1])
    prog7merm8 = re.compile(dict_regexp[microRNAname][2])
    prog8mer = re.compile(dict_regexp[microRNAname][3])
    
    prog6mer_SUPPL = re.compile(dict_regexp[microRNAname][4])
    prog7merA1_SUPPL = re.compile(dict_regexp[microRNAname][5])
    prog7merm8_SUPPL = re.compile(dict_regexp[microRNAname][6])
    prog8mer_SUPPL = re.compile(dict_regexp[microRNAname][7])
    
    #abro el archivo en el que quiero buscar la expresion regular
    ciento = open(dict_mRNARegions[mRNAregion], 'r') 
    
    
    line1 = ciento.readline()
    line2 = ciento.readline()
    
    
    
    while line1:
        
        if line1.startswith('>'): #== True:
        

            resultado6mer = re.findall(prog6mer, line2)
            resultado7merA1 = re.findall(prog7merA1, line2)
            resultado7merm8 = re.findall(prog7merm8, line2)
            resultado8mer = re.findall(prog8mer, line2)
            
            resultado6mer_SUPPL = re.findall(prog6mer_SUPPL, line2)
            resultado7merA1_SUPPL = re.findall(prog7merA1_SUPPL, line2)
            resultado7merm8_SUPP = re.findall(prog7merm8_SUPPL, line2)
            resultado8mer_SUPPL = re.findall(prog8mer_SUPPL, line2)
            
            line1 = line1.strip()[1:] + '|' + str(mRNAregion)
            line1 = line1.replace('|', ',')                
                

            
            #-------------------------
            #Primero las seeds solas
            #-------------------------
            if len(resultado6mer) != 0:
                
                for regExp in resultado6mer:
                    
                    interaction = ''.join(regExp)

                    grabarfile(line1, 'seed', microRNAname, '6mer_site', len(resultado6mer), '-', 0, archivo1)
                
            if len(resultado7merA1) != 0:
                
                for regExp in resultado7merA1:
                    
                    interaction = ''.join(regExp)

                    grabarfile(line1, 'seed' , microRNAname, '7merA1_site', len(resultado7merA1), '-', 0, archivo1)
                
            if len(resultado7merm8) != 0:
                
                for regExp in resultado7merm8:

                    interaction = ''.join(regExp)

                    grabarfile(line1, 'seed' , microRNAname, '7merm8_site', len(resultado7merm8),  '-', 0, archivo1)
                    
                    
            if len(resultado8mer) != 0:
                
                for regExp in resultado8mer:

                    interaction = ''.join(regExp)

                    grabarfile(line1, 'seed' , microRNAname, '8mer_site', len(resultado8mer),  '-', 0, archivo1)
                    
                    
                    
                    
            #------------------------
            #Ahora las suplementarias
            #------------------------
            if len(resultado6mer_SUPPL) != 0:
                
                for regExp in resultado6mer_SUPPL:                    
                    interaction = ''.join(regExp)

                    grabarfile(line1, interaction[:-1] , microRNAname, '6mer_site_SUPPL', len(resultado6mer_SUPPL), interaction[4:-7], len(interaction[4:-7]), archivo1)

              
            if len(resultado7merA1_SUPPL) != 0:
                
                for regExp in resultado7merA1_SUPPL:
                    interaction = ''.join(regExp)

                    grabarfile(line1, interaction , microRNAname, '7merA1_site_SUPPL', len(resultado7merA1_SUPPL), interaction[4:-7], len(interaction[4:-7]), archivo1)
                    
                
            if len(resultado7merm8_SUPP) != 0:
                
                for regExp in resultado7merm8_SUPP:
                    interaction = ''.join(regExp)

                    grabarfile(line1, interaction[:-1] , microRNAname, '7merm8_site_SUPPL', len(resultado7merm8_SUPP), interaction[4:-8], len(interaction[4:-8]), archivo1)
                    
                    
            if len(resultado8mer_SUPPL) != 0:
                
                for regExp in resultado8mer_SUPPL:
                    interaction = ''.join(regExp)

                    grabarfile(line1, interaction , microRNAname, '8mer_site_SUPPL', len(resultado8mer_SUPPL), interaction[4:-8], len(interaction[4:-8]), archivo1)
                    

                
        line1 = ciento.readline()
        line2 = ciento.readline()
        
        
        

    archivo1.close()
    
   


# In[17]:


#usoREGEXP('hsa-miR-301b-3p', 'CDS')


# In[ ]:





# In[ ]:





# In[12]:


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


# In[13]:


def buscarDesdeElmRNABlancos(gen_de_interes):
    
    #nombre del archivo de salida
    nombreDelARchivo = str(gen_de_interes).replace('|','_') + '_' + '_microRNA_6mer_and_7mer_plus_supplementary_13_16nt' + '.csv'
    #abro el archivo de salida
    archivo1 = open(nombreDelARchivo, 'w')
    archivo1.write('GeneName\tTranscriptID\tTranscriptType\tmRNAregion\tInteraction\tmicroRNA\tType_of_Interaction\tNumer_of_Interactions\tBridgeSequence\tBridgeLength\n')
    
    
    diccionario_mRNA_regions_def2 = diccionariomRNARegions()
    dict_regexp = generarDiccionarioMicroRNAs()
    
    
    mRNAsequence = diccionario_mRNA_regions_def2[gen_de_interes]
    
    for microRNAname in dict_regexp:
        #compilo la expresion regular del microRNA

        prog6mer = re.compile(dict_regexp[microRNAname][0])
        prog7merA1 = re.compile(dict_regexp[microRNAname][1])
        prog7merm8 = re.compile(dict_regexp[microRNAname][2])
        prog8mer = re.compile(dict_regexp[microRNAname][3])
        
        prog6mer_SUPPL = re.compile(dict_regexp[microRNAname][4])
        prog7merA1_SUPPL = re.compile(dict_regexp[microRNAname][5])
        prog7merm8_SUPPL = re.compile(dict_regexp[microRNAname][6])
        prog8mer_SUPPL = re.compile(dict_regexp[microRNAname][7])
        
        #REFINDALL
        
        resultado6mer = re.findall(prog6mer, mRNAsequence)
        resultado7merA1 = re.findall(prog7merA1, mRNAsequence)
        resultado7merm8 = re.findall(prog7merm8, mRNAsequence)
        resultado8mer = re.findall(prog8mer, mRNAsequence)

        resultado6mer_SUPPL = re.findall(prog6mer_SUPPL, mRNAsequence)
        resultado7merA1_SUPPL = re.findall(prog7merA1_SUPPL, mRNAsequence)
        resultado7merm8_SUPP = re.findall(prog7merm8_SUPPL, mRNAsequence)
        resultado8mer_SUPPL = re.findall(prog8mer_SUPPL, mRNAsequence)
        
        line1 = gen_de_interes
        line1 = line1.replace('|', '\t')


        #-------------------------
        #Primero las seeds solas
        #-------------------------
        if len(resultado6mer) != 0:

            for regExp in resultado6mer:

                interaction = ''.join(regExp)

                grabarfile(line1, 'seed', microRNAname, '6mer_site', len(resultado6mer), '-', 0, archivo1)

        if len(resultado7merA1) != 0:

            for regExp in resultado7merA1:

                interaction = ''.join(regExp)

                grabarfile(line1, 'seed' , microRNAname, '7merA1_site', len(resultado7merA1), '-', 0, archivo1)

        if len(resultado7merm8) != 0:

            for regExp in resultado7merm8:

                interaction = ''.join(regExp)

                grabarfile(line1, 'seed' , microRNAname, '7merm8_site', len(resultado7merm8),  '-', 0, archivo1)


        if len(resultado8mer) != 0:

            for regExp in resultado8mer:

                interaction = ''.join(regExp)

                grabarfile(line1, 'seed' , microRNAname, '8mer_site', len(resultado8mer),  '-', 0, archivo1)




        #------------------------
        #Ahora las suplementarias
        #------------------------
        if len(resultado6mer_SUPPL) != 0:

            for regExp in resultado6mer_SUPPL:                    
                interaction = ''.join(regExp)

                grabarfile(line1, interaction[:-1] , microRNAname, '6mer_site_SUPPL', len(resultado6mer_SUPPL), interaction[4:-7], len(interaction[4:-7]), archivo1)


        if len(resultado7merA1_SUPPL) != 0:

            for regExp in resultado7merA1_SUPPL:
                interaction = ''.join(regExp)

                grabarfile(line1, interaction , microRNAname, '7merA1_site_SUPPL', len(resultado7merA1_SUPPL), interaction[4:-7], len(interaction[4:-7]), archivo1)


        if len(resultado7merm8_SUPP) != 0:

            for regExp in resultado7merm8_SUPP:
                interaction = ''.join(regExp)

                grabarfile(line1, interaction[:-1] , microRNAname, '7merm8_site_SUPPL', len(resultado7merm8_SUPP), interaction[4:-8], len(interaction[4:-8]), archivo1)


        if len(resultado8mer_SUPPL) != 0:

            for regExp in resultado8mer_SUPPL:
                interaction = ''.join(regExp)

                grabarfile(line1, interaction , microRNAname, '8mer_site_SUPPL', len(resultado8mer_SUPPL), interaction[4:-8], len(interaction[4:-8]), archivo1)
     
        
    archivo1.close()
    
    


# In[14]:


buscarDesdeElmRNABlancos('AGO2|ENST00000220592.10|protein_coding|3UTR')


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




