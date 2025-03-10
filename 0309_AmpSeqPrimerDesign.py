#!/usr/bin/env python
# coding: utf-8

# In[18]:


import primer3
import pandas
import pybedtools as pbt
import os 
import argparse
from datetime import datetime

headerString = ['Locus']
bedtools_directory = '/opt/homebrew/bin/bedtools'               # Specify the directory containing BEDTools executables
current_path = os.environ.get('PATH', '')                       # Get the current PATH variable
os.environ['PATH'] = f"{current_path}:{bedtools_directory}"     # Append the bedtools_directory to the PATH

now = datetime.today().strftime('%Y%m%d_%H%M%S')

parser = argparse.ArgumentParser(prog = 'Primer3 Design Script', description = 'Submit sites in BED format to generate primer sets using Primer3.')
parser.add_argument('-f', '--fileName', help = 'Input file containing sites of targets in BED format.')
parser.add_argument('-o', '--outDir', default = 'Outputs', help = 'Output directory for primer excel sheet and working files.')
parser.add_argument('-s', '--ampliconSize', default = 200, help = 'Amplicon target size in base pairs.')
parser.add_argument('-b', '--buffer', default = 50, help = 'Acceptable buffer range for amplicon size. buffer-ampliconSize < amplicons < buffer+ampliconSize.')
parser.add_argument('-d', '--designSpace', default = 1.5, help = 'Size of design space for primers as a ratio of amplicon size. Must be >1.0.')
parser.add_argument('-t', '--meltingTemperature', default = 60, help = 'Ideal primer melting temperature, in ˚C.')
parser.add_argument('-n', '--numberSets', default = 3, help = 'Number of sets of primers to design.')
parser.add_argument('--FWDHandle', default = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT', help = 'Primer handle sequence for FWD primers.')
parser.add_argument('--REVHandle', default = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT', help = 'Primer handle sequence for REV primers.')
parser.add_argument('-p', '--probe01', default = 0, help = 'Binary input for whether or not to design a probe. Yes (1), No (0).')
parser.add_argument('-r', '--referencePath', default = '/Users/labuser/References/hg38.fa', help = 'Path to reference genome file used to fetch design space sequences.')
args = parser.parse_args()



# In[29]:
print('\nScript to design primers using target sites in BED file format. Relies on Primer3 \nto do the primer designs. Use the -h parameter for list of input parameters.\n\n')

isAmpSeq = True
ampSeqBuffer = int(int(args.ampliconSize)*float(args.designSpace)/2)
primersOptTM = int(args.meltingTemperature)
numSets = int(args.numberSets)
headerString = ['Name']

os.makedirs(str(args.outDir), exist_ok = True)
inputbed = pandas.read_excel(args.fileName, header = 0)

inputbed.reset_index(drop = True, inplace = True)

primerSeqs = pandas.DataFrame()


# In[20]:


print('----------DESIGNING PRIMERS FOR THESE SITES----------\n')
print(inputbed)


# In[21]:


if isAmpSeq:
    for a in range(len(inputbed)):
        inputbed.loc[a,'start'] = int(inputbed.loc[a,'start'])-ampSeqBuffer
        inputbed.loc[a,'stop'] = int(inputbed.loc[a,'stop'])+ampSeqBuffer
    inputbed.to_csv(os.path.join(args.outDir, 'DesignSpaceBED_'+now+'.txt'), sep = '\t', header = None, index = None)
    
    pbtObj = pbt.BedTool(os.path.join(args.outDir, 'DesignSpaceBED_'+now+'.txt'))                    
    seq= pbtObj.getfasta(fi= args.referencePath, bed = os.path.join(args.outDir, 'DesignSpaceBED_'+now+'.txt'), fo = os.path.join(args.outDir, 'DesignSequences_'+now+'.txt')) 
fastaSeqs = pandas.read_csv(os.path.join(args.outDir, 'DesignSequences_'+now+'.txt'), header = None, sep = '\t')

seqNames = []
dropme = []
for a in range(0, len(fastaSeqs), 2):
    seqNames.extend(fastaSeqs.loc[a, 0])
    dropme.extend([a])
    
fastaSeqs = fastaSeqs.drop(dropme)
fastaSeqs = fastaSeqs.reset_index(drop = True)
fastaSeqs


# In[71]:


print('\n\n----------BEGINNING PRIMER DESIGN----------\n')
runningList = []
for a in range(0, len(fastaSeqs), 1):
    print('Designing for: '+str(inputbed.loc[a, 'name']))
    itercount = 0
    
    primers = primer3.bindings.design_primers(
        {'SEQUENCE_ID': inputbed.loc[a, 'name'], #Sequence name from BED file Input
          'SEQUENCE_TEMPLATE':fastaSeqs.values[a][0], # Sequence of design space
          },
        {
          'PRIMER_PICK_LEFT_PRIMER':1,
          'PRIMER_PICK_RIGHT_PRIMER':1,
          'PRIMER_PICK_INTERNAL_OLIGO':args.probe01,
          'PRIMER_NUM_RETURN': numSets, 
          'PRIMER_OPT_SIZE': 20,
          'PRIMER_MIN_SIZE': 17,
          'PRIMER_MAX_SIZE': 25,
          'PRIMER_OPT_TM': primersOptTM,
          'PRIMER_MIN_TM': primersOptTM-3,
          'PRIMER_MAX_TM': primersOptTM+3,
          'PRIMER_INTERNAL_OPT_TM':primersOptTM+5,
          'PRIMER_INTERNAL_MAX_TM':(primersOptTM+5)+3,
          'PRIMER_INTERNAL_MIN_TM':(primersOptTM+5)-3,  
          'PRIMER_MIN_GC': 20.0,
          'PRIMER_MAX_GC': 80.0,
          'PRIMER_MAX_POLY_X': 5,
          'PRIMER_SALT_MONOVALENT': 50.0,
          'PRIMER_DNA_CONC': 50.0,
          'PRIMER_MAX_NS_ACCEPTED': 0,
          'PRIMER_MAX_SELF_ANY': 12,
          'PRIMER_MAX_SELF_END': 8,
          'PRIMER_PAIR_MAX_COMPL_ANY': 12,
          'PRIMER_PAIR_MAX_COMPL_END': 8,
          'PRIMER_PRODUCT_SIZE_RANGE': [int(args.ampliconSize)-int(args.buffer), int(args.ampliconSize)+int(args.buffer)]})

    primerSeqs.at[a, 0+itercount] = inputbed.iloc[a, 3]
    
    for b in range(0, int(numSets),1):
        checksArray = []
      
        checksArray.extend([primers['PRIMER_RIGHT_NUM_RETURNED'], primers['PRIMER_LEFT_NUM_RETURNED'], primers['PRIMER_INTERNAL_NUM_RETURNED']])
        rightSeq = 'No Right Primer Designed'
        rightTm = 'N/A'
        leftSeq = 'No Left Primer Designed'
        leftTm = 'N/A'
        intSeq = 'No Probe Designed'
        intTm = 'N/A'
        size = 'N/A'
        penalty = 'N/A'
        
        if checksArray[0] !=0:
            rightSeq = primers['PRIMER_RIGHT_'+str(b)+'_SEQUENCE']
            rightTm = primers['PRIMER_RIGHT_'+str(b)+'_TM']
            
        if checksArray[1]!=0:
            leftSeq = primers['PRIMER_LEFT_'+str(b)+'_SEQUENCE']
            leftTm = primers['PRIMER_LEFT_'+str(b)+'_TM']
            
        if checksArray[2]!=0:
            intSeq = primers['PRIMER_INTERNAL_'+str(b)+'_SEQUENCE']
            intTm = primers['PRIMER_INTERNAL_'+str(b)+'_TM']
            
        size = primers['PRIMER_PAIR_'+str(b)+'_PRODUCT_SIZE']
        penalty = primers['PRIMER_PAIR_'+str(b)+'_PENALTY']
        print('\tSet '+str(b+1)+':'+'   FWD - ' + str(leftSeq).upper() + '  ' + 'REV - ' + str(rightSeq).upper())
            
        runningList.append([inputbed.loc[a, 'name'], 'Set '+str(b+1), leftSeq, leftTm, rightSeq, rightTm, intSeq, intTm, size, penalty])
print('\n----------PRIMER DESIGN COMPLETE----------')
        
primerDF = pandas.DataFrame(runningList, columns = ['Site', 'Set Number', 'Left Primer', 'Left Primer Tm (C)', 'Right Primer', 'Right Primer Tm (C)', 'Probe', 'Probe Tm (C)', 'Product Size', 'Product Penalty Score'])


# In[67]:


print('\n----------APPENDING HANDLES----------')

primerDF_Adapted = primerDF.copy()
for a in range(len(primerDF_Adapted)):
    primerDF_Adapted.loc[a, 'Left Primer'] = str(args.FWDHandle) + str(primerDF_Adapted.loc[a, 'Left Primer'])
    primerDF_Adapted.loc[a, 'Right Primer'] = str(args.REVHandle) + str(primerDF_Adapted.loc[a, 'Right Primer'])


# In[70]:


with pandas.ExcelWriter(os.path.join(args.outDir, 'DesignedPrimers'+now+'.xlsx')) as writer:
    primerDF.set_index(['Site', 'Set Number']).to_excel(writer, sheet_name = 'GSP Primers', index = True)
    primerDF_Adapted.set_index(['Site', 'Set Number']).to_excel(writer, sheet_name = 'Adapted Primers', index = True)

print('\n----------COMPLETE----------\n\n\n')

