

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import os
import argparse
import re
import gffutils
import pysam
import gzip
import pyranges as pr
from natsort import natsorted, index_natsorted, order_by_index
import sys

argparser = argparse.ArgumentParser(description = 'This Software is a part of a pipeline tool to design Guiding RNA. \n This script relates back the CRISPOR output an initial genome annotation and optionnaly can provide editor specific annotation (Clinvar and/or VCF file)')
argparser.add_argument('-E','--Editor', metavar = 'file name', dest = 'Editor', nargs='+', type = str, required = False, default=None, help = 'Editor type (will be examined against reference)')
argparser.add_argument('-S','--ScoreGuide', metavar = 'file name', dest = 'scoreGuide', type = str, required = True, help = 'Crispor ouput')
argparser.add_argument('-O','--Output', metavar = 'file', dest = 'Output', type = str, required = False,default='out', help = 'prefix of the vcf')
argparser.add_argument('-X','--exclude', metavar = 'file', dest = 'excludes', type = str, required = False, default ='', help = 'List of Guides to exclude')
argparser.add_argument('-V','--Per_Variant', dest = 'per_variant', action='store_true', help = 'flag to produce a per variant VCF suitable for VEP')
argparser.add_argument('-R','--Per_sgRNA', dest = 'per_guide', action='store_true', help = 'flag to produce a per guide vcf file suitable for VEP')
argparser.add_argument('-L','--length', metavar = 'int', dest = 'length', type = int, required = False, default ='20', help = 'length of the GuideRNA without PAM')
argparser.add_argument('-B','--bed', metavar = 'file', dest = 'bed', type = str, required = True, help = 'bedFile protein per region')
argparser.add_argument('--gc', dest = 'gc', action='store_true', required = False, help = 'flag not Consider C in GC as affected')
argparser.add_argument('-G','--Genome', metavar = 'file', dest = 'Genome', type = str, required = '--gc' in sys.argv, help = 'Genome fasta')


def MutateWindow(rower):
        if args.gc and ((rower.strand=='+' and i.BE[0].upper()=='C') or (rower.strand=='-' and i.BE_RC[0].upper()=='G')) :
                returned=""
                for nuc in range(rower.editing_windowSTART-1,rower.editing_windowEND):
                        if rower.strand=='+' and str(Genome_dict[rower.Chromosome][nuc]).upper() ==i.BE[0].upper() and Genome_dict[rower.Chromosome][nuc-1].upper()!='G' :
                                returned=returned+i.BE[1]
                        elif rower.strand=='-' and str(Genome_dict[rower.Chromosome][nuc]).upper() ==i.BE_RC[0].upper() and Genome_dict[rower.Chromosome][nuc+1].upper()!='C' :
                                returned=returned+i.BE_RC[1]
                        else:
                                returned=returned+Genome_dict[rower.Chromosome][nuc]
        else :
                returned=rower.editing_windowSeq.replace(i.BE[0],i.BE[1]) if rower.strand=='+' else rower.editing_windowSeq.replace(i.BE_RC[0],i.BE_RC[1])
        return returned


if __name__ == '__main__':
        args = argparser.parse_args()
        DataDir=os.path.abspath(os.path.dirname(__file__))
        if not os.path.exists(DataDir+'/data/Editor_library.csv') :
                raise Exception('Necessary auxiliairy files are absent.')
        else :
                if not args.Editor == None :
                        data=pd.read_csv(DataDir+'/data/Editor_library.csv',index_col='Editor',sep='\t')
                        editor=data.loc[args.Editor]
                        del data
        if args.gc :
                Genome_dict = SeqIO.to_dict(SeqIO.parse(DataDir+'/project/rrg-vmooser/vct/transfert/CRISPOR/'+args.Genome+".fa", "fasta"))
        scoreGuides=pd.read_csv(args.scoreGuide,sep='\t')
        if not args.excludes == "" :
                excludes=pd.read_csv(args.excludes,header=None,sep='\t')
                scoreGuides=scoreGuides[[not(i in excludes[0].tolist()) for i in scoreGuides.targetSeq]]
        scoreGuides.drop_duplicates(subset='targetSeq', inplace=True, keep=False)
        scoreGuides['targetSeq_plusStrand']=[str(Seq(row.targetSeq).reverse_complement()) if 'rev' in  row['guideId']  else row['targetSeq'] for index, row in scoreGuides.iterrows()]
        scoreGuides['protospacer']=[j[0:args.length].upper() for j in scoreGuides.targetSeq]
        scoreGuides['PAM']=[j[args.length:len(j)].upper() for j in scoreGuides.targetSeq]
        ### Find protein corresponding
        scoreGuides['positions']=[int(i.replace('rev','').replace('forw','')) for i in scoreGuides.guideId]
        scoreGuides['Chromosome']=[re.split('-|:', j)[0] for j in scoreGuides['#seqId']]
        scoreGuides['start_seqId']=[int(re.split('-|:', j)[1]) for j in scoreGuides['#seqId']]
        scoreGuides['target_Start']=scoreGuides['start_seqId']+scoreGuides['positions']
        ### Find proiten
        bed=pr.read_bed(args.bed)
        protein=[]
        for index, row in scoreGuides.iterrows():
                df=pd.DataFrame({"Chromosome": [row['Chromosome']], "Start": [row.target_Start-args.length if 'rev' in  row.guideId  else row.target_Start ],"End":[row.target_Start if 'rev' in  row.guideId  else row.target_Start + args.length]})
                position= pr.PyRanges(df)
                names=[i for i in position.join(bed).Name]
                if len(set(names))>1:
                        raise Exception('One protein overlaps with anoter.')
                else :
                        protein.extend(list(set(names)))
        scoreGuides['Protein']=protein
        prots=set(protein)
        scoreGuides['strand']=['-'if 'rev' in j  else '+' for j in scoreGuides.guideId]
        scoreGuides.sort_values(['Protein','guideId','strand'])
        scoreGuides['num']= scoreGuides.groupby(['Protein']).cumcount()+1
        scoreGuides['ID']=scoreGuides['Protein']+'_'+ [str(j) for j in scoreGuides['num']]
        scoreGuides.index=scoreGuides['ID']
        with  open(args.Output+'_general.csv', 'wt') as general :
                general.write('ID,Protospacer,PAM,gRNA_seq_POSstrand,Chromosome,POSstart,POSend,strand,cfdSpecScore,Doench\'16-Score\n')
                for index, row in scoreGuides.iterrows():
                        general.write('{ID},{protospacer},{PAM},{gRNA_seq_POS},{chrom},{POSstart},{Strand},{cfdSpecScore},{Doench}\n'.format(
                                ID=row.ID,
                                protospacer=row['protospacer'],
                                PAM=row.PAM,
                                gRNA_seq_POS=str(Seq(row.targetSeq).reverse_complement()) if 'rev' in  row['guideId']  else row['targetSeq'],
                                chrom=row['Chromosome'],
                                Strand=row.strand,
                                POSstart=str(row.target_Start),
                                cfdSpecScore=str(row.cfdSpecScore),
                                Doench=str(row['Doench \'16-Score'])))
        ### producting Edditor specific files
        if not args.Editor == None :
                for index, i in editor.iterrows() :
                        editor_df=scoreGuides.copy()
                        i.window_start=i.window_start
                        i.window_end=i.window_end
                        editor_df['editing_windowSeq']=[j[i.window_start-1:i.window_end].upper() for j in editor_df.targetSeq]
                        editor_df['editing_windowSeq']= [str(Seq(editor_df['editing_windowSeq'][j]).reverse_complement()) if 'rev' in  editor_df['guideId'][j]  else editor_df['editing_windowSeq'][j]   for j in editor_df.index]
                        editor_df['is_empty'] = [not (editor.BE[0] in x.upper() or editor.BE_RC[0] in x.upper()) for x in editor_df["editing_windowSeq"]]
                        editor_df['editing_windowSTART']=[editor_df['target_Start'][j]-args.length+i.window_start-1 if 'forw' in editor_df.guideId[j] else editor_df['target_Start'][j] +len(editor_df.targetSeq[j])-i.window_end for j in editor_df.index]
                        editor_df['editing_windowEND']=editor_df['editing_windowSTART'] + (i.window_end-i.window_start)
                        editor_df['editing_window_mutated']=[MutateWindow(rowed) for index, rowed in editor_df.iterrows()]
                        editor_df['nchange']=[sum([not rowed.editing_windowSeq[j]==rowed.editing_window_mutated[j] for j in range(0,len(rowed.editing_window_mutated))]) for indexed, rowed in editor_df.iterrows()]
                        editor_df=editor_df.reindex(index=order_by_index(editor_df.index, index_natsorted(zip(editor_df.Chromosome, editor_df.target_Start))))
                        with open(args.Output+'.'+i.name +'_empties' + '.txt', 'w') as empties:
                                for ind, row in editor_df.iterrows():
                                        if row.editing_windowSeq ==row.editing_window_mutated :
                                                empties.write(str(row.name) +'\n')
                        if args.per_guide :
                                with gzip.open(args.Output+'.'+i.name+'_per_guides.vcf.gz','wt') as vcf :
                                        vcf.write('##fileformat=VCFv4.1\n')
                                        vcf.write('##contig=<ID=chr1,length=249250621>\n##contig=<ID=chr2,length=243199373>\n##contig=<ID=chr3,length=198022430>\n##contig=<ID=chr4,length=191154276>\n##contig=<ID=chr5,length=180915260>\n##contig=<ID=chr6,length=171115067>\n##contig=<ID=chr7,length=159138663>\n##contig=<ID=chr8,length=146364022>\n##contig=<ID=chr9,length=141213431>\n##contig=<ID=chr10,length=135534747>\n##contig=<ID=chr11,length=135006516>\n##contig=<ID=chr12,length=133851895>\n##contig=<ID=chr13,length=115169878>\n##contig=<ID=chr14,length=107349540>\n##contig=<ID=chr15,length=102531392>\n##contig=<ID=chr16,length=90354753>\n##contig=<ID=chr17,length=81195210>\n##contig=<ID=chr18,length=78077248>\n##contig=<ID=chr19,length=59128983>\n##contig=<ID=chr20,length=63025520>\n##contig=<ID=chr21,length=48129895>\n##contig=<ID=chr22,length=51304566>\n##contig=<ID=chrM,length=16571>\n##contig=<ID=chrX,length=155270560>\n##contig=<ID=chrY,length=59373566>\n')
                                        vcf.write('##INFO=<ID=Protospacer,Number=.,Type=String,Description="Protospace"> \n')
                                        vcf.write('##INFO=<ID=STRAND,Number=.,Type=String,Description="Strand"> \n')
                                        vcf.write('##INFO=<ID=PAM,Number=.,Type=String,Description="Protospacer adjacent motif"> \n')
                                        vcf.write('##INFO=<ID=Nchange,Number=.,Type=Integer,Description="Nummber of nucleotide changes"> \n')
                                        vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
                                        for ind, row in editor_df.iterrows():
                                                vcf.write('{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{END}\t{FILTER}\t{INFO}\n'.format(
                                                        chrom = str(row['Chromosome']),
                                                        pos = str(row['editing_windowSTART']) ,
                                                        vid=str(row.ID),
                                                        ref= str(row['editing_windowSeq']) ,
                                                        alt= str(row['editing_window_mutated']) ,
                                                        END='.',
                                                        FILTER ='.' ,
                                                        INFO='Protospacer='+row['protospacer']+ ';PAM='+ row['PAM'] +';STRAND='+ row.strand + ';Nchange=' + str(row.nchange)))
                        if args.per_variant :
                                ranges=editor_df.copy()[['Chromosome','editing_windowSTART','editing_windowEND','strand']]
                                ranges.rename(columns={'Chromosome':"Chromosome","editing_windowSTART":'Start', 'editing_windowEND':'End'}, inplace=True)
                                ranges['Start']=ranges.Start
                                ranges['End']= ranges.End +1      # To include the ending
                                ranges.index=range(0,len(ranges))
                                ranges['Sequence']=[j for j in editor_df.protospacer]
                                ranges['PAM']=[j for j in editor_df.PAM]
                                ranges['strand']=[j for j in editor_df.strand]
                                ranges['ID']=[j for j in editor_df.ID]
                                unmerged=pr.PyRanges(ranges)
                                with gzip.open(args.Output+'.'+i.name+'_per_variant.vcf.gz','wt') as vcf :
                                        vcf.write('##fileformat=VCFv4.1\n')
                                        vcf.write('##contig=<ID=chr1,length=249250621>\n##contig=<ID=chr2,length=243199373>\n##contig=<ID=chr3,length=198022430>\n##contig=<ID=chr4,length=191154276>\n##contig=<ID=chr5,length=180915260>\n##contig=<ID=chr6,length=171115067>\n##contig=<ID=chr7,length=159138663>\n##contig=<ID=chr8,length=146364022>\n##contig=<ID=chr9,length=141213431>\n##contig=<ID=chr10,length=135534747>\n##contig=<ID=chr11,length=135006516>\n##contig=<ID=chr12,length=133851895>\n##contig=<ID=chr13,length=115169878>\n##contig=<ID=chr14,length=107349540>\n##contig=<ID=chr15,length=102531392>\n##contig=<ID=chr16,length=90354753>\n##contig=<ID=chr17,length=81195210>\n##contig=<ID=chr18,length=78077248>\n##contig=<ID=chr19,length=59128983>\n##contig=<ID=chr20,length=63025520>\n##contig=<ID=chr21,length=48129895>\n##contig=<ID=chr22,length=51304566>\n##contig=<ID=chrM,length=16571>\n##contig=<ID=chrX,length=155270560>\n##contig=<ID=chrY,length=59373566> \n')
                                        vcf.write('##INFO=<ID=Protospacer,Number=.,Type=String,Description="Protospace"> \n')
                                        vcf.write('##INFO=<ID=STRAND,Number=.,Type=String,Description="Strand"> \n')
                                        vcf.write('##INFO=<ID=PAM,Number=.,Type=String,Description="Protospacer adjacent motif"> \n')
                                        if args.gc:
                                                vcf.write('##INFO=<ID=GC_FLAG,Number=.,Type=bool,Description="Mutation effectiveness compromised by GC pattern"> \n')
                                        vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
                                        output=set()
                                        for ind, row in editor_df.iterrows():
                                                for nuc in range(0,len(row.editing_windowSeq)):
                                                        if row.editing_windowSeq[nuc]==i.BE[0] and row.strand=='+':
                                                                output.add((row.Chromosome, row.editing_windowSTART+nuc, i.BE[0], i.BE[1],row.strand))
                                                        if row.editing_windowSeq[nuc]==i.BE_RC[0] and row.strand=='-':
                                                                output.add((row.Chromosome, row.editing_windowSTART+nuc, i.BE_RC[0], i.BE_RC[1],row.strand))
                                        df = pd.DataFrame(list(output),columns =['chrom', 'pos', 'REF','ALT','strand'],index=range(0,len(output)))
                                        output_df=df.reindex(index=order_by_index(df.index, index_natsorted(zip(df.chrom, df.pos))))
                                        for ind, row in output_df.iterrows():
                                                df=pd.DataFrame({"Chromosome": [str(row.chrom)], "Start": [row.pos],"End":[row.pos+1]})
                                                position= pr.PyRanges(df)
                                                overlap=position.join(unmerged).as_df()
                                                overlap=overlap.loc[overlap['strand']==row['strand']]
                                                if args.gc :
                                                        if row['strand']=='+' and str(Genome_dict[row.chrom][row.pos-1]).upper()=='C' and str(Genome_dict[row.chrom][row.pos-2]).upper()=='G':
                                                               gc_flag=';GC_FLAG=True'
                                                        elif  row['strand']=='-' and str(Genome_dict[row.chrom][row.pos-1]).upper()=='G' and str(Genome_dict[row.chrom][row.pos]).upper()=='C':
                                                                gc_flag=';GC_FLAG=True'
                                                        else :
                                                                gc_flag=';GC_FLAG=False'
                                                else :
                                                        gc_flag=''
                                                if not overlap.empty :
                                                        vcf.write('{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{END}\t{FILTER}\t{INFO}\n'.format(
                                                                chrom = str(row['chrom']),
                                                                pos = str(row['pos']) ,
                                                                vid=str(row.chrom)+'_'+str(row['pos'])+"_"+str(row['REF'])+"_"+str(row['ALT']),
                                                                ref= str(row['REF']) ,
                                                                alt= str(row['ALT']) ,
                                                                END='.',
                                                                FILTER ='.' ,
                                                                INFO='STRAND='+ row.strand+';GuideId='+'|'.join(overlap.ID)+';PAM='+ '|'.join(overlap['PAM'])+';protospacers='+'|'.join(overlap.Sequence)+gc_flag))
