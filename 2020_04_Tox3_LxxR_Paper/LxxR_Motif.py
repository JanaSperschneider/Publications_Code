import os 
import sys
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from subprocess import Popen, PIPE
# -----------------------------------------------------------------------------------------------------------
def signalp_results(CLEAVAGE_SITES):

    f = open('signalp3.out', 'r')
    content = f.readlines()
    f.close()

    for line in content: 
        if line.startswith('#'):
            pass
        else:
            result = line.split(' ')
            result = [x for x in result if x != '']
            protein = result[0]
            cleavage = result[5]
            d_score = result[13].split('\t')[0]
            if d_score == 'Y':
                CLEAVAGE_SITES[protein] = int(cleavage)
            else:
                CLEAVAGE_SITES[protein] = 0

    return CLEAVAGE_SITES
# -----------------------------------------------------------------------------------------------------------
def get_seqs_ids_fasta(FASTA_FILE):
    """ Function: get_seqs_ids_fasta()
        Purpose:  Given a FASTA format file, this function extracts
                  the list of identifiers and the list of sequences 
                  in the order in which they appear in the FASTA file.
              
        Input:    Path to FASTA format file.
    
        Return:   List of identifiers and list of sequences in the order 
                  in which they appear in the FASTA file.
    """ 
    identifiers = []
    sequences = []

    with open(FASTA_FILE) as f: 
        content = f.readlines()

        for position, line in enumerate(content):
            if '>' in line:
                identifiers.append(line)
                seq = []
                following_lines = content[position + 1:]
                for next_line in following_lines:
                    if '>' not in next_line:
                        seq.append(next_line.strip())
                    else:
                        break
                sequence = "".join(seq)
                sequence = sequence.replace('*', '')
                sequences.append(sequence)

    return identifiers, sequences
# -----------------------------------------------------------------------------------------------------------
FASTA_FILE = sys.argv[1]

identifiers, sequences = get_seqs_ids_fasta(FASTA_FILE)
identifiers_clean = [ident.replace('>','').split()[0].strip() for ident in identifiers]
sequences = [seq.upper() for seq in sequences]
# -----------------------------------------------------------------------------------------------------------
count = 0
cys_count = 0
predictions = []
# -----------------------------------------------------------------------------------------------------------
CLEAVAGE_SITES = {}
# -----------------------------------------------------------------------------------------------------------
for identifier, seq in zip(identifiers, sequences):
    FILE_TMP = 'temp.fasta'
    f = open(FILE_TMP, 'w')
    f.writelines(identifier)
    f.writelines(seq.replace('*','') + '\n')
    f.close()

    # Call SignalP 3.0
    with open('signalp3.out', 'wb') as out:    
        Process = Popen(['../Software/signalp-3.0/signalp', '-t', 'euk', '-f', 'short', 'temp.fasta'], shell=False, stdout=out)
        sts = Process.wait()

    CLEAVAGE_SITES = signalp_results(CLEAVAGE_SITES)
# -----------------------------------------------------------------------------------------------------------
for ident, long_ident, seq in zip(identifiers_clean, identifiers, sequences):
    full_seq = seq
    signalp_ident = ident[:20]

    mature_protein_start = CLEAVAGE_SITES[signalp_ident]
    if mature_protein_start == 0:
        mature_protein_start = 15
    else:
        seq = seq[mature_protein_start:]

    positions = [(m.start(), m.end()) for m in re.finditer('L..R', seq)]
    positions += [(m.start(), m.end()) for m in re.finditer('[KR]R', seq)]
    positions = sorted(positions)

    site_found_for_protein = False
    cys_rich = False
    
    if positions:
        for pos_start, pos_end in positions:
            first_part = seq[:pos_start]
            remainder = seq[pos_end:]
            

            if len(first_part) >= 5 and 100.0*len(remainder)/len(seq) > 50.0:
                Nterminus = ProteinAnalysis(first_part)
                Cterminus = ProteinAnalysis(remainder)

                disorder_associated_AAS = ['K', 'E', 'N', 'S', 'P', 'G', 'R', 'D', 'Q', 'M']
                order_associated_AAS = ['W', 'Y', 'F', 'I', 'C', 'L', 'V', 'H']

                disorder_AAs_Nterminus = 0
                for aa in disorder_associated_AAS:
                    disorder_AAs_Nterminus += first_part.count(aa)
                nondisorder_AAs_Nterminus = 0
                for aa in order_associated_AAS:
                    nondisorder_AAs_Nterminus += first_part.count(aa)

                disorder_AAs_Cterminus = 0
                for aa in disorder_associated_AAS:
                    disorder_AAs_Cterminus += remainder.count(aa)
                nondisorder_AAs_Cterminus = 0
                for aa in order_associated_AAS:
                    nondisorder_AAs_Cterminus += remainder.count(aa)                        

                proportion_disorder_Nterminus = 100.0*disorder_AAs_Nterminus/(disorder_AAs_Nterminus+nondisorder_AAs_Nterminus)
                proportion_disorder_Cterminus = 100.0*disorder_AAs_Cterminus/(disorder_AAs_Cterminus+nondisorder_AAs_Cterminus)

                if proportion_disorder_Nterminus > 66.6 and proportion_disorder_Nterminus > proportion_disorder_Cterminus:
                    site_found_for_protein = True   

                    print('----------')
                    print(ident, mature_protein_start, len(full_seq))
                    print(seq[pos_start:pos_end])
                    print(mature_protein_start+pos_start+1, mature_protein_start+pos_end+1, proportion_disorder_Nterminus, proportion_disorder_Cterminus)
                                
                    if remainder.count('C') % 2 == 0 and remainder.count('C') >= 4:
                        cys_rich = True   
                                
    if site_found_for_protein:
        count += 1
        predictions.append((long_ident, full_seq))
    if cys_rich == True:
        cys_count += 1

#for pred in predictions:
#    print(pred)
#print(predictions)
print('----------')
print(FASTA_FILE),
if count > 0:
    print(len(identifiers), count, round(100.0*count/len(identifiers), 2), round(100.0*cys_count/float(count), 2))
else:
    print('None predicted')
