from Bio import SeqIO
import sys

FASTA_FILE_PATH = "/mnt/scratch/linpaul1/1uniprot_trembl.fasta"

def analyze_fasta_file(file_path):
    sequence_count = 0
    total_length = 0
    gc_count = 0
    at_count = 0
    
    first_two_records = []

    try:
        for record in SeqIO.parse(file_path, "fasta"):
            sequence_count += 1
            seq_len = len(record.seq)
            total_length += seq_len
            
            gc_count += record.seq.upper().count('G')
            gc_count += record.seq.upper().count('C')

            at_count += record.seq.upper().count('A')
            at_count += record.seq.upper().count('T')
            
            if len(first_two_records) < 30:
                print(record.seq)
            

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        return None

    avg_length = total_length / sequence_count if sequence_count > 0 else 0
    gc_content = (gc_count / total_length) * 100 if total_length > 0 else 0
    at_content = (at_count / total_length) * 100 if total_length > 0 else 0

    return {
        "sequence_count": sequence_count,
        "total_length": total_length,
        "average_length": avg_length,
        "gc_content": gc_content,
        "at_content": at_content,
        "first_two_records": first_two_records
    }

analysis_results = analyze_fasta_file(FASTA_FILE_PATH)

print(f"Total Sequences: {analysis_results['sequence_count']:,}")
print(f"Total Length: {analysis_results['total_length']:,} bases")
print(f"Average Sequence Length: {analysis_results['average_length']:.2f} bases")
print(f"Overall GC Content: {analysis_results['gc_content']:.2f}%")
print(f"Overall AT Content: {analysis_results['at_content']:.2f}%")

for i, record in enumerate(analysis_results["first_two_records"]):
    print(f"Sequence {i+1} ID: {record.id}")
    print(f"Sequence {i+1} Length: {len(record.seq)}")
    print(f"Sequence {i+1} Snippet: {str(record.seq)[:50]}...")

