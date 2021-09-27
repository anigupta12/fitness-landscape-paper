import os
from Bio import SeqIO

fwd_seq = 'GATAAACGGTACGCTGAGGG'
rev_seq = 'CTCAGACCACGCTGATGCCC'

bad_read = 'bad read'

def extract_combined_region(r1,r2):
    # find primer regions, remove primers and everything before, take rev_complement of r2, concatenate reads
    library_r1_cutoff = 72 # could adjust these, just have to be after last mutant position
    library_r2_cutoff = 72

    r1_pos = r1.find(fwd_seq)
    if r1_pos >=0:
        found_fwd = True
    else: found_fwd = False

    r2_pos = r2.find(rev_seq)
    if r2_pos >=0:
        found_rev = True
    else: found_rev = False

    if found_fwd and found_rev:
        lib_start = r1_pos+len(fwd_seq)
        lib_end = lib_start + library_r1_cutoff
        library_r1 = r1[lib_start:lib_end]

        lib_start = r2_pos+len(rev_seq)
        lib_end = lib_start + library_r2_cutoff
        library_r2 = r2[lib_start:lib_end]
        # print 'L1', library_r1
        # print 'L2', library_r2
        # print 'L2 RC', library_r2.reverse_complement()
        combined=library_r1+'-'+library_r2.reverse_complement

        # todo: check quality score are above threshold!

    else:
        combined=bad_read
        # print 'r1 bad', found_fwd, r1
        # print 'r2 bad', found_rev, r2
    return combined

def tabulate(key, combined_counts):
    # add 1 to count for this key, or add it to the dictionary
    if key in combined_counts:
        combined_counts[key] += 1
    else:
        combined_counts[key] = 1

def parse_both_files(f1, f2):
    # read 2 fastq files, parse the reads and return counts
    combined_counts = {}; # holds the counts for each library region
    read1_gen = SeqIO.parse(f1, 'fastq')
    read2_gen = SeqIO.parse(f2, 'fastq')
    
    reads_processed_counts = 0
    good_read_counts = 0
    bad_read_counts = 0

    while True:
        try:
            read1= read1_gen.next()
            read2= read2_gen.next()
            reads_processed_counts += 1
            # print 'read 1', read1.seq
            # print 'read 2', read2.seq
            combined = extract_combined_region(read1.seq,read2.seq) # todo: pass in quality scores
            if combined is not bad_read:
                good_read_counts += 1
                tabulate(str(combined), combined_counts)
            else:
                bad_read_counts += 1

            # if reads_processed_counts % 10000 is 0:
            #     print reads_processed_counts
            #     # break

        except StopIteration:
            # print 'end of files'
            break
    print 'processed {p} reads, good = {g}, bad = {b}'.format(
        p = reads_processed_counts, 
        g = good_read_counts, 
        b = bad_read_counts)
    print 'found {n} unique sequences\n'.format(n=len(combined_counts))
    return combined_counts

def write_counts(dict, outputfile):
    # write sequences and counts to a text file
    file = open(outputfile, 'w')
    for w in sorted(dict, key=dict.get, reverse=True):
        file.write('{seq}, {num}\n'.format(seq=w, num=dict[w]))
    file.close()

def parse_folder(folderpath):
    # parse all files in a given folder, output to 'folderpath_parsed' folder
    library_names = []
    # get all library names from R1 files
    for fname in os.listdir(folderpath):
        pos = fname.find('_R1_001.fastq')
        if pos != -1:
            lib_name = fname[:pos]
            # print lib_name
            if lib_name not in library_names:
                library_names.append(lib_name)
            # else:
            #     print 'ERROR: duplicate library name:', lib_name
    # parse each library
    for lib_name in library_names:
        r1_file = '{folder}/{lib}_R1_001.fastq'.format(folder=folderpath, lib=lib_name)
        r2_file = '{folder}/{lib}_R2_001.fastq'.format(folder=folderpath, lib=lib_name)
        outputfile = '{folder}_parsed/{lib}.txt'.format(folder=folderpath, lib=lib_name[:-9]) # remove last 8 characters of name

        print "parsing {lib} files".format(lib=lib_name)
        counts = parse_both_files(r1_file, r2_file)
        write_counts(counts, outputfile)

 #### MAIN BITS #####

foldername = '27jun18_parsed' # Enter Folder name
parse_folder(foldername)