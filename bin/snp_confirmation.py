#!/usr/bin/env python3

import csv
import os
import sys
import re

##  USAGE:python snp_confirmation.py S1.sam amr_analytic_matrix.csv snp_metadata_example.csv indel_count.csv long long_HMM.csv

def header_collect(long_file, snp_metadata):
    """This method reads in the headers from a given AMR_matrix and a SNP_metadata file and
    pulls out the common headers between the two along with the info associated with the gene"""
    matrix_headers = {}
    snp_headers = {}
    temp_list = []
    with open(long_file, 'r') as csvfile:  # pulls the headers from the AMR_matrix, delimiting by tab
        reader = csv.reader(csvfile, delimiter='\t')
        csvfile.readline()
        for row in reader:
            #print(row)
            if row != []:
                matrix_headers[row[1]] = ''
    csvfile.close()
    with open(snp_metadata, 'r') as csvfile: # pulls the headers from the SNP_metadata, delimiting by comma
        reader = csv.reader(csvfile, delimiter=',')
        csvfile.readline()
        for row in reader:
            #print(row)
            if row != []:
                if matrix_headers.__contains__(row[0]):  # pulls out the headers if they are contained in the AMR_matrix
                    snp_headers[row[0]] = []
        csvfile.seek(0)
        for gene in snp_headers:
            for row in reader:
                #print(row)
                if row != []:
                    if gene == row[0]:
                        temp_list.append([row[2], row[3], row[4], row[5], row[10]])
            snp_headers[gene] = temp_list
            temp_list = []
            csvfile.seek(0)
    csvfile.close()
    return snp_headers


def indel_count_create(output_file, indel_counts):
    """This method creates the necessary indel count table"""
    try:
        #add given sample's indel counts since the file exists
        read_indelfile = open(output_file)
        read_indelfile.close()
        with open(output_file, 'a') as append_indelfile:
            appender = csv.writer(append_indelfile, delimiter=',')
            print(indel_counts)
            for row in indel_counts:
                print(row)
                appender.writerow(row)
        append_indelfile.close()
    except IOError:
        #create indel count table from scratch
        indel_header = ['Sample', 'Gene', 'indel_no_snp_yes', 'indel_no_snp_no', 'indel_yes']
        with open(output_file, 'wb',) as write_indelfile:
            edit = csv.writer(write_indelfile, delimiter=',')
            edit.writerow(indel_header)
            for row in indel_counts:
                edit.writerow(row)
        write_indelfile.close()


def long_table_zero(long_file, new_long_file, indel_counts):
    """This method sets the appropriate counts to zero"""
    copy = []
    with open(long_file, 'r') as long_file:
        reader = csv.reader(long_file, delimiter=',')
        for row in reader:
            copy.append(row)
    long_file.close()

    for row in indel_counts:
        if row[2] == 0:
            for element in copy:
                if element[1] == row[1]:
                    element[2] = 0

    with open(new_long_file, 'wb') as new_long_file:
        writer = csv.writer(new_long_file, delimiter=',')
        for row in copy:
            writer.writerow(row)
    new_long_file.close()


def get_sample_name(long_file):
    """Get the sample name for later use"""
    with open(long_file, 'r') as long_file:
        reader = csv.reader(long_file, delimiter='\t')
        long_file.readline()
        for row in reader:
            if row[0] != []:
                sample_name = row[0]
                break
    long_file.close()
    return sample_name


def start_stop_to_one_based(dictionary):
    """Convert start/stop indices given by the SNP metadata to int instead of str,
    and add one to make them 1 based instead of 0 based"""
    for lists in dictionary:
        for each_list in range(len(dictionary[lists])):
            for start_stop in range(0, 2):
                dictionary[lists][each_list][start_stop] = (int(dictionary[lists][each_list][start_stop]))+1
    return dictionary


def cigar_str_list(cig_list):
    """Put str in the cigar_list into their own list"""
    str_list = []
    for x in cig_list[1::2]:
        str_list.append(x)
    return str_list


def cigar_int_list(cig_list):
    """Convert int strings in the cigar_list to int and put them into a list"""
    int_list = []
    for x in cig_list[::2]:
        int_list.append(int(x))
    return int_list


def cigar_count(cig_str, cig_int):
    """Count of how long the sequence should be based on the given cigar string"""
    count = 0
    y = 0
    for x in cig_str:
        if x == 'M':
            count += (cig_int[y])
        elif x == 'I':
            count += (cig_int[y])
        elif x == 'S':
            count += (cig_int[y])
        elif x == '=':
            count += (cig_int[y])
        elif x == 'X':
            count += (cig_int[y])
        y += 1
    return count


def cig_edit_read(read, cigar_str, cigar_int):
    for operator in cigar_str:
        if operator == 'D':
            loops = cigar_int[cigar_str.index('D')]
            index = 0
            for x in range(0, cigar_str.index('D')):
                index += cigar_int[x]
            while loops > 0:
                read.insert(index, '')
                loops -= 1
        elif operator == 'N':
            loops = cigar_int[cigar_str.index('N')]
            index = 0
            for x in range(0, cigar_str.index('N')):
                index += cigar_int[x]
            while loops > 0:
                read.insert(index, '')
                loops -= 1
    return read


def cig_edit_data_base(data_base, cigar_str, cigar_int):
    for operator in cigar_str:
        if operator == 'I':
            loops = cigar_int[cigar_str.index('I')]
            index = 0
            for x in range(0, cigar_str.index('I')):
                index += cigar_int[x]
            while loops > 0:
                data_base.insert(index, '')
                loops -= 1
    return data_base


def protein_identifier(ref_codon):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    codon = ''
    for i in ref_codon:
        codon += i
    return table[codon]


class SamParser:
    """This object takes as input a SAM file path and constructs an iterable that outputs
     sequence information.  Only one line will be held in memory at a time using this method.
    """
    def __init__(self, filepath):
        """
        constructor
        @param filepath: filepath to the input raw SAM file.
        """
        if os.path.exists(filepath):  # if file is a file, read from the file
            self.sam_file = str(filepath)
            self.stdin = False
        elif not sys.stdin.isatty():  # else read from standard in
            self.stdin = True
        else:
            raise ValueError("Parameter filepath must be a SAM file")
        self.current_line = None
        self.reads_mapping = 0
        self.reads_total = 0
        # Allowable bitflags for SAM file -> reads with both mates mapping, regardless of other flags
        self.true_flags = (99, 147, 83, 163, 67, 131, 115, 179, 81, 161, 97, 145, 65, 129, 113, 177)

    def __iter__(self):
        return self

    def _iterate(self):
        # Skip all leading whitespace
        while True:
            if self.stdin:
                sam_line = sys.stdin.readline()  # read from stdin
            else:
                sam_line = self.sam_file.readline()  # read from file
            if not sam_line:
                return  # End of file
            if sam_line.startswith("@SQ"):  # these lines contain refseq length information
                temp = sam_line.split()
                return temp[1][3:], temp[2][3:]
            elif sam_line[0] != "@":  # these lines are the actual reads
                self.reads_total += 1
                if self.reads_total % 100000 == 0:  # update the counter on stdout every 100000 reads
                    sys.stdout.write("\rReads processed: {}".format(self.reads_total))
                    sys.stdout.flush()
                temp = sam_line.split()
                if int(temp[1]) in self.true_flags and temp[2] is not "*" and int(temp[3]) is not 0:
                    self.reads_mapping += 1
                    return temp[1], temp[2], temp[3], temp[5], temp[9]

    def next(self):
        if not self.stdin and type(self.sam_file) is str:  # only open file here if sam_file is a str and not file
            self.sam_file = open(self.sam_file, "r")
        value = self._iterate()
        if not value:  # close file on EOF
            if not self.stdin:
                self.sam_file.close()
            # Steven had this print but I do not believe it necessary
            # print("{0} with both mates mapped out of {1} total reads\n".format(self.reads_mapping, self.reads_total))
            raise StopIteration()
        else:
            return value


if __name__ == '__main__':
    S = SamParser(sys.argv[1])
    snp_data = header_collect(sys.argv[2], sys.argv[3])  # dictionary containing the snp_data
    snp_data = start_stop_to_one_based(snp_data)
    final_indel = []
    set()
    values = S.next()
    while values:
        #print(snp_data)
        print(values[1])
        if snp_data.__contains__(values[1]):  # checks if the dictionary contains the given header
            cigar_list = list(filter(None, re.split('(\d+)', values[3])))  # get cigar string and divide it into numbers/letters
            cigar_int = cigar_int_list(cigar_list)
            cigar_str = cigar_str_list(cigar_list)
            read_length = cigar_count(cigar_str, cigar_int)
            for x in snp_data:
                #print(values[1])
                if x == values[1]:
                    indel_counts = [0, 0, 0]
                    for y in range(0, len(snp_data[x])):
                        start_index_both = int(values[2])
                        stop_index_both = int(values[2]) + read_length
                        start_snp = snp_data[x][y][0]
                        stop_snp = snp_data[x][y][1]

                        if (start_index_both <= start_snp) and (stop_index_both >= stop_snp):
                            data_base = list(snp_data[x][y][4][start_index_both:stop_index_both])
                            read = list(values[4])

                            data_base = cig_edit_data_base(data_base, cigar_str, cigar_int)
                            read = cig_edit_read(read, cigar_str, cigar_int)

                            indel_present = False

                            for base in range((start_snp-start_index_both)+(stop_snp-start_snp)):
                                if (data_base[base] == '') or (read[base] == ''):
                                    indel_present = True
                                    #indel yes
                                    indel_counts[2] += 1
                                    break

                            if indel_present == True:
                                continue

                            data_base_snp = data_base[(start_snp-start_index_both):(start_snp-start_index_both)+(stop_snp-start_snp)]
                            read_snp = read[(start_snp-start_index_both):(start_snp-start_index_both)+(stop_snp-start_snp)]

                            if data_base_snp != read_snp:
                                if str(protein_identifier(read_snp)) == str(snp_data[x][y][2]):
                                    #indel no snp no
                                    indel_counts[1] += 1
                                else:
                                    #indel no snp yes
                                    indel_counts[0] += 1
                            else:
                                #indel no snp no
                                indel_counts[1] += 1

                    one_same = False
                    if indel_counts != [0, 0, 0]:
                        indel_counts.insert(0, x)
                        if final_indel == []:
                            final_indel.append(indel_counts)
                            continue
                        else:
                            for entry in final_indel:
                                if entry[0] == indel_counts[0]:
                                    one_same = True
                                    for element in range(1, len(indel_counts)):
                                        entry[element] += indel_counts[element]
                                    break

                        if one_same == False:
                            final_indel.append(indel_counts)

        try:
            values = S.next()
        except StopIteration:
            sample_name = get_sample_name(sys.argv[2])
            for entry in final_indel:
                entry.insert(0, sample_name)
            indel_count_create(sys.argv[4], final_indel)
            long_table_zero(sys.argv[2], sys.argv[5], final_indel)
            break

