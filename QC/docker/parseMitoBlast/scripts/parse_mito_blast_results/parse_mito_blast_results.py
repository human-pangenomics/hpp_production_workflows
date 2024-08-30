#! /usr/bin/env python

'''
    Slightly adjusted from 2020 work by Marcela Uliano-Silva's MitoHiFi scrip parse_blast.py
    See: https://github.com/marcelauliano/MitoHiFi

    This script is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    '''
import pandas as pd
import sys

## Read in sample file info from command line
try:
    blast_result_fn = sys.argv[1]
    parsed_blast_fn = sys.argv[2]

except IndexError:
    raise SystemExit(f"{sys.argv[0]} ERROR: Check that sample name is provided")


my_names = ["qseqid", "sseqid", "pident", "alilength" , "mismatch", "gapopen", "qstart", "qend",
            "sstart", "send", "evalue" , "bitscore", "leng_query", "s_length",]


blast_cov = pd.read_csv(blast_result_fn,
                        sep="\t", names = my_names, )

#Get the percentage of the query in the blast aligment
blast_cov['%q_in_match'] = blast_cov['alilength']*100 / (blast_cov['leng_query'])

#sum percentages of query sequence in blast match based on column id
a= blast_cov.groupby('qseqid')['%q_in_match'].sum().to_frame().rename(columns={'qseqid':'%q_in_match'}).reset_index()

#get size of query and subject and drop duplicates
seqsizes = blast_cov[['qseqid', 'leng_query', 's_length']].drop_duplicates(subset='qseqid')

#merge 'a' and 'seqsizes' dataframes by 'qseqid'
result = pd.merge(a, seqsizes, on='qseqid')

# Now let's filter the blast matches
# if the lenght of the query is 10x the size of the subject (close-related mitogenome), drop it. (
# As its likely the match belongs to a NUMT
len_cutoff = (result['s_length'] * 10)
result1 = result[(result['leng_query'] < len_cutoff)].sort_values(by='%q_in_match', ascending=False)

# if the lenght of the query is 80% smaller than the length of the subject, drop it. 
# It's unlikely you will have a complete mitogenome.
result1['perc'] = result1["leng_query"]*100/(result1["s_length"])
ac=result1[result1['perc'] > 80].sort_values(by='%q_in_match')

# if the % of the query in the blast match is smaller than 70%, drop it
ac[(ac['%q_in_match'] > 70)].sort_values(by='%q_in_match', ascending=False).to_csv(parsed_blast_fn, index=False, sep="\t")


print("Parsing of blast done.")