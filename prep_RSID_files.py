import os
import sys, traceback
import argparse
import time
import six
import pandas as pd
import numpy as np
import subprocess


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def sec_to_str(t):
    '''Convert seconds to days:hours:minutes:seconds'''
    [d, h, m, s, n] = six.moves.reduce(lambda ll, b : divmod(ll[0], b) + ll[1:], [(t, 1), 60, 60, 24])
    f = ''
    if d > 0:
        f += '{D}d:'.format(D=d)
    if h > 0:
        f += '{H}h:'.format(H=h)
    if m > 0:
        f += '{M}m:'.format(M=m)

    f += '{S}s'.format(S=s)
    return f

class Logger(object):
    '''
    Lightweight logging.
    '''
    def __init__(self, fh, mode):
        self.fh = fh
        self.log_fh = open(fh, mode) if (fh is not None) else None

        # remove error file from previous run if it exists
        try:
            os.remove(fh + '.error')
        except OSError:
            pass

    def system(self, command, timeout=None):  # timeout in seconds, None means no timeout.
        start_time = time.time()
        log.log('>command start time: {T}'.format(T=time.ctime()) )
        self.log(command)        

        # Required to use bash shell
        subprocess.call('/bin/bash -c "$PROC"', shell=True, env={'PROC': command})
 
        time_elapsed = round(time.time()-start_time,2)
        log.log('=command elapsed time: {T}'.format(T=sec_to_str(time_elapsed)))
        log.log('<command end time: {T}'.format(T=time.ctime()) )        

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.
        '''
        eprint(msg)
        if self.log_fh:
            self.log_fh.write(str(msg).rstrip() + '\n')
            self.log_fh.flush()

    def error(self, msg):
        '''
        Print to log file, error file and stdout with a single command.
        '''
        eprint(msg)
        if self.log_fh:
            self.log_fh.write(str(msg).rstrip() + '\n')
            with open(self.fh + '.error', 'w') as error_fh:
                error_fh.write(str(msg).rstrip() + '\n')

reduce_151 = """snp_file="snp151.txt.gz"
snp_out="snp151_rsID_TEMPORARY.txt"

zcat ${{snp_file}} | awk -F'\t' '{{ print $2":"$4"\t"$5 }}' > ${{snp_out}}
rm ${{snp_file}}"""

def split_to_chroms(log):
    log.log('Seperating snp151_rsID_TEMPORARY.txt file into seperate RSID files')
    snps = pd.read_csv('snp151_rsID_TEMPORARY.txt', sep='\t', header=None)
    for chrom in np.arange(1, 23).tolist() + ['X', 'Y']:
        print('Processing ' + str(chrom))
        chrm_start = 'chr' + str(chrom) + ':'
        chrom_ind = snps.iloc[:, 0].str.startswith(chrm_start).values
        out_file = './chr' + str(chrom) + '_snp151_rsID_TEMPORARY.txt'
        snps.iloc[chrom_ind, :].to_csv(out_file, sep='\t', index=False, header=False)

remove_dups = """echo "Removing duplicate SNPs from rsID files"

chrom_files=($(ls chr*_snp151_rsID_TEMPORARY.txt))

# Remove Duplicates
for chrfile in ${{chrom_files[@]}}
do
    echo "Processing ${{chrfile}}"
    awk '{k=($1"\t"$2)} {{a[$1]++;b[$1]=k}}END{for (x in a) if (a[x]==1)print b[x]}' ${{chrfile}} > "${{chrfile/snp151_rsID/nodups}}"
done

# concatonate to single file
chrom_dups=($(ls chr*_nodups_TEMPORARY.txt))

rm -f AllChr_Sorted_Tabdelim.txt # Remove all chroms file if they exist
for chrom_dup in ${{chrom_dups[@]}}; do cat $chrom_dup >> AllChr_Sorted_Tabdelim.txt; done"""


def main(snp_dir, build, keep_temp, log):
    
    os.chdir(snp_dir)
    
    # Download
    log.log(f'Download SNP file for {build}')
    log.system(f'wget http://hgdownload.soe.ucsc.edu/goldenPath/{build}/database/snp151.txt.gz')
    
    # Reduce downloaded file
    log.log('Reducing columns of snp151')
    log.system(reduce_151.format())
    
    # Split chroms
    split_to_chroms(log)
    
    # Remove Duplicates
    log.system(remove_dups.format())
    
    # Remove temporary files
    if not keep_temp:
        log.log('Removing TEMPORARY files')
        log.system(f'rm {snp_dir}/*TEMPORARY*')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This function downloads and prepares RSID files for updating SNP names after imputation')
    parser.add_argument('-snp_dir', help='Directory to download SNP RSID file', type=str, default=None)
    parser.add_argument('-build', help='Genome build for RSIDs to be downloaded (hg38 or hg19, for TOPMED and HRC respectively)', type=str, default='hg38', choices=['hg19', 'hg38'])
    parser.add_argument('-keep_temp', help='Keep temporary files uses up disk space (used for debugging)', action='store_true')

    opt = parser.parse_args()

    try:
        # Logging
        start_time = time.time()
        defaults = vars(parser.parse_args([]))
        opts = vars(opt)
        non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
        log = Logger(opt.snp_dir + '/prep_RSID_files.log', 'w')
        header = "Call: \n"
        header += './prep_RSID_files.py \\\n'
        options = ['\t--'+x.replace('_','-')+' '+str(opts[x]).replace('\t', '\\t')+' \\' for x in non_defaults]
        header += '\n'.join(options).replace('True','').replace('False','')
        header = header[0:-1]+'\n'
        log.log(header)


        main(opt.snp_dir, opt.build, opt.keep_temp, log)

    except Exception:
        ex_type, ex, tb = sys.exc_info()
        log.error( traceback.format_exc(ex) )
        raise

    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()) )
        time_elapsed = round(time.time()-start_time,2)
        log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))