import sys, traceback
import os
import argparse
import time
import subprocess
from glob import glob
import numpy as np
import six
import pandas as pd


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
        command = 'source ~/.bashrc; ' + command # preppend sourcing bashrc to command
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

best_guess_conversion="""echo "Converting vcf files to plink using best guess (prob greater than 0.9)"
for CHR in {{1..22}}
do
  plink2 --vcf {post_imput_dir}/chr${{CHR}}.dose.vcf.gz --double-id --import-dosage-certainty 0.9 --make-bed --out {post_imput_dir}/chr${{CHR}}_imputed_plink_TEMPORARY
done

plink2 --vcf {post_imput_dir}/chrX.dose.vcf.gz --double-id --import-dosage-certainty 0.9 --make-bed --out {post_imput_dir}/chr23_imputed_plink_TEMPORARY
"""

def split_a1_a2(post_imput_dir, log):
    log.log('Splitting a1 and a2 labels from coordiantes in .bim file')
    for chrom in np.arange(1, 24):
        log.log('Processing chrom ' + str(chrom))
        chrom_file = "%s/chr%i_imputed_plink_TEMPORARY.bim" % (post_imput_dir, chrom)
        bim = pd.read_csv(chrom_file, sep='\t', header=None)
        chr_bp = bim.iloc[:, 1].str.split(':', expand=True)
        bim.iloc[:, 1] = (chr_bp.iloc[:, 0] + ':' + chr_bp.iloc[:, 1]).values
        bim.to_csv(chrom_file, sep='\t', header=None, index=False)

RemoveDup_UpdateNames = """new_name="{snp_dir}/AllChr_Sorted_Tabdelim.txt"

echo "Removing duplicate snps and renaming RSIDs"
chrom_files=($(ls {post_imput_dir}/chr*_imputed_plink_TEMPORARY.bim))
for CHR in {{1..23}}
do
    echo "Processing {post_imput_dir}/chr${{CHR}}_imputed_plink_TEMPORARY.bim"
    fileprefix={post_imput_dir}/chr${{CHR}}_imputed_plink_TEMPORARY
    plink2 --bfile ${{fileprefix}} --write-snplist --out ${{fileprefix}}_allsnps
    #Find duplciates
    snplist=${{fileprefix}}_allsnps.snplist
    dupfile=${{fileprefix}}_duplicatedsnps.snplist
    cat $snplist | sort | uniq -d > $dupfile
    # Remove duplicates
    nodup=${{fileprefix}}_NoDuplicates
    plink2 --bfile $fileprefix --exclude $dupfile --make-bed --out ${{nodup}}
    # Update Name to RSIDs
    rsidout={post_imput_dir}/chr${{CHR}}_imputed_plink_RSID
    plink2 --bfile ${{nodup}} --update-name ${{new_name}} --make-bed --out ${{rsidout}} &
done"""

def fix_fam_IDs(post_imput_dir, log):
    # TOPMED appends IID and FID into one ID, this function checks if IDS are of format IID_FID and fixes them to be IID
    rsidfam = f'{post_imput_dir}/chr{{chrom}}_imputed_plink_RSID.fam'
    for chrom in np.arange(1, 24):
        fam_temp = rsidfam.format(chrom=chrom)
        fam = pd.read_csv(fam_temp, sep='\t', header=None)
        log.log('Checking fam IDs from ' + fam_temp)
        if fam.iloc[0, 0].count('_')==3:
            log.log('IDs of wrong format (IID_FID) from TOPMED: ' +  fam.iloc[0, 0])
            fix_IDs = lambda x: x.split('_')[0] + '_' + x.split('_')[1]
            fam.iloc[:, 0] = fam.iloc[:, 0].apply(fix_IDs)
            fam.iloc[:, 1] = fam.iloc[:, 1].apply(fix_IDs)
            log.log('Writing corrected IDs to ' + fam_temp)
            fam.to_csv(fam_temp, sep='\t', header=False, index=False)

def merge_files(post_imput_dir, log):
    # Function to merge files
    merge_file = f'{post_imput_dir}/merge.list'
    rsidfile = f'{post_imput_dir}/chr{{chrom}}_imputed_plink_RSID'
    files = pd.Series(np.arange(1, 24)).apply(lambda x: rsidfile.format(chrom=x))
    files.to_csv(merge_file, index=False, header=False)
    merge_out = f'{post_imput_dir}/merged_chroms'
    log.log('Concatonating files to ' + merge_out)
    log.system(f'plink  --merge-list {merge_file} --make-bed --out {merge_out}')

def main(post_imput_dir, snp_dir, keep_temp, log):
    
    # Convert to plink using best guess (prob>0.9)
    log.system(best_guess_conversion.format(post_imput_dir=post_imput_dir))
    
    # Split a1 and a2 labels from the bim co-ordinate columns
    split_a1_a2(post_imput_dir, log)
    
    # Remove duplicates and update coordiantes to RSID numbers
    log.system(RemoveDup_UpdateNames.format(post_imput_dir=post_imput_dir, snp_dir=snp_dir))
    
    # Remove TEMPORARY files
    if not keep_temp:
        log.log('Removing TEMPORARY files')
        temp_files = pd.Series(glob(f'{post_imput_dir}/*TEMPORARY*'))  
        temp_files = temp_files.loc[np.logical_not(temp_files.str.endswith('.log'))] # Retain TEMPORARY log files
        for temp_file in temp_files:
            log.system(f'rm {temp_file}')
    
    # Fix fam IDs
    fix_fam_IDs(post_imput_dir=post_imput_dir, log=log)

    # Merged files
    merge_files(post_imput_dir=post_imput_dir, log=log)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This function converts vcf files dowloaded from TOPMED imputation server to plink files and updates names to RSIDs (if availible)')
    parser.add_argument('-post_imput_dir', help='Filepath to directory containing imputation results downloaded from TOPMED server', type=str, default=None)
    parser.add_argument('-snp_dir', help='Directory containing output from prep_RSID_files (hg19 for HRC and hg38 TOPMED reference)', type=str, default='TOPMED')
    parser.add_argument('-keep_temp', help='Keep temporary files uses up disk space (used for debugging)', action='store_true')

    opt = parser.parse_args()
    start_time = time.time()

    try:
        # Logging    
        defaults = vars(parser.parse_args([]))
        opts = vars(opt)
        non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
        log = Logger(opt.post_imput_dir + '/plink_conversion.log', 'w')
        header = "Call: \n"
        header += './plink_conversion.py \\\n'
        options = ['\t--'+x.replace('_','-')+' '+str(opts[x]).replace('\t', '\\t')+' \\' for x in non_defaults]
        header += '\n'.join(options).replace('True','').replace('False','')
        header = header[0:-1]+'\n'
        log.log(header)


        main(opt.post_imput_dir, opt.snp_dir, opt.keep_temp, log)

    except Exception:
        ex_type, ex, tb = sys.exc_info()
        log.error( traceback.format_exc(ex) )
        raise

    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()) )
        time_elapsed = round(time.time()-start_time,2)
        log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))