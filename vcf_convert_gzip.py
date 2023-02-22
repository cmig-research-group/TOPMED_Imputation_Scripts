import os
import sys, traceback
import argparse
import time
import six

## Function to seperate plink file into seperate chromosomes and gzip for upload to Michigan Imputation Server
## Usage
# python vcf_convert_gzip.py -bed_file /path/to/merged_bed_file_prefix -out_file /path/to_output/vcf/ABCD_release_2020_10_chr

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

        os.system(command)
 
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


def main(bed_file, out_file, log):
    
    log.log('Converting plink files to vcf format')
    for chrom in list(range(1, 23)) + ['x']:
        # Seperate into chromomosome files
        log.log(f'Creating seperate file for chrom {chrom} at {out_file + str(chrom)}')
        log.system('plink --bfile ' + bed_file + ' --chr ' + str(chrom) + ' --recode --out ' + out_file + str(chrom)) 
        # log.system('plink --bfile ' + bed_file + ' --chr ' + str(chrom) + ' --output-chrom M --recode --out ' + new_text_fileset + str(chrom)) # Not sure why this had output-chrom M
        ped_file = out_file + str(chrom) + '.ped'
        map_file = out_file + str(chrom) + '.map'
        vcf = out_file + str(chrom)
        # Convert to VCF
        log.log(f'Creating vcf file at {vcf}')
        log.system('plink --ped ' + ped_file + ' --map ' + map_file + ' --recode vcf --out ' + vcf)
        script_dir = os.path.dirname(os.path.abspath(__file__))
        log.system('perl ' + script_dir + '/vcf-sort ' + vcf + '.vcf | bgzip -c > ' + vcf + '.vcf.gz')
        log.log('Removing temporary files')
        nosex_file = map_file.replace('map', 'nosex')
        log.system(f'rm {ped_file} {map_file} {vcf + ".vcf"} {nosex_file}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This function converts a merged bed file to seperate chromosome gziped vcf files for upload to Michigan Imptuation Server')
    parser.add_argument('-bed_file', help='Filepath to merged bed file', type=str, default=None)
    parser.add_argument('-out_file', help='Output file, chromosome number will be appended to this for each file', type=str, default=None)

    # raw_genotype = '/space/gwas-syn2/1/data/GWAS/ABCD/genotype/release3.0/genotype_QCed'
    # bed_file = os.path.join(raw_genotype, 'ABCD_release_3.0_QCed')
    # out_dir = '/space/gwas-syn2/1/data/GWAS/ABCD/genotype_proc/imputation/HRC_2020_10/pre_proc'
    # out_file = out_dir + '/ABCD_release_2020_10_chr'

    opt = parser.parse_args()

    try:
        # Logging
        start_time = time.time()
        defaults = vars(parser.parse_args([]))
        opts = vars(opt)
        non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
        log = Logger(opt.out_file + '.log', 'w')
        header = "Call: \n"
        header += './vcf_convert_gzip.py \\\n'
        options = ['\t--'+x.replace('_','-')+' '+str(opts[x]).replace('\t', '\\t')+' \\' for x in non_defaults]
        header += '\n'.join(options).replace('True','').replace('False','')
        header = header[0:-1]+'\n'
        log.log(header)


        main(opt.bed_file, opt.out_file, log)

    except Exception:
        ex_type, ex, tb = sys.exc_info()
        log.error( traceback.format_exc(ex) )
        raise

    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()) )
        time_elapsed = round(time.time()-start_time,2)
        log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))