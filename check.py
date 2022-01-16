import os
import signal
import csv
import re
import argparse
import tempfile
from subprocess import Popen, PIPE
import sys

c_blue   = '\033[49;0;34m'
c_red    = '\033[49;0;31m'
c_green  = '\033[49;0;32m'
c_cyan   = '\033[49;0;36m'
c_yellow = '\033[49;1;33m'
c_reset  = '\033[0m'

def green(string):
   return c_green + string + c_reset

def red(string):
   return c_red + string + c_reset

def yellow(string):
   return c_yellow + string + c_reset

def log_progress(msg):
   sys.stdout.write(msg)


def run(args):
   os.system('rm -rf ' + tempfile.gettempdir() + '/*')
   path = os.walk(args.benchmark_dir)
   i = 0;
   for root, _, files in path:
      for file in files:
         if (args.file_name != None):
            if (file != args.file_name): continue

         with open('check/mydo', 'w') as f:
            if (i == 0):
               f.write('cirr ' + os.path.join(root, file) + '\n')
            else:
               f.write('cirr ' + os.path.join(root, file) + ' -r\n')
            f.write('cirrar\n')
            if args.check:
               f.write('cirrarw -d ' + tempfile.gettempdir()+'\n')
            i += 1
            f.write('q -f\n')
         
         cmd = './cirTest -F check/mydo'
         process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE, preexec_fn=os.setsid)
         try:
            output, err = process.communicate(timeout=600)
            if (err != b''):
               log_progress(red('ERROR: '))
            else:
               log_progress(green('SUCCESS: '))
            with open(os.path.join(args.result_dir, file[:-4]+'.out'), 'w') as f:
               f.write(output.decode('UTF-8'))
         except:
            os.killpg(process.pid, signal.SIGINT)
            log_progress(yellow('TIMEOUT: '))
         log_progress(os.path.join(root, file) + '\n')

         if (args.check):
            os.chdir('check')
            tmp_path = os.walk(tempfile.gettempdir())
            for tmp_root, _, tmp_files in tmp_path:
               for tmp_file in tmp_files:
                  log_progress('    Correctness check on '+ tmp_file + '... ')
                  tmp_cmd = './aigmiter -o ' + 'miter.aag ' + os.path.join('..', root, file) + ' ' + os.path.join(tmp_root, tmp_file) + '\n'
                  tmp_cmd += './aigtocnf miter.aag miter.dimacs\n'
                  tmp_cmd += 'minisat miter.dimacs'
                  process_tmp = Popen(tmp_cmd, shell=True, stdout=PIPE, stderr=PIPE, preexec_fn=os.setsid)
                  output, err = process_tmp.communicate(timeout=20)
                  text = output.decode('UTF-8')
                  if (text.find('UNSATISFIABLE') != -1):
                     log_progress(green('success\n'))
                  elif (text.find('SATISFIABLE') != -1):
                     log_progress(red('fail'))
                     log_progress(', miter returns SAT\n')
                  else:
                     log_progress(red('fail'))
                     log_progress(', other errors\n')
            os.chdir('..')

         os.system('rm -rf ' + tempfile.gettempdir() + '/*')


if __name__ == "__main__":
   parser = argparse.ArgumentParser()
   parser.add_argument('--check', action='store_true', default=False)
   parser.add_argument('--benchmark_dir', type=str, required=True)
   parser.add_argument('--result_dir', type=str, default=None)
   parser.add_argument('--file_name', type=str, default=None)
   args = parser.parse_args()

   run(args)