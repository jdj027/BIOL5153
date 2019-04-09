#! /usr/bin/env python3
# def py_sub_ahpcc(job_name,queue,ppn,walltime):
# convert this to a function at some point
# Ex:   $ ./py_sub_AHPCC kiwi -q=aja -p=2 -w=3


import argparse
# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description = "generates a PBS job script for the AHPCC Trestles cluster")

# add positional arguments
parser.add_argumen("job_name", help="the name of your job")

# add optional arguments
parser.add_argument("-q", "--queue", help='name of the queue in Trestles (default is q06h32c)', default='q06h32c')
parser.add_argument("-p", "--ppn", help='number of processors (default is 32)', type=int, default=32)
parser.add_argument("-w", "--walltime", help='amount of time needed for job (in hours) (default is 6)', type=int, default='6')

# parse the command line
args = parser.parse_args()

print('#PBS -N ' +  args.job_name)
print('#PBS -q ' + args.queue)
print('#PBS -j oe')
print('#PBS -o ' +  args.job_name + '.$PBS_JOBID')
print('#PBS -l ' + ' nodes=1:ppn=' + str(args.ppn))
print('#PBS -l ' + ' walltime=' + str(args.walltime) + ':00:00')
print()
print('cd $PBS_O_WORKDIR')
print('module purge')
print('module load python/3.6.0-anaconda')

# job_name = 'kiwi' # name the job
# queue = 'aja'	    # name of queue
# ppn = '1'	    # number of processors
# walltime = '1'    # in hours

# print('#PBS -N ' +  job_name)
# print('#PBS -q ' +  queue)
# print('#PBS -j oe')
# print('#PBS -o ' +  job_name + '.$PBS_JOBID')
# print('#PBS -l ' + ' nodes=1:ppn=' + ppn)
# print('#PBS -l ' + ' walltime=' + walltime + ':00:00')
# print()
# print('cd $PBS_O_WORKDIR')
# print('module purge')
# print('module load python/3.6.0-anaconda')
