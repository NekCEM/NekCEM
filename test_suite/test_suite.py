from subprocess import call
import numpy as np
import sys
import re
from os import path,chdir,getcwd
import json

#Check if a string is a number or not
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#Grep a file for the CPU lines, store numeric values in return array
def grep_for_CPU(file_to_read):
    opened_file = open(file_to_read,'r')

    res = []
    restmp = []
    for line in opened_file:
        if re.search("CPU",line):
            split_line = line.split()
           #Some lines (the last one) are text, not numbers. Don't include those
            if is_number(split_line[0]):
                restmp.append([(x) for x in line.split()]) 
    
    for restmp1 in restmp:
        res.append(restmp1)

    return np.array(res)

#Do a numerical diff of test and ref
def ndiff(test, ref, eps):
    tests = grep_for_CPU(test)
    refs  = grep_for_CPU(ref)

    assert tests.shape == refs.shape, "Arrays are wrong size, comparison not valid"

    max_diff = np.max(np.abs(tests[:,8:13].astype(float)-refs[:,8:13].astype(float)))

    test_results.write("Final two lines from test:\n")
    test_results.write(' '.join(tests[-2]) + '\n')
    test_results.write(' '.join(tests[-1]) + '\n')
    test_results.write("\n")
    
    if max_diff > eps:
        fail_string = ("Comparison in " + test + 
                       " larger than eps, with maximum difference of %.2g \n" %max_diff)
        test_results.write(fail_string)
        test_results.write("Final two lines from ref:\n")
        test_results.write(' '.join(refs[-2]) + '\n')
        test_results.write(' '.join(refs[-1]) + '\n')
#        assert max_diff < eps, fail_string
        return 1
    else:
        test_results.write(test + " was successful with a maximum difference of %.2g \n"%max_diff)
        return 0


#Assume we are already in the working directory
def change_rea_parameters(test,params):
    #Read in rea file
    lines = open(test + ".rea").readlines()
    
    #Loop through parameters, changing them one by one
    for param in params:
        #Because python strings are immutable, we temporarily convert to a list
        #param_num+3 is because the file is offset by 4 header lines, and python counts from 0
        tmp_line = list(lines[int(param["param_num"])+3])

        #change the parameter value
        #NOTE: THIS ASSUMES INTEGER PARAMETERS LESS THAN 9, WHICH ARE IN THE 3rd COLUMN
        #since we are changing a single
        #character in the list. A possible improvement would be supporting floats, etc
        tmp_line[2] =param["param_val"]
    
        #copy the list back to the lines
        lines[int(param["param_num"])+3] = ''.join(tmp_line)

        #write the new, changed rea file
        out = open(test + ".rea","w")
        out.writelines(lines)
        out.close()


#Clean, make, and run, followed by ndiff with reference
def make_and_run_test(test):
    original_dir = getcwd()
    #Move to test directory
    
    chdir('../' + path.dirname(test['name']))
    test_results.write("\n========================================================================\n")
    test_results.write("Running test " + test['name'] + " on arch " + test['arch'] + " with np="+str(test['np'])+"\n")
    param_string = "With .rea params: "
    for param in test["rea_params"]:
        param_string += param["param_num"] + " = " + param["param_val"] + " "

    test_results.write(param_string + "\n\n")

    #Clean all
    #call(["../../bin/cleanall"])

    #Change rea parameters
    change_rea_parameters(path.basename(test['name']),test["rea_params"])


    #Make example
    #../../bin/makenekmpi -a arch example number_of_procs
    errors = call(["../../bin/makenekmpi","-a",test['arch'],path.basename(test['name'])])
    #Check for errors

    if errors!=0:
        fail_string = "Error while making " + test['name'] + " on arch " + test['arch'] + " with np="+str(test['np'])+"\n"
        test_results.write(fail_string)

#    assert errors==0, fail_string

    #Run example
    #../../bin/nek example number_of_procs
    errors = call(["../../bin/nek",path.basename(test['name']),str(test['np'])])

 
    if errors!=0:
        fail_string = "Error while running " + test['name'] + " on arch " + test['arch'] + " with np="+str(test['np'])+"\n"
        test_results.write(fail_string)
#   assert errors==0, fail_string
    
    test_output = path.basename(test['name']) + ".np=" + str(test['np']) + ".output"
    #Diff with reference
    ref = original_dir + '/refs/' + test['ref']
    fail = ndiff(test_output,ref,float(test['eps']))

    #move back to original dir
    chdir(original_dir)
    return fail

test_results = open("results","w")

with open('list_of_tests') as data_file:    
    all_tests = json.load(data_file)

failures = 0
for test in all_tests:
    failures += make_and_run_test(test)

#Print summary information
test_results.write("========================================================================\n\n")
test_results.write(str(len(all_tests)) + " test(s) ran with " + str(failures) + " failures.\n")
