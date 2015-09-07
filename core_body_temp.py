# Mouse CBT analysis
#
# Lee Organick 
# Steiner Lab
#
#

##import cProfile, pstats, StringIO
##pr = cProfile.Profile()
##pr.enable()
##import matplotlib.ticker as plticker
##import matplotlib.dates as md
##import dateutil
##import copy
##import pylab as pl
##import matplotlib.patches as mpatches

import os
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import math
pi = math.pi
import scipy
from scipy import stats
import collections

def clean_all_csv_files(all_csv_data_files):
    """Given a csv file, creates a new .csv file with "clean_" preceding the old .csv filename.
    Newly created file is identical to the original, but all rows before the occurence of the phrase
    "Deg. C Date" are deleted. This function is necessary for later uses of dictreader."""
    for filename in all_csv_data_files:
        if "clean" in filename:
            pass
        else:
            count = 0
            f = csv.reader(open(filename, "rU")) #opens the original file given
            write = csv.writer(open("clean_" + filename, "wb"),delimiter=',') #creates new, empty file
            for row in f:
                #Once the phrase "Veh Deg. C Date" occurs, takes that row and all rows below it
                if any("Deg. C Date" in s for s in row):
                    count += 1
                if count > 0:
                    write.writerow(row)

def get_data_file_names(): 
    """Returns list of strings of .csv/.TXT data files in the directory the code is being run in that
    don't have "user" or "test" in their title but do have '.csv', 'CBT ' (NOTE space after 'CBT'!!)
    or 'Proper' in their title."""
    data_files = []
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    for f in files:
        if 'user' in f:
            pass
        elif 'test' in f:
            pass
        elif '.csv' in f:
            data_files.append(f)
        elif 'Proper' in f:
            data_files.append(f)
        elif 'CBT ' in f:###This seems to change, make variable that user changes in user_modify?
            data_files.append(f)
    return data_files

def get_clean_data_file_names():
    """Returns a list of strings of .csv files in the directory that have the phrase "clean_" in
    the title. This is needed to collect data only from the clean .csv data files so we can then
    utilize dictreader."""
    clean_data_files = []
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    for f in files:
        if 'clean_' in f:
            clean_data_files.append(f)
    return clean_data_files

def get_all_mouse_ids_csv(filenames):
    """Given a list of csv files in which the first row of each file contains the mouse ids in the
    format "# Deg. C Data", returns a set of mouse ids as they appear in the file."""
    mouse_ids = set() #type is set() to ensure each id only appears once
    for f in filenames:
        count = 0   #counts ensure no time wasted iterating thru whole file, only 1st row analyzed
        data = csv.reader(open(f, 'rU'), quotechar='"', delimiter = ',')
        for line in data:
            count += 1
            if count == 1:
                for string in line:
                    m = re.search('[0-9]{1} [a-zA-Z]+ Deg. C Data', string)
                    found = ""
                    if m:
                        #string is originally in format like "10 Veh Deg. C Data"
                        #.split() splits on space and makes list of string, takes 1st item
                        mouse_ids.add(string)
            elif count <1:
                break   #exits this for loop and goes to first for loop
    #doesn't truly sort, puts '10 nnn' before '2 nnn' for example, but doesn't need to be sorted
    return sorted(mouse_ids)    #note: type is now list

def clean_mouse_ids(mouse_ids):
    """Given a list of mouse IDs in the form of "10 Veh Deg. C Data" (the number and the word
    "Veh" may be different depending on the experiment)."""
    first_items = set([])
    for mouse in mouse_ids:
        first_items.add(mouse.split()[0])
    return sorted(first_items) #returns sorted list (not truly sorted, as '10' is before '2'
        
def get_csv_day_labels(lst_clean_csv_files):
    """Given a list of strings that are the names of the cleaned csv data files, returns a sorted
    list of each date in the experiment (found by looking at the left-most column cells in each
    file given. There are no repeated days in the list: each only occurs once."""
    day_labels = set()   #type is set() to ensure each label only appears once
    for f in lst_clean_csv_files:
        count = 0
        data = csv.reader(open(f, 'rU'), quotechar='"', delimiter = ',')
        for line in data:
            count += 1
            if count > 1:
                day_labels.add(line[0])
    return sorted(day_labels)   #note: type is now list
    

def day_label(day):
    """Given a string that is the filename, isolates the date and returns date in string format.
    NOTE: this assumes the date is before the first comma, after the phrase "clean_"  """
    split_file_name = [x.strip() for x in day.split(',')]
    file_date = split_file_name[0]
    return file_date.split('_')[1]

def sort_day_labels(unsorted_day_labels):
    """Given a list of strings that are day labels (where the first date also has an entry
    that looks like, 'first_day Light Only'), switches the 1st and 2nd item once sorted(). This
    function is needed because sorted() sorts 'Light Only' day as 2nd item, not 1st."""
    first_sort = sorted(unsorted_day_labels)
    template = list(first_sort) #saves a copy of first_sort that is subsequently not modified
    first_sort.remove(first_sort[0]) #removes first full day from its 1st position
    first_sort.insert(1, template[0]) #inserts first full day into 2nd position; Light Only is 1st
    return first_sort

def mouse_label(mouse):
    """Given a string that is the mouse ID in the csv file, isolates the mouse number while
    disregarding the treatment"""
    split_file_name= [x.strip() for x in mouse.split(' ')]
    ID = split_file_name[0]
    return ID
                  
def extract(filename, data_dict, mouse_ids):
    """Given a csv data file and a dictionary template mapping mouse ids to an empty list, this
    function returns a dictionary where the mouse IDs are the keys, mapped to a list of tuples
    that contain (time, temp). An example line would look like:
    {mouse1: [(time1, temp1), (time2, temp2)] , mouse2: [(time1, temp1), (time2, temp2)]...}"""
    
    data = csv.DictReader(open(filename, 'rU'), quotechar='"', delimiter = ',')  
    for row in data:
        for mouse in mouse_ids:     
            if mouse in row:
                if row[mouse] != 'NaN':#cleaning data by ignoring points with no temp
                    if "2 Veh Deg. C Time" in row:   ####make a variable!!!! TEMPORARY
                        pt = (row['2 Veh Deg. C Time'] , float(row[mouse]) )   #makes (time, temp)
                        data_dict[mouse].append((pt))
                    elif "2 Acyline Deg. C Time" in row:
                        pt = (row['2 Acyline Deg. C Time'] , float(row[mouse]) ) #makes (time, temp)
                        data_dict[mouse].append((pt))
    return data_dict

def extract_treatment_start_date(filename):
    """Given a .csv file where the first cell in one of the rows contains the phrase
    'treatment started', returns the string in the cell to the right."""
    data = csv.reader(open(filename, 'rU'), quotechar='"', delimiter = ',')
    for line in data:
        if 'Date treatment started' in line[0]: #need line[0], not line, or code won't find label
            extracted_treatment_start_date = line[1]
            date_pieces = extracted_treatment_start_date.split('/')
            
            #variables good to go IF they are in the correct mm/dd/yyyy
            month = date_pieces[0]
            day = date_pieces[1]
            yr =  date_pieces[2]
            
            #if not in correct format, will at a 0 before m and d, 20 before yy
            if len(date_pieces[0]) < 2:
                month = '0'+date_pieces[0]
            if len(date_pieces[1]) < 2:
                day = '0'+date_pieces[0]
            if len(date_pieces[2]) == 2:
                yr = '20'+date_pieces[2]
                
            date = month + '-' + day + '-' + yr
            return date

def extract_last_n_day_pre(filename):
    """Given a .csv file where the first cell in one of the rows contains the phrase
    'Analyze last n days pre treatment", returns the integer to the right of that cell."""
    data = csv.reader(open(filename, 'rU'), quotechar='"', delimiter = ',')
    for line in data:
        if 'Analyze last n days pre treatment' in line[0]:
            return int(line[1])
        
def extract_last_n_day_post(filename):
    """Given a .csv file where the first cell in one of the rows contains the phrase
    "Analyze last n days post treatment", returns the integer to the right of that cell."""
    data = csv.reader(open(filename, 'rU'), quotechar='"', delimiter = ',')
    for line in data:
        if 'Analyze last n days post treatment' in line[0]:
            return int(line[1])
        
def extract_ints_in_mavg(filename):
    """Given a .csv file where the first cell in one of the rows contains the phrase
    'moving average number of points', returns the integer in the cell to the right."""
    data = csv.reader(open(filename, 'rU'), quotechar='"', delimiter = ',')
    for line in data:
        if 'moving average number of points' in line:
            return int(line[1])

def extract_ints_in_moving_stdev(filename):
    """Given a .csv file where the first cell in one of the rows contains the phrase
    'moving standard deviation number of points', returns the integer in the cell to the right."""
    data = csv.reader(open(filename, 'rU'), quotechar='"', delimiter = ',')
    for line in data:
        if 'moving standard deviation number of points' in line:
            return int(line[1])

def extract_mouse_nums_to_use(filename):
    """Given a .csv file where the first cell in at least one of the rows contains the phrase
    'Treatment ' (NOTE the space after the word 'Treatment'), returns a sorted list of strings of
    all mouse numbers given in those rows. Assumes mouse numbers have no other characters."""
    mouse_nums = []
    data = csv.reader(open(filename, 'rU'), quotechar='"', delimiter = ',')
    for line in data:
        if 'Treatment ' in line[0]:
            for string in line:
                match = re.search('^\d+', string)
                #matches strings that begin with a digit
                found =""
                if match:
                    found = match.group()
                    mouse_nums.append(string)
    #sorts mouse nums (assumes mouse num is ONLY digits)
    return sorted(mouse_nums, key=lambda x:float(x))

def extract_tx1_mice(filename):
    """Given a .csv file where the first cell in one of the rows contains the phrase
    'Treatment 1' (NOTE the space after the word 'Treatment'), returns a sorted list of strings of
    all mouse numbers given in that row. Assumes mouse numbers have no other characters."""
    mouse_nums = []
    data = csv.reader(open(filename, 'rU'), quotechar='"', delimiter = ',')
    for line in data:
        if 'Treatment 1' in line[0]:
            for string in line:
                match = re.search('^\d+', string)
                #matches strings that begin with a digit
                found =""
                if match:
                    found = match.group()
                    mouse_nums.append(string)
    #sorts mouse nums (assumes mouse num is ONLY digits)
    return sorted(mouse_nums, key=lambda x:float(x))

def extract_tx2_mice(filename):
    """Given a .csv file where the first cell in one of the rows contains the phrase
    'Treatment 1' (NOTE the space after the word 'Treatment'), returns a sorted list of strings of
    all mouse numbers given in that row. Assumes mouse numbers have no other characters."""
    mouse_nums = []
    data = csv.reader(open(filename, 'rU'), quotechar='"', delimiter = ',')
    for line in data:
        if 'Treatment 2' in line[0]:
            for string in line:
                match = re.search('^\d+', string)
                #matches strings that begin with a digit
                found =""
                if match:
                    found = match.group()
                    mouse_nums.append(string)
    #sorts mouse nums (assumes mouse num is ONLY digits)
    return sorted(mouse_nums, key=lambda x:float(x))

def separate_light_dark(data_dict, mouse):
    """Given data_dict and a given mouse, will return a dictionary of two lists of (time, temp)
    tuples, one list for the dark cycle (defined as 18:00:00 to 5:59:59) and one list for the
    light cycle (defined as 6:00:00 to 17:59:59).
    Dict looks like {"Light Cycle": [(t0,temp0), (t1,temp(1)...], "Dark Cycle": [(t0, temp0)...]}"""
    parsed_data = {}
    light_data = []
    dark_data = []
    for time_temp_tuple in data_dict[mouse]:
        if '10:00:00'<= time_temp_tuple[0] <= '17:59:59'or time_temp_tuple[0][0]>= '6':
            light_data.append(time_temp_tuple) 
        else:
            dark_data.append(time_temp_tuple)
    parsed_data["Light Cycle"] = light_data
    parsed_data["Dark Cycle"] = dark_data
    return parsed_data

def make_master_tt_dic(filenames, mouse_ids):
    """Makes a dictionary of each day label mapped to each mouse mapped to a dictionary where
    Light/Dark cycle are the keys mapping to a list of time/temp tuples."""
    master_tt_dic = {}
    for day in filenames:
        data_dict = {}
        organized_data = {}
        for mouse in mouse_ids:
            data_dict[mouse]= []  #Makes empty list for mouse if there's no data for mouse that day
        extract(day, data_dict, mouse_ids)   #fleshes out data_dict
        for mouse in mouse_ids:
            num_mouse = mouse.split()[0]
            if len(data_dict[mouse]) > 0:   #ensures mice with no data don't get included
##                if mouse_label(mouse) != '8':
##                    if mouse_label(mouse) != '19':
##                        if mouse_label(mouse) != '20':
##                            organized_data[mouse_label(mouse)] = separate_light_dark(data_dict, mouse)
                organized_data[num_mouse] = separate_light_dark(data_dict, mouse)
        master_tt_dic[day_label(day)] = organized_data
    return master_tt_dic

def extract_all_data(filename, data_dict, mouse_ids):
    """Just like the function "extract", except this function does not filter out 'NaN'."""
    data = csv.DictReader(open(filename, 'rU'), quotechar='"', delimiter = ',')  
    for row in data:
        for mouse in mouse_ids:     
            if mouse in row:
                if "2 Veh Deg. C Time" in row:
                    pt = (row['2 Veh Deg. C Time'] , float(row[mouse]) )   #makes (time, temp)
                    data_dict[mouse].append((pt))
                elif "2 Acyline Deg. C Time" in row:
                    pt = (row['2 Acyline Deg. C Time'] , float(row[mouse]) ) #makes (time, temp)
                    data_dict[mouse].append((pt))
    return data_dict



#delete below
##def extract(filename, data_dict, mouse_ids):
##    """Given a csv data file and a dictionary template mapping mouse ids to an empty list, this
##    function returns a dictionary where the mouse IDs are the keys, mapped to a list of tuples
##    that contain (time, temp). An example line would look like:
##    {mouse1: [(time1, temp1), (time2, temp2)] , mouse2: [(time1, temp1), (time2, temp2)]...}"""
##    
##    data = csv.DictReader(open(filename, 'rU'), quotechar='"', delimiter = ',')  
##    for row in data:
##        for mouse in mouse_ids:     
##            if mouse in row:
##                if row[mouse] != 'NaN':#cleaning data by ignoring points with no temp
##                    if "2 Veh Deg. C Time" in row:   ####make a variable!!!! TEMPORARY
##                        pt = (row['2 Veh Deg. C Time'] , float(row[mouse]) )   #makes (time, temp)
##                        data_dict[mouse].append((pt))
##                    elif "2 Acyline Deg. C Time" in row:
##                        pt = (row['2 Acyline Deg. C Time'] , float(row[mouse]) ) #makes (time, temp)
##                        data_dict[mouse].append((pt))
##    return data_dict
#delete above



def make_all_times_dic(filenames, mouse_ids):
    """Just like 'make_master_tt_dic', but this function includes 'NaN' in tt tuple."""
    all_times_dic = {}
    for day in filenames:
        data_dict = {}
        organized_data = {}
        for mouse in mouse_ids:
            data_dict[mouse]= []  #Makes empty list for mouse if there's no data for mouse that day
        extract_all_data(day, data_dict, mouse_ids)
        for mouse in mouse_ids:
            if len(data_dict[mouse]) > 0:   #ensures mice with no data don't get included
##                if mouse_label(mouse) != '8':
##                    if mouse_label(mouse) != '19':
##                        if mouse_label(mouse) != '20':
##                            organized_data[mouse_label(mouse)] = separate_light_dark(data_dict, mouse)
                organized_data[mouse_label(mouse)] = separate_light_dark(data_dict, mouse)
        all_times_dic[day_label(day)] = organized_data
    return all_times_dic




#delete below
##def make_master_tt_dic(filenames, mouse_ids):
##    """Makes a dictionary of each day label mapped to each mouse mapped to a dictionary where
##    Light/Dark cycle are the keys mapping to a list of time/temp tuples."""
##    master_tt_dic = {}
##    for day in filenames:
##        data_dict = {}
##        organized_data = {}
##        for mouse in mouse_ids:
##            data_dict[mouse]= []  #Makes empty list for mouse if there's no data for mouse that day
##        extract(day, data_dict, mouse_ids)   #fleshes out data_dict
##        for mouse in mouse_ids:
##            num_mouse = mouse.split()[0]
##            if len(data_dict[mouse]) > 0:   #ensures mice with no data don't get included
####                if mouse_label(mouse) != '8':
####                    if mouse_label(mouse) != '19':
####                        if mouse_label(mouse) != '20':
####                            organized_data[mouse_label(mouse)] = separate_light_dark(data_dict, mouse)
##                organized_data[num_mouse] = separate_light_dark(data_dict, mouse)
##        master_tt_dic[day_label(day)] = organized_data
##    return master_tt_dic
#delete above




def extract_light_cycle_times(filename):
    """Given a .csv file where the first cell in one of the rows contains the phrase
    'Light Cycle', and the first cell of another row contains the phrase, "Dark Cycle",
    returns a dictionary where the cycle is mapped to a list that has the string of the
    lower bound of the cycle and the upper bound of the cycle."""
    raw_cycle_bounds = {}
    data = csv.reader(open(filename, 'rU'), quotechar='"', delimiter = ',')
    for line in data:
        if 'Light Cycle' in line:
            raw_cycle_bounds['Light Cycle'] = [line[1], line[2]]
        elif 'Dark Cycle' in line:
            raw_cycle_bounds['Dark Cycle'] = [line[1], line[2]]
    return raw_cycle_bounds

def extract_txt_mouse_ids(filenames):
    """Returns a list of strings that are mouse_ids from a list of strings that are appropriate
    files from the directory the code is run in.
    Takes the first two characters of the third item in the title (items separated by a space)
    and appends that string to a list of mouse_ids."""
    mouse_ids = []
    for data_file in filenames:
        if 'CBT ' in data_file:
            if '.TXT' in data_file:
                mouse_id = extract_one_txt_mouse_id(data_file)
                if mouse_id[0] == '0':
                    mouse_ids.append(mouse_id[1])  ##NOT flexible, this is temporary work around
                else:
                    mouse_ids.append(mouse_id)
    return mouse_ids

def extract_one_txt_mouse_id(filename):
    """Given a string that is a filename (.txt format with mouse ID as the fourth string when split
    on a space), returns a string that is the mouse ID number."""
    split_mouse_id = filename.split(' ')
    if split_mouse_id[3][0] == '0': 
        return split_mouse_id[3][1:2] ##NOT flexible, this is temporary
    else:
        return split_mouse_id[3][0:2] ##NOT flexible, this is temporary

def extract_raw_data_txt(files):
    raw_tt_dict = {}
    for data_file in files:
        if 'CBT ' in data_file:
            mouse_id = extract_one_txt_mouse_id(data_file)
            opened_file = open(data_file, 'r')        
            data = opened_file.readlines()
            one_mouse_dict = {}
            for line in data:
                match = re.search(r'^\s\d\d\/\d\d\/\d\d\d\d\s', line)
                # ^ is new line
                # \s is white space
                # \d is any number 0-9
                if match:
                    splt_ln = line.split(' ')
                    raw_date = str(line.split(' ')[1])
                    date_pieces = raw_date.split('/')
                    date = date_pieces[0]+'-'+date_pieces[1]+'-'+date_pieces[2]
                    #date reformatting needed to prevent errors in plotting code interpreting
                    # / as making new directory
                    time = splt_ln[3]+':00' #needed because original code needs seconds
                    temp = float(splt_ln[4][0:4])
                    tt_tuple = (time, temp ) #time, float(temp)
                    one_mouse_dict.setdefault(date, []).append(tt_tuple)
            for day in one_mouse_dict:
                raw_tt_dict.setdefault(day, {})[mouse_id] = one_mouse_dict[day]
    return raw_tt_dict

def cycle_bounds_txt(raw_cycle_bounds):
    """Given a dictionary containing {'Light Cycle': [early time, later time]}, returns a
    dictionary where times are in format nn:nn, not nn:nn:nn or n:nn etc."""
    cycle_bounds = {}
    for cycle in raw_cycle_bounds:
        for time in raw_cycle_bounds[cycle]:
            reformatting_time = time.split(":")
            if len(reformatting_time[0]) < 2:
                reformatting_time[0] = '0'+reformatting_time[0]
                #concatination needed, or .csv reads '6:00' while .txt reads '06:00'...MUST be same
            usable_time = str(reformatting_time[0]+':'+reformatting_time[1])
            cycle_bounds.setdefault(cycle, []).append(usable_time)
    return cycle_bounds

def get_cycle_bounds_txt(user_input):
    """Given a .csv file where the first cell in one of the rows contains the phrase
    'Light Cycle', and the first cell of another row contains the phrase, "Dark Cycle",
    returns a dictionary where the cycle is mapped to a list that has the string of the lower
    bound of the cycle and the upper bound of the cycle.Formats the times to be nn:nn"""
    raw_cycle_bounds = extract_light_cycle_times(user_input)
    cycle_bounds = cycle_bounds_txt(raw_cycle_bounds)
    return cycle_bounds

def get_last_n_pre_days(treatment_start_date, day_labels, filename):
    """Given a list of strings that are the days (day_labels), and a string that is the date the
    treatments began (according to tt_dic, NOT actual days since tt_dic days start on the given
    day's dark cycle and end on the end of the next day's light cycyle, making everything slightly
    off), and the name of the .csv file that the user modifies, returns a list of strings that are
    the last n day of the pre treatment. See the file, "read_me_for_user_modify" under the section
    "Analyze last n days" for more information """
    n_days_to_analyze = extract_last_n_day_pre(filename) #this is an int
    print
    print n_days_to_analyze
    print type(n_days_to_analyze)
    count = 0
    for day in day_labels: #makes list of all days between expt. start and treatment start date
        count +=1
        if day == treatment_start_date:
            all_pre_days = day_labels[0:count-1]
            n_pre_days = all_pre_days[-n_days_to_analyze:len(all_pre_days)]
            #makes new list of days counting back from treatment day that's n_days_to_analyze long
            #list does not include treatment_start_date
            return n_pre_days

def get_last_n_post_days(treatment_start_date, day_labels, filename):
    """Given a list of strings that are the days (day_labels), and a string that is the date the
    treatments began (according to tt_dic, NOT actual days since tt_dic days start on the given
    day's dark cycle and end on the end of the next day's light cycyle, making everything slightly
    off), and the name of the .csv file that the user modifies, returns a list of strings that are
    the last n day of the post treatment. See the file, "read_me_for_user_modify" under the section
    "Analyze last n days" for more information """
    n_days_to_analyze = extract_last_n_day_post(filename) #this is an int
    all_post_days = list(day_labels)
    for day in day_labels: #makes list of all days btw expt. start and treatment start date
        if day != treatment_start_date:
            all_post_days.remove(day)
        elif day == treatment_start_date:
            #keep in mind treatments usually start mid-light cycle on the treatment start date
            n_post_days = all_post_days[-n_days_to_analyze-1:len(all_post_days)-1]
            #makes new list of days counting back n_days_to_analyze long from end of expt, NOT
            #including the last day in the expt, since this is only about half a cycle typically,
            #or at least only the light cycle, not the dark cycle
            return n_post_days

def separate_light_dark_txt(data_dict, mouse, user_input):
    """Given data_dict (a dict where keys are mouse IDs mapped to a list of tt tuples), and a
    given mouse, will return a dictionary of two lists of (time, temp) tuples, one list for the
    dark cycle (defined as 18:00:00 to 5:59:59) and one list for the light cycle (defined as
    6:00:00 to 17:59:59).Dict looks like {"Light Cycle": [(t0,temp0), (t1,temp(1)...],
    "Dark Cycle": [(t0, temp0)...]}"""
    parsed_data = {}
    light_data = []
    dark_data = []
    cycle_bounds = get_cycle_bounds_txt(user_input)
    earliest_light = cycle_bounds['Light Cycle'][0]
    latest_light = cycle_bounds['Light Cycle'][1]

    for time_temp_tuple in data_dict[mouse]:
        if earliest_light <= time_temp_tuple[0] <= latest_light:
            light_data.append(time_temp_tuple) 
        else:
            dark_data.append(time_temp_tuple)
    parsed_data["Light Cycle"] = light_data
    parsed_data["Dark Cycle"] = dark_data
    return parsed_data

def make_raw_master_tt_dic_txt(filenames, user_input):
    """Makes a dictionary of each day label mapped to each mouse mapped to a dictionary where
    Light/Dark cycle are the keys mapping to a list of time/temp tuples."""
    raw_master_tt_dic = {}
    raw_tt_dict = extract_raw_data_txt(filenames)
    #raw_tt_dict is day: {mouse: [(t,t)...], mouse2: [(t,t)...]...}, day2: {mouse:[],..}, ...}
    for day in raw_tt_dict:
        for mouse in raw_tt_dict[day]:
            separated_cycles_dict = separate_light_dark_txt(raw_tt_dict[day], mouse, user_input)
            raw_master_tt_dic.setdefault(day, {})[mouse] = separated_cycles_dict
    return raw_master_tt_dic
        
def calibration_data_dict(data_files):
    """Given a list of files in the directory, takes the file with "Calibration Document" in the
    title (must be a .csv for this to work) and return a dictionary where the datalogger id is
    mapped to a list where the first string is slope and the second string is intercept."""
    cal_line = []
    calibration_dict = {}
    for f in data_files:
        if "Calibration Document" in f:
            filename = f
    data = csv.reader(open(filename, 'rU'), quotechar='"', delimiter = ',')
    for line in data:
        for cell in line:
            if '0000' in cell: ##NOT very flexible, change
                if line[1].split('-')[1][0] == '0':
                    datalog_id = line[1].split('-')[1][1] ##NOT very flexible, change
                else:
                    datalog_id = line[1].split('-')[1]
                calibration_dict[datalog_id] = [line[11], line[12]]
                 ##NOT very flexible
                 #line[11] is slope
                 #line[12] is intercept
    return calibration_dict

def calibrate_data(data_files, master_tt_dic):
    """Given a list of files in the directory, finds the Calibration Document and uses that
    information to convert each data logger's temperature using the given constants. Returns a
    dictionary identical to given master_tt_dic (day: {mouse: {cycle: [tt, tt...]}}} but the
    temperature has now been computed with the following equation: (raw_temp - intercept) * slope
    and a new dictionary has been made so the temperature is this new calibrated temperature (the
    float is rounded to the nearest hundredth, though the technology is only accurate to the tenth."""
    calibrated_raw_tt_dict = {}
    calibration_dict = calibration_data_dict(data_files)
    for day in master_tt_dic:
        calibrated_raw_tt_dict[day] = {}
        for mouse in master_tt_dic[day]:
            calibrated_raw_tt_dict[day][mouse] = {}
            slope = float(calibration_dict[mouse][0])
            intercept = float(calibration_dict[mouse][1])
            for cycle in master_tt_dic[day][mouse]:
                calibrated_raw_tt_dict[day][mouse][cycle] = []
                for tt in master_tt_dic[day][mouse][cycle]:
                    calibrated_temp = (float(tt[1]) - intercept)*slope
                    rounded_cal_temp = round(calibrated_temp, 2)
                    hundredths_cal_temp = '%.2f' % rounded_cal_temp
                    #rounds the calibrated_temp to the nearest hundredth and stores it as nn.nn
                    time_temp = (tt[0], float(hundredths_cal_temp))
                    calibrated_raw_tt_dict[day][mouse][cycle].append(time_temp)
                #as given by manufacturer: Corrected Temp = (Measured temp - intercept)*(Slope)
    return calibrated_raw_tt_dict

##def refit_to_master_tt_dic(calibrated_tt_dic, mouse_ids, user_input): 
##    """Takes a tt_dic still mapped to its original days (what appears in the raw .txt data) and
##    shifts the data to conform to the .csv way of organizing data, where the light cycle now belongs
##    to the previous day's entry, and the dark cycle between midnight and morning belongs to the
##    previous day as well. The only mapping that stays the same is the dark cycle between evening and
##    midnight. This function assumes that the first day begins recording in the light cycle and the
##    last day stops recording in the light cycle."""
##    master_tt_dic = {}
##
##    cycle_bounds = get_cycle_bounds_txt(user_input)           
##    earliest_light = cycle_bounds['Light Cycle'][0]
##    latest_light = cycle_bounds['Light Cycle'][1]
##    
##    #determines what days the raw data recorded, and makes new list of days for master_tt_dic to use
##    all_days = sorted(calibrated_tt_dic.keys())
##    
##    #building master_tt_dic structure
##    #this is preferable to setdefault due to different day keys btw the two dictionaries
##    ##CHECK IF THIS IS STILL NECESSARY (temporary)
##    for day in all_days:
##        master_tt_dic[day] = {}
##        for mouse in mouse_ids:
##            master_tt_dic[day][mouse]={}
##            master_tt_dic[day][mouse]['Dark Cycle'] = [] 
##            master_tt_dic[day][mouse]['Light Cycle'] = []                
##    #sort into correct bins        
##    count=0
##    for day in all_days:
##        for mouse in calibrated_tt_dic[day]:
##            for cycle in calibrated_tt_dic[day][mouse]:
##                for tt in calibrated_tt_dic[day][mouse][cycle]:
##                    if earliest_light <= tt[0] <= latest_light:
##                        master_tt_dic[day][mouse]['Light Cycle'].append(tt)
##                    elif latest_light < tt[0] <= '23:59':
##                        master_tt_dic[day][mouse]['Dark Cycle'].append(tt)
##                    elif '00:00' <= tt[0] < earliest_light:
##                        #this block needed, or count-1 adds 1st day's data to last day
##                        if count == 0:
##                            master_tt_dic[all_days[count]][mouse]['Dark Cycle'].append(tt)
##                        elif count > 0:
##                            master_tt_dic[all_days[count-1]][mouse]['Dark Cycle'].append(tt)
##    return master_tt_dic

def refit_to_master_tt_dic(calibrated_tt_dic, mouse_ids, user_input): 
    """Takes a tt_dic still mapped to its original days (what appears in the raw .txt data) and
    shifts the data to conform to the .csv way of organizing data, where the light cycle now belongs
    to the previous day's entry, and the dark cycle between midnight and morning belongs to the
    previous day as well. The only mapping that stays the same is the dark cycle between evening and
    midnight. This function assumes that the first day begins recording in the light cycle and the
    last day stops recording in the light cycle."""
    master_tt_dic = {}

    cycle_bounds = get_cycle_bounds_txt(user_input)           
    earliest_light = cycle_bounds['Light Cycle'][0]
    latest_light = cycle_bounds['Light Cycle'][1]
    
    #determines what days the raw data recorded, and makes new list of days for master_tt_dic to use
    all_days = sorted(calibrated_tt_dic.keys())
    
    #building master_tt_dic structure
    #this is preferable to setdefault due to different day keys btw the two dictionaries
    ##CHECK IF THIS IS STILL NECESSARY (temporary)
    for day in all_days:
        master_tt_dic[day] = {}
        for mouse in mouse_ids:
            master_tt_dic[day][mouse]={}
            master_tt_dic[day][mouse]['Dark Cycle'] = [] 
            master_tt_dic[day][mouse]['Light Cycle'] = []                
    #sort into correct bins        
    count=0
    for day in all_days:
        count += 1
        for mouse in calibrated_tt_dic[day]:
            latest_premidnight = []
            for cycle in calibrated_tt_dic[day][mouse]:
                for tt in calibrated_tt_dic[day][mouse][cycle]:
                    if earliest_light <= tt[0] <= latest_light:
                        master_tt_dic[day][mouse]['Light Cycle'].append(tt)
                    elif latest_light < tt[0] <= '23:59':
                        latest_premidnight.append(tt)
 ##                       master_tt_dic[day][mouse]['Dark Cycle'].append(tt)
                    elif '00:00' <= tt[0] < earliest_light:
                        if count < 2:
                            master_tt_dic[all_days[count-1]][mouse]['Dark Cycle'].append(tt)
                        elif count >= 2:
                            master_tt_dic[all_days[count-2]][mouse]['Dark Cycle'].append(tt)
            t_cnt = 0
            for time_temp in latest_premidnight:
                t_cnt +=1
                master_tt_dic[day][mouse]['Dark Cycle'].insert(t_cnt-1, time_temp)
                        
                        #this block needed, or count-1 adds 1st day's data to last day
##                        if count == 0:
##                            master_tt_dic[all_days[count]][mouse]['Dark Cycle'].append(tt)
##                        elif count > 0:
##                            master_tt_dic[all_days[count-1]][mouse]['Dark Cycle'].append(tt)
    return master_tt_dic     

def list_CBT(day, mouse, cycle, master_tt_dic):
    """Returns a list of CBTs for the given day, mouse and light cycle"""
    CBT_list = []
    for time_temp_tuple in master_tt_dic[day][mouse][cycle]:
        CBT_list.append(time_temp_tuple[1])
    return CBT_list

def list_times(day, mouse, cycle, master_tt_dic):
    """Returns a list of times in seconds) for the given day, mouse, and light cycle"""
    time_list = []
    for time_temp_tuple in master_tt_dic[day][mouse][cycle]:
        time_list.append(hms_to_secs(time_temp_tuple[0]) / 3600.0) #3600.0 is secs/hr
    return time_list

def stder(CBT_list):
    """Returns the standard error of the mean for the list of CBTs"""
    stdev_temp = np.std(CBT_list)
    stder_temp = stdev_temp / math.sqrt(len(CBT_list))
    return stder_temp

def stdev(CBT_list):
    """Returns the population standard deviation of the list of CBTs"""
    stdev_temp = np.std(CBT_list)
    return stdev_temp

def find_all_avgs_ers(day_labels, mouse_nums, times, master_tt_dic): #Ultimately into excel, not print
    """Prints the mean, standard error and standard deviation for each mouse for each day's light
    cycle. If there is no data for the given day and cycle, the function prints 'There is no data
    for this cycle' """
    file = open("mean_stder_stdev.txt", "w")
    for day in day_labels:
        file.write("\n")
        file.write(day + "\n")
        for mouse in mouse_nums:
            file.write("\n")
            file.write("Mouse" + mouse + "\n")
            for cycle in times:
                file.write(cycle + "\n")
                if len(master_tt_dic[day][mouse][cycle]) > 0: #needed for files w/o a cycle
                    CBT_lst = list_CBT(day, mouse, cycle, master_tt_dic)
                    mean = str(np.mean(CBT_lst))
                    std_er = str(stder(CBT_lst))
                    std_dev = str(stdev(CBT_lst))
                    file.write("Mean:" + mean + "\n")
                    file.write("Population STD error:" + std_er + "\n")
                    file.write("Population STD deviation:"+ std_dev + "\n")
                else:
                    file.write("There is no data for this cycle\n")
    file.close()

def n_pt_mavg(CBT_list, n_ints_in_mavg):
    """Given a list of CBTs (floats) and an integer, returns a new list of averaged pts (n point
    moving averages). First and last edge points simply have fewer points averaged.
    NOTE- even numbered n_ints_in_mavg will take (n_ints_in_mavg/2) pts before and after the
    selected point."""
    left_norm = n_ints_in_mavg/2 + 1
    right_norm = n_ints_in_mavg/2
    new_list = []
    count = 0
    for temp in CBT_list:
        count += 1
        if count <= n_ints_in_mavg/2:
            new_list.append(np.mean(CBT_list[count-count:count+right_norm]))
        elif count >= len(CBT_list)-right_norm:
            new_list.append(np.mean(CBT_list[count-left_norm:len(CBT_list)]))
        else:
            new_list.append(np.mean(CBT_list[count-left_norm:count+right_norm]))
    return new_list

def n_moving_stdev(CBT_list, n_stdev):
    """Given a list of CBTs (floats) and an integer, returns a new list of standard deviation of
    selected pts (number of points used is n_stdev). First and last edge points simply have fewer points averaged.
    NOTE- even numbered n_stdev will take (n_stdev/2) pts before and after the selected point."""
    left_norm = n_stdev/2 + 1
    right_norm = n_stdev/2
    new_list = []
    count = 0
    for temp in CBT_list:
        count += 1
        if count <= n_stdev/2:
            new_list.append(scipy.stats.tstd(CBT_list[count-count:count+right_norm]))
        elif count >= len(CBT_list)-right_norm:
            new_list.append(scipy.stats.tstd(CBT_list[count-left_norm:len(CBT_list)]))
        else:
            new_list.append(scipy.stats.tstd(CBT_list[count-left_norm:count+right_norm]))
    return new_list

def extract_n_moving_stdv_axis(filename):
    """Given a .csv file where the first cell in one of the rows contains the phrase
    'Analyze last n days pre treatment", returns the integer to the right of that cell."""
    data = csv.reader(open(filename, 'rU'), quotechar='"', delimiter = ',')
    for line in data:
        if 'Moving stdev plot y axis range' in line[0]:
            min_bound = float(line[1])
            max_bound = float(line[2])
            return [min_bound, max_bound]
        
def plot_n_moving_stdv(day_list, mouse_list, cycle_list, master_tt_dic, n_stdev, filename):
    """Given a list of days (strings), mouse numbers (strings), cycles (strings), the master_tt_dic,
    and the number of points to be used in calculating the standard deviation, saves a plot
    of time of day vs.standard deviation of range given to analyze around that time. Saves one plot
    per day per mouse."""
    ylims = extract_n_moving_stdv_axis(filename)
    for day in day_list:
        for mouse in mouse_list:
            CBT_list = []
            time_list = []
            for cycle in cycle_list:
                CBT_list_first = list_CBT(day, mouse, cycle, master_tt_dic)
                time_list_first = list_times(day, mouse, cycle, master_tt_dic)
                #this combines dark and light cycle for each day
                CBT_list.extend(CBT_list_first)
                time_list.extend(time_list_first)
            stdev_lst = n_moving_stdev(CBT_list, n_stdev)
                
            #assert len(stdev_lst) == len(CBT_list)
            #assert len(stdev_lst) == len(time_list)
            
            fig = plt.figure()
            ax = fig.gca()
            
            plt.figure()
            plt.plot(time_list, stdev_lst, 'b.-')
            plt.title("%s pt. Moving Standard Deviation %s, mouse %s" %(str(n_stdev),day, mouse))
            plt.ylabel("%s point sample std. deviation" %(str(n_stdev))) 
            plt.xlabel("Time of Day (hrs)")
            plt.xlim(0, 24)
            plt.ylim(0, 0.6, 0.1)
            ax.set_xticks(np.arange(0,24,1)) 
            ax.set_yticks(np.arange(0,0.6, 0.1))  
            plt.grid()
            plt.savefig(str(n_stdev)+"_moving stdv graphs/%s_pt_stdv_plot_%s_mouse_%s.png" %(str(n_stdev),day, mouse))
            plt.close()
    
def make_mav_master_dic(day_labels, mouse_nums, times, master_tt_dic, n_ints_in_mavg):
    mav_master_dic = {}
    for day in day_labels:
        mav_master_dic[day] = {}
        for mouse in mouse_nums:
            mav_master_dic[day][mouse] = {}
            for cycle in times:
                mav_master_dic[day][mouse][cycle] = n_pt_mavg(list_CBT(day, mouse, cycle, master_tt_dic), n_ints_in_mavg)
    return mav_master_dic
                                                                
def hms_to_secs(t):
    """Given a time (ie '18:00:28'), returns that time in seconds as an int (ie 64828)"""
    h, m, s = [int(i) for i in t.split(':')]
    return 3600*h + 60*m + s
    
def get_x_axis(tt_dic, mouse, pre_post):
    """Returns a list of lists of times or CBTs per mouse per cycle"""
    x_data = []
    for tt in tt_dic[pre_post][mouse]:
        x_data.append(hms_to_secs(tt[0]))
    return x_data

def get_y_axis(tt_dic, mouse, pre_post):
    """Returns a list of lists of times or CBTs per mouse per cycle"""
    y_data = []
    for tt in tt_dic[pre_post][mouse]:
        y_data.append(tt[1])
    return y_data
        
def plot_mouse(x_pre, y_pre, x_post, y_post, tx1_tx2, mouse):
    """Saves a plot of CBT vs. time for given mouse in given day and cycle"""
    plt.figure()
    
    plt.subplot(211)
    plt.scatter(x_pre, y_pre)
    plt.ylabel('CBT in deg. C')
    plt.title("Pre")
    plt.xlim(0, 86400)
    plt.ylim(34.5, 39)
    
    plt.subplot(212)
    plt.scatter(x_post, y_post)
    plt.xlabel('Time of Day in seconds')
    plt.title("Post")
    plt.xlim(0, 86400)
    plt.ylim(34.5, 39)
    
    plt.suptitle("Mouse "+ mouse)
    plt.savefig("2_day_pre_post_" + mouse + ".png")


def daily_temps_dic(master_tt_dic):
    """Given master_tt_dic, returns a new dictionary with each day mapped to all the times in that
    file. Each time is then mapped to a list of each mouse's CBT for that day and time."""
    daily_avgs = {}
    for day in master_tt_dic:
        daily_avgs[day] = {}
        for mouse in master_tt_dic[day]:
            for cycle in master_tt_dic[day][mouse]:
                for i in range(len(master_tt_dic[day][mouse][cycle])):
                    #maps each day to an empty list initially
                    daily_avgs[day].setdefault(master_tt_dic[day][mouse][cycle][i][0] , [])
                    daily_avgs[day][master_tt_dic[day][mouse][cycle][i][0]].append(master_tt_dic[day][mouse][cycle][i][1])
    for day in daily_avgs:
        for time in daily_avgs[day]:
            avg = np.mean(daily_avgs[day][time])
            daily_avgs[day][time] = avg
    return daily_avgs

def x_times(daily_avgs, day): #####combine this function with get_x_axis?
    """Returns list of times"""
    x_data = []
    for time in daily_avgs[day]:
        x_data.append(hms_to_secs(time) / 3600.0) #3600.0 is secs/hr
    return x_data

def y_avgs(daily_avgs, day):  #####combine this function with get_x_axis?
    """Returns list of times"""
    y_data = []
    for time in daily_avgs[day]:
        y_data.append(daily_avgs[day][time])
    return y_data

def extract_avg_plot_axis(filename):
    """Given a .csv file where the first cell in one of the rows contains the phrase
    'Analyze last n days pre treatment", returns the integer to the right of that cell."""
    data = csv.reader(open(filename, 'rU'), quotechar='"', delimiter = ',')
    for line in data:
        if 'Avg plot y axis range' in line[0]:
            min_bound = int(line[1])
            max_bound = int(line[2])
            return [min_bound, max_bound]
        
def avg_plot(daily_avgs, day, ylims):
    """Given master_tt_dic and a day, plots an average of all mice CBTs for that day"""
    plt.figure()
    x_data = x_times(daily_avgs, day)
    y_data = y_avgs(daily_avgs, day)

    fig = plt.figure()
    ax = fig.gca()
        
    plt.scatter(x_data, y_data)
    plt.ylabel("CBT in deg. C")
    plt.xlabel("Time of day in hours")
    plt.title("Mouse CBT averaged for " + day)
    plt.xlim(0, 24)
    plt.ylim(ylims[0], ylims[1])
    
    ax.set_xticks(np.arange(0,24,1))
    ax.set_yticks(np.arange(ylims[0], ylims[1], 0.5))
    plt.grid()
    #saves to directory 'avg_plot graphs'
    plt.savefig(os.path.join('avg_plot graphs', day + '_mouse_avgs.png'))
    plt.close()

def make_a_directory(directory_name):
    """Given a string that you want to be the directory name, makes a directory with that name in
    the current directory the code is running in"""
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    
def all_avg_plots(master_tt_dic, filename):
    """Saves all averaged daily plots"""
    daily_avgs = daily_temps_dic(master_tt_dic)
    ylims = extract_avg_plot_axis(filename)
    for day in daily_avgs:
        avg_plot(daily_avgs, day, ylims)
        
####################################################

def get_last_two_cycles_moving_stdev(master_tt_dic, mouse, n_stdev, last_two_cycles): 
    """Prints each cycle's average moving standard deviation for given mouse.
    n_stdev is extracted from user_modify file, dictating the number of points to use in the stdev calculation
    last_two_cycles is a list of the last three days"""
    dark_cycle = []
    light_cycle = []
    dark_cycle.extend(n_moving_stdev(list_CBT(last_two_cycles[1], mouse, 'Dark Cycle', master_tt_dic), n_stdev))
    dark_cycle.extend(n_moving_stdev(list_CBT(last_two_cycles[2], mouse, 'Dark Cycle', master_tt_dic), n_stdev))
    light_cycle.extend(n_moving_stdev(list_CBT(last_two_cycles[0], mouse, 'Dark Cycle', master_tt_dic), n_stdev)) 
    light_cycle.extend(n_moving_stdev(list_CBT(last_two_cycles[1], mouse, 'Dark Cycle', master_tt_dic), n_stdev)) 
    return [str(np.mean(dark_cycle)), str(np.mean(light_cycle))]

def get_all_last_2_cycles_moving_stdev(master_tt_dic, tx1_mice, tx2_mice, n_stdev, last_two_cycles):
    """Prints each cycle's average moving standard deviation for all mice, also prints the
    treatment group the mouse belonged to"""
    file = open("all_last_2_cycles_moving_stdev.txt", "w")
    
    file.write("Treatment 1 Mice-Avg of last two full days' %s pt. moving stdev" %(str(n_stdev)) + '\n')
    for mouse in tx1_mice:
        file.write("Mouse "+ mouse + '\n')
        mean_lst = get_last_two_cycles_moving_stdev(master_tt_dic, mouse, n_stdev, last_two_cycles)
        file.write("Dark Cycle: " + mean_lst[0] + '\n')
        file.write("Light Cycle: " + mean_lst[1] + '\n')
        file.write('\n')
    file.write('\n')
    file.write("Treatment 2 Mice-Avg of last two full days' %s pt. moving stdev" %str((n_stdev)) + '\n')
    for mouse in tx2_mice:
        file.write("Mouse" + mouse + '\n')
        tx2_mean_lst = get_last_two_cycles_moving_stdev(master_tt_dic, mouse, n_stdev, last_two_cycles)
        file.write("Dark Cycle: " + tx2_mean_lst[0] + '\n')
        file.write("Light Cycle: " + tx2_mean_lst[1] + '\n')
        file.write('\n') 
        
def make_organized_time_temp_list(last_four_days, times, mouse_list, all_times_dic):
    """Given a list of days, cycles, mice, and master/all_times dict, returns a list of
    (time point, avg CBT for given mice) tuples organized by time point."""
    day_to_tt_dict = make_time_to_temps_dict(last_four_days, times, mouse_list, all_times_dic)
    time_pt_tuples_list = avg_each_time_CBT(day_to_tt_dict, last_four_days)
    organized_tts_list = sorted(time_pt_tuples_list)
    return organized_tts_list

def make_time_to_temps_dict(last_four_days, times, mouse_list, all_times_dic):
    """Given a list of days, a list of light cycles, a list of mouse IDs, and all_times or master_tt
    dictionary, returns the following dic---
    day:[(time point, [mouse1temp, mouse2temp...]),(time point2, [mouse1temp, mouse2temp...])...]
    """
    avged_list = []
    pre_time_temps_list = []
    pre_time_temps_d = {}
    day_tally = 0
    
    for day in last_four_days:
        day_tally += 1
        pre_time_to_temps_dic = collections.defaultdict(list)
        for mouse in mouse_list:
            if day_tally == 1:
                count = 0
            if day_tally > 1:
                count = saved_count
            for cycle in times:
##                if day=='8-14-14' and cycle=='Light Cycle':
##                    #8-14 light is first treatment cycle
##                    continue #exits this for loop
##                if day=='8-11-14 Light only' and cycle=='Dark Cycle':
##                    # 8-11 dark doesn't have any data
##                    continue #do not use break, doesn't work with Dark Cycle
##                if day =='8-20-14' and cycle=='Light Cycle':
##                    #8-20 light cycle is incomplete
##                    continue
                for key, value in all_times_dic[day][mouse][cycle]:
                    count+=1  #count allows for easier sorting later
                    pre_time_to_temps_dic[count].append(value)
        saved_count = int(count) #not sure if int() is necessary
        pre_time_temps_d[day] = pre_time_to_temps_dic.items()
    return pre_time_temps_d

def avg_each_time_CBT(pre_time_temps_d, day_list):
    """Given a dictionary of day mapped to list of (time point, [mouse1temp, mouse2temp...])...
    returns a new list of tuples where time pt. is coupled w/the avg CBT of all of selected mice"""
    avg_tt_lst = []
    for day in day_list:
        for t in pre_time_temps_d[day]:
            avg_tt_lst.append( (t[0], np.mean(t[1])) )
    return avg_tt_lst

def parse_list(any_list, sample_frequency):
    """Given a list, returns a new list where every nth tuple from the initial list is included.
    N is determined by the given sample_frequency integer."""
    #NOTE- unless the given list is organized, parsing will be meaningless because it can only
    #parse the list in the order the list is given.
    parsed_list = any_list[::sample_frequency]
    return parsed_list

def plot_each_treatment_last_days(last_four_days_pre, last_four_days_post, times, tx2_mice,
                                  tx1_mice, sample_frequency, all_times_dic):
    """Given lists of pre and post treatment days, a list of light cycles, lists of treatment 1 and
    treatment 2 mice, and all_time or master_tt dictionary, saves two plots. One of pre-treatment
    temperature averages every sample_frequency, and one of post-treatment temperature averages every
    sample_frequency. Each plot has two lines-one of treatment 1 and treatment 2 averages."""
    pre_tx2_lst =  make_organized_time_temp_list(last_four_days_pre, times, tx2_mice, all_times_dic)
    pre_tx1_lst =  make_organized_time_temp_list(last_four_days_pre, times,tx1_mice, all_times_dic)
    post_tx2_lst =  make_organized_time_temp_list(last_four_days_post, times, tx2_mice, all_times_dic)
    post_tx1_lst =  make_organized_time_temp_list(last_four_days_post, times, tx1_mice,all_times_dic)
    #necessary for controling axis
    fig = plt.figure()
    ax = fig.gca()
        
    x_times = parse_list( make_last_days_x_list(pre_tx2_lst), sample_frequency)
    y_tx2 = parse_list( make_last_days_y_list(pre_tx2_lst), sample_frequency)
    y_tx1 = parse_list( make_last_days_y_list(pre_tx1_lst), sample_frequency)
    #blue is tx2, red is tx1
    plt.figure()
    plt.plot(x_times, y_tx2, 'b-', label='Treatment 2')
    plt.plot(x_times, y_tx1, 'r-', label='Treatment 1')
    plt.legend()
    plt.ylim(35, 38.5, .5)
    plt.xticks(np.arange(min(x_times), max(x_times)+120, 1440))
    plt.title('Avg. across treatment every time point-pre')
    plt.xlabel('Time in data points')
    plt.ylabel('Average CBT in deg C')
    plt.savefig('last_pre_days_avg.png')
    #plt.show()

    plt.figure()
    x_times = parse_list( make_last_days_x_list(post_tx2_lst), sample_frequency)
    y_tx2 = parse_list( make_last_days_y_list(post_tx2_lst), sample_frequency)
    y_tx1 = parse_list( make_last_days_y_list(post_tx1_lst), sample_frequency)
    #blue is tx2, red is tx1
    plt.figure()
    plt.plot(x_times, y_tx2, 'b-',label = 'Treatment 2')
    plt.plot(x_times, y_tx1, 'r-', label = 'Treatment 1')
    plt.legend()
    plt.xticks(np.arange(min(x_times), max(x_times)+120, 1440))
    plt.ylim(35, 38.5, 0.5)
    plt.title('Avg. across treatment every 30 seconds-post')
    plt.xlabel('Time (every 1440 is 12 hrs)')
    plt.ylabel('Average CBT in deg C')
    plt.savefig('last_post_days_avg.png')
    #plt.show()
#####################################################################################################

def stdev_each_time_CBT(pre_time_temps_d, day_list):
    """Given a dictionary of day mapped to list of (time point, [mouse1temp, mouse2temp...])...
    returns a new list of tuples where time pt. is coupled w/the standard deviation for the CBT
    of all of selected mice"""
    stdev_tt_lst = []
    for day in day_list:
        for t in pre_time_temps_d[day]:
            stdev_tt_lst.append( (t[0], stdev(t[1])) )
    return stdev_tt_lst 

def make_stdev_organized_time_temp_list(days, times, mouse_list, pre_time_temps_dic):
    """Given a list of days, cycles, mice, and master/all_times dict, returns a list of
    (time point, avg CBT for given mice) tuples organized by time point."""
    day_to_tt_dict = make_time_to_temps_dict(days, times, mouse_list, pre_time_temps_dic)
    time_pt_tuples_list = stdev_each_time_CBT(day_to_tt_dict, days)
    organized_tts_list = sorted(time_pt_tuples_list)
    return organized_tts_list

def plot_stdev_each_treatment_last_days(last_four_days_pre, last_four_days_post, times, tx2_mice,
                                  tx1_mice, sample_frequency, all_times_dic):
    """Given lists of pre and post treatment days, a list of light cycles, lists of tx1 and
    tx2 mice, and all_time or master_tt dictionary, saves two plots. One of pre-treatment
    temperature averages every sample_frequency, and one of post-treatment temperature averages every
    sample_frequency. Each plot has two lines-one of tx1 and tx2 averages."""
    pre_tx2_lst =  make_stdev_organized_time_temp_list(last_four_days_pre, times, tx2_mice, all_times_dic)
    pre_tx1_lst =  make_stdev_organized_time_temp_list(last_four_days_pre, times, tx1_mice, all_times_dic)
    post_tx2_lst =  make_stdev_organized_time_temp_list(last_four_days_post, times, tx2_mice, all_times_dic)
    post_tx1_lst =  make_stdev_organized_time_temp_list(last_four_days_post, times,tx1_mice,all_times_dic)
    #necessary for controling axis
    fig = plt.figure()
    ax = fig.gca()
        
    x_times = parse_list( make_last_days_x_list(pre_tx2_lst), sample_frequency)
    y_tx2 = parse_list( make_last_days_y_list(pre_tx2_lst), sample_frequency)
    y_tx1 = parse_list( make_last_days_y_list(pre_tx1_lst), sample_frequency)
    #blue is tx2, red is tx1
    plt.figure()
    plt.plot(x_times, y_tx2, 'b-', label='Treatment 2')
    plt.plot(x_times, y_tx1, 'r-', label='Treatment 1')
    plt.legend()
    plt.ylim(0, 3, .25)
    plt.xticks(np.arange(min(x_times), max(x_times)+120, 1440))
    plt.title('Stdev across treatment every 30 seconds-pre')
    plt.xlabel('Time (every 1440 is 12 hrs)')
    plt.ylabel('Stdev of CBT in deg C')
    plt.savefig('last_pre_days_stdev.png')
    #plt.show()

    plt.figure()
    x_times = parse_list( make_last_days_x_list(post_tx2_lst), sample_frequency)
    y_tx2 = parse_list( make_last_days_y_list(post_tx2_lst), sample_frequency)
    y_tx1 = parse_list( make_last_days_y_list(post_tx1_lst), sample_frequency)
    #blue is tx2, red is tx1
    plt.figure()
    plt.plot(x_times, y_tx2, 'b-',label = 'Treatment 2')
    plt.plot(x_times, y_tx1, 'r-', label = 'Treatment 1')
    plt.legend()
    plt.xticks(np.arange(min(x_times), max(x_times)+120, 1440))
    plt.ylim(0, 3, 0.25)
    plt.title('Stdev across treatment every 30 seconds-post')
    plt.xlabel('Time (every 1440 is 12 hrs)')
    plt.ylabel('Stdev CBT in deg C')
    plt.savefig('last_post_days_stdev.png')
    #plt.show()
    
####################################################################################################       
def make_last_days_x_list(time_tuple_lst):
    """Given a list of time temp tuples, returns a list of times"""
    time_list = []
    for tpl in time_tuple_lst:
        time_list.append(tpl[0])
    return time_list
    
def make_last_days_y_list(time_tuple_lst):
    """Given a list of time temp tuples, returns a list of temperatures"""
    temp_list = []
    for tpl in time_tuple_lst:
        temp_list.append(tpl[1])  #do not use "extend"...says "numpy.float64 object is not iterable" 
    return temp_list

def overall_expt_plot(day_labels, times, tx2_mice, tx1_mice, sample_frequency, all_times_dic):
    """Given lists of all the days, all the times, all tx1 mice, all tx2 mice, an integer of
    how often to plot the data points, and the all_times_dic, plots the average CBT of each
    treatment at each time point for the entire experiment."""
    tx2_lst =  make_organized_time_temp_list(day_labels, times, tx2_mice, all_times_dic)
    tx1_lst =  make_organized_time_temp_list(day_labels, times,tx1_mice,all_times_dic) 
    plt.figure()
    x_times = parse_list( make_last_days_x_list(tx2_lst), sample_frequency)
    y_tx2 = parse_list( make_last_days_y_list(tx2_lst), sample_frequency)
    y_tx1 = parse_list( make_last_days_y_list(tx1_lst), sample_frequency)
    #blue is tx2, red is tx1
    plt.figure()
    plt.plot(x_times, y_tx2, 'b-',label = 'Treatment 2')
    plt.plot(x_times, y_tx1, 'r-', label = 'Treatment 1')
    plt.legend()
    plt.xticks(np.arange(min(x_times), max(x_times)+120, 2880), rotation=45)
    plt.ylim(35, 38.5, 0.5)
    plt.title('Avg. across treatment every 30 seconds-entire experiment')
    plt.xlabel('Time (every 2880 is 24 hrs)')
    plt.ylabel('Average CBT in deg C')
    plt.savefig('avg_temp_per_pt_entire_expt.png')

def overall_expt_plot_stdev(day_labels, times, tx2_mice, tx1_mice, sample_frequency, all_times_dic):
    """Given lists of all the days, all the times, all treatment 2 mice, all treatment 1 mice, an integer of
    how often to plot the data points, and the all_times_dic, plots the stdev of the CBT of each
    treatment at each time point for the entire experiment."""
    tx2_lst =  make_stdev_organized_time_temp_list(day_labels, times, tx2_mice, all_times_dic)
    tx1_lst =  make_stdev_organized_time_temp_list(day_labels, times,tx1_mice,all_times_dic) 
    plt.figure()
    x_times = parse_list( make_last_days_x_list(tx2_lst), sample_frequency)
    y_tx2 = parse_list( make_last_days_y_list(tx2_lst), sample_frequency)
    y_tx1 = parse_list( make_last_days_y_list(tx1_lst), sample_frequency)
    #blue is tx2, red is tx1
    plt.figure()
    plt.plot(x_times, y_tx2, 'b-',label = 'Treatment 2')
    plt.plot(x_times, y_tx1, 'r-', label = 'Treatment 1')
    plt.legend()
    plt.xticks(np.arange(min(x_times), max(x_times)+120, 2880), rotation=45)
    plt.ylim(0, 2.5, 0.25)
    plt.title('Stdev across treatment every 30 seconds-entire experiment')
    plt.xlabel('Time (every 2880 is 24 hrs)')
    plt.ylabel('Stdev CBT in deg C')
    plt.savefig('stdev_temp_per_pt_entire_expt.png')
    
#################
#### MAIN
#################
    
    
def main():
    """This is the main python code that is run in this program."""
    
    user_input = 'user_modify.csv'
    
    filenames = get_data_file_names()
    # gets all .csv files in the directory without the words 'test' or 'user' in the file name,
    # also returns files with the term "Proper" or "CSV " as long as they don't have 'test'/'user'

    mouse_nums = extract_mouse_nums_to_use(user_input)
    #mouse_nums is a sorted list of strings of all mouse numbers to be used in analysis
    #as dictated by user

    tx1_mice = extract_tx1_mice(user_input)
    tx2_mice = extract_tx2_mice(user_input)
    #list of strings of mouse numbers in the treatment group

    n_ints_in_mavg = extract_ints_in_mavg(user_input)
    #this defines how many points to be used in calculating moving averages
    
    n_stdev = extract_ints_in_moving_stdev(user_input)
    #this defines how many points to be used in calculating moving standard deviation 

    times = ["Dark Cycle", "Light Cycle"]


    ########
    ######## This determines how to get data & certain variables dependant on data format
    ########
    print
    if '.TXT' in filenames[0]:
        print "This program is expecting .txt data"
        print
        print "Analyzing the following files for experiment data:"
        for data_file in filenames:
            print data_file
        mouse_ids = extract_txt_mouse_ids(filenames)
        raw_master_tt_dic = make_raw_master_tt_dic_txt(filenames, user_input)
        calibrated_tt_dict = calibrate_data(filenames, raw_master_tt_dic)
        master_tt_dic = refit_to_master_tt_dic(calibrated_tt_dict, mouse_ids, user_input)
##        day_labels = sort_day_labels(master_tt_dic.keys()) ##change to properly order
        day_labels = sorted(master_tt_dic.keys())
        
        
        ##NOTE that function to create last_four_days has not been created yet, so not plotted
        
    elif '.csv' in filenames[0]:
        print "This program is expecting .csv data"
        print
        print "Analyzing the following files for experiment data:"
        clean_all_csv_files(filenames)
        clean_csv_data_files = get_clean_data_file_names()
        for clean_file in clean_csv_data_files:
            print clean_file
            
        mouse_ids = get_all_mouse_ids_csv(clean_csv_data_files)
        #mouse ids is a list of strings (that are digits) from the csv files with data
        
        raw_days_tt_dic = make_master_tt_dic(clean_csv_data_files, mouse_ids)
        mouse_ids = clean_mouse_ids(mouse_ids)
        master_tt_dic = refit_to_master_tt_dic(raw_days_tt_dic, mouse_ids, user_input)
        day_labels = sorted(master_tt_dic.keys())

    ########
    ######## The rest of the code is not perturbed by different data formats
    ########
    last_two_cycles = day_labels[-3:len(day_labels)] #this makes a list of the last three days
    
    find_all_avgs_ers(day_labels, mouse_nums, times, master_tt_dic) #modify to make excel doc, NOT print
    mav_master_dic = make_mav_master_dic(day_labels, mouse_nums, times, master_tt_dic, n_ints_in_mavg)
    
    #makes avg_plot graphs directory and all_avg_plots places generated graphs in there
    make_a_directory('avg_plot graphs')
    all_avg_plots(master_tt_dic, user_input)

    #makes n_moving_stdev graphs directory and plot_n_moving_stdv places generated graphs in there
    make_a_directory(str(n_stdev)+'_moving stdv graphs')
    plot_n_moving_stdv(day_labels, mouse_nums, times, master_tt_dic, n_stdev, user_input)
    
    get_all_last_2_cycles_moving_stdev(master_tt_dic, tx1_mice, tx2_mice, n_stdev, last_two_cycles)
    
    ##################################################################################################
    #broken beyond here
    all_times_dic = make_all_times_dic(filenames, mouse_nums)
    overall_expt_plot(day_labels, times, tx2_mice, tx1_mice, 1, all_times_dic) #modify for flexibility!!
    overall_expt_plot_stdev(day_labels, times, tx2_mice, tx1_mice, 1, all_times_dic) #modify for flexibility!!
##
##
##
    tx_start_date = extract_treatment_start_date(user_input)
    #execute the following if given a treatment start date
    if len(tx_start_date) > 0:
        last_n_pre_days = get_last_n_pre_days(tx_start_date, day_labels, user_input)
        last_n_post_days = get_last_n_post_days(tx_start_date, day_labels, user_input)
        plot_each_treatment_last_days(last_n_pre_days, last_n_post_days, times, tx2_mice,
                                      tx1_mice, 1, all_times_dic)

        plot_each_treatment_last_days(last_four_pre_days, last_four_post_days, times, tx2_mice,
                                     tx1_mice, 1, all_times_dic)
        plot_stdev_each_treatment_last_days(last_four_pre_days, last_four_post_days, times, tx2_mice,
                                            tx1_mice, 1, all_times_dic)
    #Make sample frequency a variable
    #Make another function that does moving avg/stdev
    
if __name__ == "__main__":
    main()
    print
    print "All done! Your data awaits you."

##pr.disable()
##s = StringIO.StringIO()
##sortby = 'percall'
##ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
##ps.print_stats()
##print s.getvalue()
