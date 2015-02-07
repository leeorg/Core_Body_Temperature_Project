# Mouse CBT analysis
#
# Lee Organick 
# Steiner Lab
#
#


##import cProfile, pstats, StringIO
##pr = cProfile.Profile()
##pr.enable()

import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import matplotlib.dates as md
import dateutil
import math
pi = math.pi
import scipy
from scipy import stats
import copy
import pylab as pl
import matplotlib.patches as mpatches
import collections


def day_label(day):
    """Given a string that is the filename, isolates the date and returns date in string format"""
    split_file_name = [x.strip() for x in day.split(',')]
    file_date = split_file_name[0]
    return file_date

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
                    if "2 Veh Deg. C Time" in row:
                        pt = (row['2 Veh Deg. C Time'] , float(row[mouse]) )   #makes (time, temp)
                        data_dict[mouse].append((pt))
                    elif "2 Acyline Deg. C Time" in row:
                        pt = (row['2 Acyline Deg. C Time'] , float(row[mouse]) ) #makes (time, temp)
                        data_dict[mouse].append((pt))
    return data_dict

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
        extract(day, data_dict, mouse_ids)
        for mouse in mouse_ids:
            if len(data_dict[mouse]) > 0:   #ensures mice with no data don't get included
                if mouse_label(mouse) != '8':
                    if mouse_label(mouse) != '19':
                        if mouse_label(mouse) != '20':
                            organized_data[mouse_label(mouse)] = separate_light_dark(data_dict, mouse)
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
                if mouse_label(mouse) != '8':
                    if mouse_label(mouse) != '19':
                        if mouse_label(mouse) != '20':
                            organized_data[mouse_label(mouse)] = separate_light_dark(data_dict, mouse)
        all_times_dic[day_label(day)] = organized_data
    return all_times_dic

def list_CBT(day, mouse, cycle, master_tt_dic):
    """Returns a list of CBTs for the given day, mouse and light cycle"""
    CBT_list = []
    for time_temp_tuple in master_tt_dic[day][mouse][cycle]:
        CBT_list.append(time_temp_tuple[1])
    return CBT_list

def list_times(day, mouse, cycle, master_tt_dic):
    """Returns a list of times (sec / 3600.0) for the given day, mouse, and light cycle"""
    time_list = []
    for time_temp_tuple in master_tt_dic[day][mouse][cycle]:
        time_list.append(hms_to_secs(time_temp_tuple[0]) / 3600.0)
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
    for day in day_labels:
        print
        print day
        for mouse in mouse_nums:
            print
            print "Mouse", mouse
            for cycle in times:
                print cycle
                if len(master_tt_dic[day][mouse][cycle]) > 0: #needed for files w/o a cycle
                    CBT_lst = list_CBT(day, mouse, cycle, master_tt_dic)
                    print "Mean:",np.mean(CBT_lst) , \
                          " STD error:", stder(CBT_lst) , \
                          " STD dev:", stdev(CBT_lst)
                else:
                    print "There is no data for this cycle"

##def mav(CBT_list):
##    """Given a list of CBTs, returns a new list of averaged points (9 point moving averages). First
##    and last four points simply have fewer points averaged."""
##    count = 0
##    new = []
##    for temp in CBT_list:
##        count += 1
##        if count == 1:
##            new.append(np.mean(CBT_list[count-1:count+4]))
##        elif count == 2:
##            new.append(np.mean(CBT_list[count-2:count+4]))
##        elif count == 3:
##            new.append(np.mean(CBT_list[count-3:count+4]))
##        elif count == 4:
##            new.append(np.mean(CBT_list[count-4:count+4]))
##            
##        elif count == len(CBT_list):
##            new.append(np.mean(CBT_list[count-5:count]))
##        elif count == len(CBT_list) - 1:
##            new.append(np.mean(CBT_list[count-5:count+1]))
##        elif count == len(CBT_list) - 2:
##            new.append(np.mean(CBT_list[count-5:count+2]))
##        elif count == len(CBT_list) - 3:
##            new.append(np.mean(CBT_list[count-5:count+3]))
##        else:
##            new.append(np.mean(CBT_list[count-5:count+4]))
##    return new

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
    selected pts (number of points used is n_stdev). First and last ten points simply have fewer points averaged.
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

def plot_n_moving_stdv(day_list, mouse_list, cycle_list, master_tt_dic, n_stdev):
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
            plt.savefig("%s_pt_stdv_plot_%s_mouse_%s.png" %(str(n_stdev),day, mouse))
            
    
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
        
def plot_mouse(x_pre, y_pre, x_post, y_post, acyline_veh, mouse):
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
        x_data.append(hms_to_secs(time) / 3600.0)
    return x_data

def y_avgs(daily_avgs, day):  #####combine this function with get_x_axis?
    """Returns list of times"""
    y_data = []
    for time in daily_avgs[day]:
        y_data.append(daily_avgs[day][time])
    return y_data
                   
def avg_plot(daily_avgs, day):
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
    plt.ylim(34.5, 39)
    
    ax.set_xticks(np.arange(0,24,1))
    ax.set_yticks(np.arange(34.5, 39, 0.5))
    plt.grid()
    
    plt.savefig(day + "_mouse_avgs.png" )   
    
def all_avg_plots(master_tt_dic):
    """Saves all averaged daily plots"""
    daily_avgs = daily_temps_dic(master_tt_dic)
    for day in daily_avgs:
        avg_plot(daily_avgs, day)
        
####################################################
        
def write_dic_to_csv(dic_to_map, first_label, second_label, third_label):
    """Given a dictionary, writes keys down one column maps the value in the cell to the right"""
    with open('%s_%s_%s.csv' %(first_label,second_label,third_label), 'wb') as f:
         w = csv.writer(f)
         w.writerows(dic_to_map.items())

def get_last_two_cycles_moving_stdev(master_tt_dic, mouse, n_stdev): ##modify for future use (dates)
    """Prints each cycle's average moving standard deviation for given mouse"""
    dark_cycle = []
    light_cycle = []
    dark_cycle.extend(n_moving_stdev(list_CBT('8-19-14', mouse, 'Dark Cycle', master_tt_dic), n_stdev))
    dark_cycle.extend(n_moving_stdev(list_CBT('8-20-14', mouse, 'Dark Cycle', master_tt_dic), n_stdev))
    light_cycle.extend(n_moving_stdev(list_CBT('8-18-14', mouse, 'Light Cycle', master_tt_dic), n_stdev))
    light_cycle.extend(n_moving_stdev(list_CBT('8-19-14', mouse, 'Light Cycle', master_tt_dic), n_stdev))
    print "Dark Cycle:", np.mean(dark_cycle)
    print "Light Cycle:", np.mean(light_cycle)

def get_all_last_2_cycles_moving_stdev(master_tt_dic, acyline_mice, veh_mice, n_stdev):#later--into excel, not print
    """Prints each cycle's average moving standard deviation for all mice, also prints the
    treatment group the mouse belonged to"""
    print "Acyline Mice-Avg of last two full days' %s pt. moving stdev" %(str(n_stdev))
    for mouse in acyline_mice:
        print "Mouse", mouse
        get_last_two_cycles_moving_stdev(master_tt_dic, mouse, n_stdev)
        
    print "Vehicle Mice-Avg of last two full days' %s pt. moving stdev" %str((n_stdev))
    for mouse in veh_mice:
        print "Mouse", mouse
        get_last_two_cycles_moving_stdev(master_tt_dic, mouse, n_stdev)
        
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
                if day=='8-14-14' and cycle=='Light Cycle':
                    #8-14 light is first treatment cycle
                    continue #exits this for loop
                if day=='8-11-14 Light only' and cycle=='Dark Cycle':
                    # 8-11 dark doesn't have any data
                    continue #do not use break, doesn't work with Dark Cycle
                if day =='8-20-14' and cycle=='Light Cycle':
                    #8-20 light cycle is incomplete
                    continue
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

def plot_each_treatment_last_days(last_four_days_pre, last_four_days_post, times, veh_mice,
                                  treatment_mice, sample_frequency, all_times_dic):
    """Given lists of pre and post treatment days, a list of light cycles, lists of vehicle and
    treatment mice, and all_time or master_tt dictionary, saves two plots. One of pre-treatment
    temperature averages every sample_frequency, and one of post-treatment temperature averages every
    sample_frequency. Each plot has two lines-one of veh and treatment averages."""
    pre_veh_lst =  make_organized_time_temp_list(last_four_days_pre, times, veh_mice, all_times_dic)
    pre_treatment_lst =  make_organized_time_temp_list(last_four_days_pre, times,
                                                       treatment_mice, all_times_dic)
    post_veh_lst =  make_organized_time_temp_list(last_four_days_post, times, veh_mice, all_times_dic)
    post_treatment_lst =  make_organized_time_temp_list(last_four_days_post, times,
                                                        treatment_mice,all_times_dic)
    #necessary for controling axis
    fig = plt.figure()
    ax = fig.gca()
        
    x_times = parse_list( make_last_days_x_list(pre_veh_lst), sample_frequency)
    y_veh = parse_list( make_last_days_y_list(pre_veh_lst), sample_frequency)
    y_treatment = parse_list( make_last_days_y_list(pre_treatment_lst), sample_frequency)
    #blue is vehicle, red is treatment
    plt.figure()
    plt.plot(x_times, y_veh, 'b-', label='Vehicle')
    plt.plot(x_times, y_treatment, 'r-', label='Acyline')
    plt.legend()
    plt.ylim(35, 38.5, .5)
    plt.xticks(np.arange(min(x_times), max(x_times)+120, 1440))
    plt.title('Avg. across treatment every 30 seconds-pre')
    plt.xlabel('Time (every 1440 is 12 hrs)')
    plt.ylabel('Average CBT in deg C')
    plt.savefig('last_pre_days_avg.png')
    #plt.show()

    plt.figure()
    x_times = parse_list( make_last_days_x_list(post_veh_lst), sample_frequency)
    y_veh = parse_list( make_last_days_y_list(post_veh_lst), sample_frequency)
    y_treatment = parse_list( make_last_days_y_list(post_treatment_lst), sample_frequency)
    #blue is vehicle, red is treatment
    plt.figure()
    plt.plot(x_times, y_veh, 'b-',label = 'Vehicle')
    plt.plot(x_times, y_treatment, 'r-', label = 'Acyline')
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

def plot_stdev_each_treatment_last_days(last_four_days_pre, last_four_days_post, times, veh_mice,
                                  treatment_mice, sample_frequency, all_times_dic):
    """Given lists of pre and post treatment days, a list of light cycles, lists of vehicle and
    treatment mice, and all_time or master_tt dictionary, saves two plots. One of pre-treatment
    temperature averages every sample_frequency, and one of post-treatment temperature averages every
    sample_frequency. Each plot has two lines-one of veh and treatment averages."""
    pre_veh_lst =  make_stdev_organized_time_temp_list(last_four_days_pre, times, veh_mice, all_times_dic)
    pre_treatment_lst =  make_stdev_organized_time_temp_list(last_four_days_pre, times,
                                                       treatment_mice, all_times_dic)
    post_veh_lst =  make_stdev_organized_time_temp_list(last_four_days_post, times, veh_mice, all_times_dic)
    post_treatment_lst =  make_stdev_organized_time_temp_list(last_four_days_post, times,
                                                        treatment_mice,all_times_dic)
    #necessary for controling axis
    fig = plt.figure()
    ax = fig.gca()
        
    x_times = parse_list( make_last_days_x_list(pre_veh_lst), sample_frequency)
    y_veh = parse_list( make_last_days_y_list(pre_veh_lst), sample_frequency)
    y_treatment = parse_list( make_last_days_y_list(pre_treatment_lst), sample_frequency)
    #blue is vehicle, red is treatment
    plt.figure()
    plt.plot(x_times, y_veh, 'b-', label='Vehicle')
    plt.plot(x_times, y_treatment, 'r-', label='Acyline')
    plt.legend()
    plt.ylim(0, 3, .25)
    plt.xticks(np.arange(min(x_times), max(x_times)+120, 1440))
    plt.title('Stdev across treatment every 30 seconds-pre')
    plt.xlabel('Time (every 1440 is 12 hrs)')
    plt.ylabel('Stdev of CBT in deg C')
    plt.savefig('last_pre_days_stdev.png')
    #plt.show()

    plt.figure()
    x_times = parse_list( make_last_days_x_list(post_veh_lst), sample_frequency)
    y_veh = parse_list( make_last_days_y_list(post_veh_lst), sample_frequency)
    y_treatment = parse_list( make_last_days_y_list(post_treatment_lst), sample_frequency)
    #blue is vehicle, red is treatment
    plt.figure()
    plt.plot(x_times, y_veh, 'b-',label = 'Vehicle')
    plt.plot(x_times, y_treatment, 'r-', label = 'Acyline')
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

def overall_expt_plot(day_labels, times, veh_mice, treatment_mice, sample_frequency, all_times_dic):
    """Given lists of all the days, all the times, all vehicle mice, all acyline mice, an integer of
    how often to plot the data points, and the all_times_dic, plots the average CBT of each
    treatment at each time point for the entire experiment."""
    veh_lst =  make_organized_time_temp_list(day_labels, times, veh_mice, all_times_dic)
    treatment_lst =  make_organized_time_temp_list(day_labels, times,treatment_mice,all_times_dic) 
    plt.figure()
    x_times = parse_list( make_last_days_x_list(veh_lst), sample_frequency)
    y_veh = parse_list( make_last_days_y_list(veh_lst), sample_frequency)
    y_treatment = parse_list( make_last_days_y_list(treatment_lst), sample_frequency)
    #blue is vehicle, red is treatment
    plt.figure()
    plt.plot(x_times, y_veh, 'b-',label = 'Vehicle')
    plt.plot(x_times, y_treatment, 'r-', label = 'Acyline')
    plt.legend()
    plt.xticks(np.arange(min(x_times), max(x_times)+120, 2880), rotation=45)
    plt.ylim(35, 38.5, 0.5)
    plt.title('Avg. across treatment every 30 seconds-entire experiment')
    plt.xlabel('Time (every 2880 is 24 hrs)')
    plt.ylabel('Average CBT in deg C')
    plt.savefig('avg_temp_per_pt_entire_expt.png')

def overall_expt_plot_stdev(day_labels, times, veh_mice, treatment_mice, sample_frequency, all_times_dic):
    """Given lists of all the days, all the times, all vehicle mice, all acyline mice, an integer of
    how often to plot the data points, and the all_times_dic, plots the stdev of the CBT of each
    treatment at each time point for the entire experiment."""
    veh_lst =  make_stdev_organized_time_temp_list(day_labels, times, veh_mice, all_times_dic)
    treatment_lst =  make_stdev_organized_time_temp_list(day_labels, times,treatment_mice,all_times_dic) 
    plt.figure()
    x_times = parse_list( make_last_days_x_list(veh_lst), sample_frequency)
    y_veh = parse_list( make_last_days_y_list(veh_lst), sample_frequency)
    y_treatment = parse_list( make_last_days_y_list(treatment_lst), sample_frequency)
    #blue is vehicle, red is treatment
    plt.figure()
    plt.plot(x_times, y_veh, 'b-',label = 'Vehicle')
    plt.plot(x_times, y_treatment, 'r-', label = 'Acyline')
    plt.legend()
    plt.xticks(np.arange(min(x_times), max(x_times)+120, 2880), rotation=45)
    plt.ylim(0, 2.5, 0.25)
    plt.title('Stdev across treatment every 30 seconds-entire experiment')
    plt.xlabel('Time (every 2880 is 24 hrs)')
    plt.ylabel('Stdev CBT in deg C')
    plt.savefig('stdev_temp_per_pt_entire_expt.png')

##def get_all_mouse_ids():
##    """Looks at all .csv data files and extracts the mouse ids from them"""
##    ###ADD THIS FUNCTION###
##
##def extract_mouse_ids(filename):
##    """Given a string that is a name of a data file, extracts the mouse ids from the first row of
##    the document (by making sure they have the word "Data" in the cell). Returns a list of each
##    unique mouse id as it appears in the file."""
##    mouse_id_set = set()
##    data = csv.DictReader(open(filename, 'rU'), quotechar='"', delimiter = ',')
##    count = 0
##    for row in data:
##        count += 1
##        if count > 1:
##            break
##        print "row:", row
##        for key in row:
##            if 'Data' in key:
##                mouse_id_set.add(key)
##    return sorted(mouse_id_set)    
    
#################
#### MAIN
#################
    
def main():
    """This is the main python code that is run in this program."""
    
    mouse_ids = ["2 Veh Deg. C Data","2 Acyline Deg. C Data","3 Veh Deg. C Data","3 Acyline Deg. C Data",
                 "4 Veh Deg. C Data","4 Acyline Deg. C Data","6 Veh Deg. C Data", "6 Acyline Deg. C Data",
                 "7 Veh Deg. C Data","7 Acyline Deg. C Data","8 Veh Deg. C Data","8 Acyline Deg. C Data",
                 "9 Veh Deg. C Data","9 Acyline Deg. C Data","10 Veh Deg. C Data",
                 "10 Acyline Deg. C Data","11 Veh Deg. C Data","11 Acyline Deg. C Data",
                 "12 Veh Deg. C Data","12 Acyline Deg. C Data","13 Veh Deg. C Data",
                 "13 Acyline Deg. C Data","14 Veh Deg. C Data","14 Acyline Deg. C Data",
                 "16 Veh Deg. C Data","16 Acyline Deg. C Data","17 Veh Deg. C Data",
                 "17 Acyline Deg. C Data", "18 Veh Deg. C Data", "18 Acyline Deg. C Data",
                 "19 Veh Deg. C Data","19 Acyline Deg. C Data","20 Veh Deg. C Data",
                 "20 Acyline Deg. C Data",]
    
##    mouse_ids = #get_all_mouse ids ##MAKE THIS FUNCTION##

    mouse_nums = ['2', '3', '4', '6', '7', '9', '10', '11', '12', '13', '14', '16', '17', '18']

    # each file for the Acyline project (except the first) starts at 18:00:00 on the day of the
    # file label, then goes until 17:59:59 the next day
    filenames = ["8-11-14 Light only, Light Cycle.csv",
                 "8-11-14, Dark Cycle.csv",
                 "8-12-14, Dark Cycle.csv",
                 "8-13-14, Dark Cycle.csv",
                 "8-14-14, Dark Cycle.csv",
                 "8-15-14, Dark Cycle.csv",
                 "8-16-14, Dark Cycle.csv",
                 "8-17-14, Dark Cycle.csv",
                 "8-18-14, Dark Cycle.csv",
                 "8-19-14, Dark Cycle.csv",
                 "8-20-14, Dark Cycle.csv",]

    day_labels = [ '8-11-14 Light only', '8-11-14', '8-12-14', '8-13-14', '8-14-14', '8-15-14',
                   '8-16-14', '8-17-14', '8-18-14', '8-19-14', '8-20-14']

    times = ["Dark Cycle", "Light Cycle"]

    treatments = ["Acyline", "Vehicle"]

    acyline_mice = ['2', '4', '9', '11', '14', '17', '18']
    veh_mice = ['3', '6', '7', '10', '12', '13', '16']
    all_mice = ['2', '3', '4', '6', '7', '9', '10', '11', '12', '13', '14', '16', '17', '18']


    treatment_days = ['8-14-14', '8-15-14', '8-16-14', '8-17-14', '8-18-14', '8-19-14', '8-20-14']
    veh_days = ['8-11-14 Light only','8-11-14','8-12-14','8-13-14']

    last_four_pre_days = ['8-11-14 Light only', '8-11-14', '8-12-14', '8-13-14', '8-14-14']
    #throw out dark cycle of 8-11 light only because it has no data
    #throw out light cycle of 8-14 file because this is first round of treatment

    last_four_post_days = ['8-16-14', '8-17-14', '8-18-14', '8-19-14', '8-20-14']
    #throw out light cycle of 8-20 because it's not a full cycle

    n_ints_in_mavg = extract_ints_in_mavg('user_modify.csv')
    #this defines how many points to be used in calculating moving averages
    
    n_stdev = extract_ints_in_moving_stdev('user_modify.csv')
    #this defines how many points to be used in calculating moving standard deviation

    master_tt_dic = make_master_tt_dic(filenames, mouse_ids)
    find_all_avgs_ers(day_labels, mouse_nums, times, master_tt_dic) #modify to make excel doc, NOT print
    mav_master_dic = make_mav_master_dic(day_labels, mouse_nums, times, master_tt_dic, n_ints_in_mavg)
    all_avg_plots(master_tt_dic)
    
    plot_n_moving_stdv(day_labels, all_mice, times, master_tt_dic, n_stdev)

    get_all_last_2_cycles_moving_stdev(master_tt_dic, acyline_mice, veh_mice, n_stdev)
    
    ##################################################################################################
    all_times_dic = make_all_times_dic(filenames, mouse_ids)

    plot_each_treatment_last_days(last_four_pre_days, last_four_post_days, times, veh_mice,
                                 acyline_mice, 1, all_times_dic)
    overall_expt_plot(day_labels, times, veh_mice, acyline_mice, 1, all_times_dic) 
    plot_stdev_each_treatment_last_days(last_four_pre_days, last_four_post_days, times, veh_mice,
                                        acyline_mice, 1, all_times_dic)
    overall_expt_plot_stdev(day_labels, times, veh_mice, acyline_mice, 1, all_times_dic)
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
