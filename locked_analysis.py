import os, sys, re, random, unittest, glob
import argparse
from string import Template #** Unsure if being used?
from math import exp, log #** currently not used
from pprint import pprint
import numpy as np
from scipy import linalg
from copy import deepcopy #**currently not used-- but normally used for copying dictionaries
import xml.etree.ElementTree as ET
from numpy import linalg as LA
from pprint import pprint
import MDAnalysis as mda

GRID_TRANSITION = "SEEKR: GRID TRANSITION: " #this is the string we are looking for in the output files
GRID_TRANSITION_COMPILE = re.compile(GRID_TRANSITION) #here we are compiling this a regular expression so we can search for it
FORWARD_OUTPUT_GLOB = "fwd_rev*.out.*" # used to find all the forward phase simulation output


def parse_md_transitions(): #This function readsthe output files find trajectories with successful transitions
  'find all forward phase simulation output files'
  forward_output_filenames= glob.glob(FORWARD_OUTPUT_GLOB) #gathers all of the files that match the FORWARD_OUTPUT_GLOB defined above
 # read files and sift out the transition lines
  transition_lines = [] # create a list containing all lines of the output that mark transitions
  transitions= [] #create a list to store all of the transitions
  for filename in forward_output_filenames: #go through each file in the glob and look for transitions
    print 'output filename', filename #a test so we can make sure it is reading the right files
    for line in open(filename,'r'): #open the output file for reading and go through each line
      if re.match(GRID_TRANSITION_COMPILE, line): #if the line matches the GRID_TRANSITION expression we defined above...
        transition_lines.append(line) #add (the text) of that line to the transition_lines list
  # feed the lines into the Transition object-- we are creating an object here because this may come in handy for future uses of the code
    for line in transition_lines: #go through each entry in the transition_lines list
      transitions.append(Transition(line)) #Add each of the successful transitions to the transitions list-- take a look at the Transition class to see what's going on here
  return transitions #pass the transitions list back to the main function where it was called

####This is extra stuff we probably don't need anymore (but let's hang on to it just in case)####
#  line=line.strip()
#  linetail = line[len(GRID_TRANSITION)-1:] # just take the last bit of the line, the important part, but not the endline
#  linelist = linetail.split(',') # split line into a list of elements
#  lineID = linelist.split(', ID:':-1)
#  return lineID
#forward_filename = forward..0.dcd
#forward_dcd = forward_filename[:8] + lineID + foward_filename[8:]
#forward_dir_dcd = os.path.join(self.directory,'md','fwd_rev',forward_dcd)
########################

class Transition():
  ''' object representing a transition event between one milestone to another. Used to construct transition statistics'''
  def __init__(self, line): #we are going to store all of the information below into an object so that we can access it later- The line is what is passed from parse_md_transitions
    line = line.strip() # we have to parse the line
    self.line = line #keep track of which line number we are on
    linetail = line[len(GRID_TRANSITION)-1:] # just take the last bit of the line, the important part, but not the endline
    linelist = linetail.split(',') # split line into a list of elements
    dictlist = map(lambda a: a.strip().split(': '), linelist) # map the line to a list of lists for dictionary conversion
    linedict = dict(dictlist) # convert the list of lists into a dictionary
    #pprint(linedict)
### now we will store each of the valeus we took from the line as entries in a dictionary###
    self.src = int(linedict['source'].strip()) #the milestone we started from
    self.dest = int(linedict['destination'].strip()) #milestone where we ended
    self.cur_step = float(linedict['stepnum'].strip()) #how many steps the simulation ran before stopping
    self.time = float(linedict['incubation time'].strip().split()[0]) #simulation time to transition between milestones
    self.ligand_com = linedict['ligand COM'].strip().split() #center of mass of the ligand
    self.receptor_com = linedict['receptor COM'].strip().split() #center of mass of the receptor
    self.receptor_start_com = linedict['receptor start COM'].strip().split() #center of mass of the receptor when we started the simulation
    self.ID= str(linedict['ID'].strip()) #which forward trajectory (number) 

  def print_status(self): #this will print these values for whatever transition line it is given -- a good check
    print "src:", self.src
    print "dest:", self.dest
    print "cur_step:", self.cur_step
    print "time:", self.time
    print "ligand_com:", self.ligand_com
    print "receptor_com:", self.receptor_com
    print "receptor_start_com:", self.receptor_start_com


###we are not using this right now, but may need it soon###
'''
def get_md_transition_statistics(transitions): #we must give this function our list of Transition objects that are stored in the list transitions
  'parse the transition data to obtain transition counts'
  counts = {} # all the sources and their destinations
  total_counts = {} # keeps track of all counts to any destination
  total_times = {} # the total time of all counts, to be averaged later
  avg_times = {} # the times to transition out of each source
  #site_index = self.site
  for transition in transitions: #for each transition, we will extract the following information
    source = transition.src
    dest = transition.dest
    time = transition.time
    src_key = '%d_%d' % (site_index, source)
    dest_key = '%d_%d' % (site_index, dest)
    if src_key in counts.keys(): ##here we are counting all of the transitions and to which milestones they go
      if dest_key in counts[src_key].keys():
        counts[src_key][dest_key] += 1
      else:
        counts[src_key][dest_key] = 1
      total_times[src_key] += time
      total_counts[src_key] += 1
    else:
      counts[src_key] = {dest_key:1}
      total_counts[src_key] = 1
      total_times[src_key] = time
  for src_key in total_times.keys():
    avg_times[src_key] = total_times[src_key] / total_counts[src_key]

    return counts, total_counts, total_times, avg_times
'''

def extract_info(struct_filename,traj_filename): #here is where we will get the coordinates of the ligand and receptor
  ref = mda.Universe(struct_filename) #creating our MDAnalysis universe for the RMSD reference, giving it the prmtop file
  mobile = mda.Universe(struct_filename,traj_filename,topology_format='PRMTOP') #universe for our trajectory of interest 
  first_frame = mda.coordinates.DCD.DCDReader(dcd[0]) #coordinates for the first frame of the trajectory
  last_frame = mda.coordinates.DCD.DCDReader(dcd[-1]) #coordinates for the last frame of the trajectory 
  mda.analysis.align.alignto(first_frame, ref, select="resname BCD and name CA", mass_weighted=True) ##be careful here what is the variable "frame", you have frame1 and frame2  
  mda.analysis.align.alignto(last_frame, ref, select="resname BCD and name CA", mass_weighted=True)

  return(first_frame,last_frame)

def get_coords():
#this is from the other script-- we are just extracting the coordinates we need and calculating the plane and its normal vector-- using the aligned trajectory from above
  plane_normal = []
  lig_com = []
  rec_com = []
  O3_A= u.select_atoms("resname BCD and name O3").center_of_mass()
  O13_B= u.select_atoms("resname BCD and name O13").center_of_mass()
  O23_C= u.select_atoms("resname BCD and name O23").center_of_mass()
  BA = O3_A -O13_B
  BC = O23_C - O13_B
  plane_normal.append(np.cross(BA , BC))
  lig_com.append(u.select_atoms("resname APN").center_of_mass())
  rec_com.append(u.select_atoms("resname BCD").center_of_mass())

  plane_normal = np.array(plane_normal)
  lig_com = np.array(lig_com)
  rec_com = np.array(rec_com)
  #print 'plane_normal:'
  #pprint(plane_normal)
  #print 'Ligand COM Coordinates:'
  #pprint(lig_com)
  #print 'Receptor COM Coordinates:'
  #pprint(rec_com)

  return (plane_normal_first, lig_com_first, rec_com_first, plane_normal_last, lig_com_last, plane_normal_last)

def measure_angle(): ##straight from the other code-- a function to calculate the orientation angle
  theta = []
  lig_vec = np.subtract(lig_com, rec_com)
    #print 'lig_vec:'
    #pprint(lig_vec)
  for i in range(len(plane_normal)):
    theta.append(np.arccos(np.dot(lig_vec[i],plane_normal[i])/(LA.norm(lig_vec[i])*LA.norm(plane_normal[i]))))
  theta = np.array(theta)
  theta = np.degrees(theta)
  print 'Theta:', theta.shape
  pprint(theta)
  return theta

## should there be two measure_angle functions so that the angles can be compared as they are below? or is there a way to pass the functions and still compare?

## we will probably want to use this to determine which face the ligand is on!

#def compare_angles(measure_angle):
  #same_secondary=[]
   # for angle in compare_angles(theta):
   # if (first > 90.0 && last > 90.0):
    #  same_secondary.append(angle)
 # same_primary=[]
   # for angle in compare_angles(theta):
   # if (first < 90.0 && last < 90.0):
    #   same_primary.append(angle)
 # primary_secondary_trans=[]
   # for angle in compare_angles(theta):
   # if (first < 90.0 && last > 90.0):
     #  primary_secondary_trans.append(angle)
  #secondary_primary_trans=[]
   # for angle in compare_angles(theta):
   # if (first > 90.0 && last < 90.0):
    #   secondary_primary.append(angle)
 # return (same_secondary, same_primary, primary_secondary_trans, secondary_primary_trans)


def main():
#parser allows uss to provide arguments when we run the code from the command line
  parser= argparse.ArgumentParser(description='Analyzes forward reverse output for systems with locked milestones')
  parser.add_argument('anchor', metavar='ANCHOR', type=str, help="the anchor that will be analyzed") #the first argument is the anchor number we want to analyze
  parser.add_argument('-m', '--milestones', dest="milestones", type=str, help="Milestones file") # This should contain most of what the user needs regarding the milestones
  args = parser.parse_args() # parse all the arguments
  args = vars(args) 

  anchor= args['anchor'] #store the anchor number as a variable
  print "Performing analysis for anchor: " , anchor
  transitions= parse_md_transitions()  #calls the parse_md_transitions function and stores whatever is returned as the variable "transitions"
  Transition.print_status(transitions[0]) #just printing to check that we read the file properly
  #get_md_transitions(transitions)
  struct_filename= '../building/holo.prmtop' #the name of the structure we wil give to MDAnalysis
  for transition in transitions: #for each transition that was returned from parse_md_transitions...
    #print "ID:", str(transition.ID)
    traj_filename= str("forward."+transition.ID+".0.dcd") #select the trajectory file that corresponds to that transition-- transition.ID comes from the file we parsed
    print traj_filename #print to check
##Everything appears to be working correctly up to this point in main### 

  plane_normal_first, lig_com_first, rec_com_first= get_coords(first_frame)
  plane_normal_last, lig_com_last, rec_com_last= get_coords(last_frame)

  theta_first = measure_angle(plane_normal_first, rec_com_first, lig_com_first)
  theta_last = measure_angle(plane_normal_last, rec_com_last, lig_com_last)
#### The next thing we need is to use the trajectory above to get our coordinates... which functions (from above) will we need to call to do this and what variables to they requirei??#####



if __name__ == "__main__": main()
