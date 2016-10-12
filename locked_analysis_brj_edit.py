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

GRID_TRANSITION = "SEEKR: GRID TRANSITION: "
GRID_TRANSITION_COMPILE = re.compile(GRID_TRANSITION)
FORWARD_OUTPUT_GLOB = "fwd_rev*.out.*" # used to find all the forward phase simulation output


def parse_md_transitions():
  'find all forward phase simulation output files'
  #forward_dir_glob = os.path.join(self.directory,'md','fwd_rev',FORWARD_OUTPUT_GLOB)
  #forward_output_filenames = glob.glob(forward_dir_glob)
  forward_output_filenames= glob.glob(FORWARD_OUTPUT_GLOB)
 # read files and sift out the transition lines
  transition_lines = [] # a list containing all lines of the output that mark transitions
  transitions= []
  for filename in forward_output_filenames:
    print 'output filename', filename
    for line in open(filename,'r'):
      if re.match(GRID_TRANSITION_COMPILE, line):
        transition_lines.append(line)
  # feed the lines into the Transition object
    for line in transition_lines:
      transitions.append(Transition(line))
  return transitions
#  line=line.strip()
#  linetail = line[len(GRID_TRANSITION)-1:] # just take the last bit of the line, the important part, but not the endline
#  linelist = linetail.split(',') # split line into a list of elements
#  lineID = linelist.split(', ID:':-1)
#  return lineID
#forward_filename = forward..0.dcd
#forward_dcd = forward_filename[:8] + lineID + foward_filename[8:]
#forward_dir_dcd = os.path.join(self.directory,'md','fwd_rev',forward_dcd)


class Transition():
  ''' object representing a transition event between one milestone to another. Used to construct transition statistics'''
  def __init__(self, line):
    line = line.strip() # we have to parse the line
    self.line = line
    linetail = line[len(GRID_TRANSITION)-1:] # just take the last bit of the line, the important part, but not the endline
    linelist = linetail.split(',') # split line into a list of elements
    dictlist = map(lambda a: a.strip().split(': '), linelist) # map the line to a list of lists for dictionary conversion
    linedict = dict(dictlist) # convert the list of lists into a dictionary
    #pprint(linedict)
    self.src = int(linedict['source'].strip())
    self.dest = int(linedict['destination'].strip())
    self.cur_step = float(linedict['stepnum'].strip())
    self.time = float(linedict['incubation time'].strip().split()[0])
    self.ligand_com = linedict['ligand COM'].strip().split()
    self.receptor_com = linedict['receptor COM'].strip().split()
    self.receptor_start_com = linedict['receptor start COM'].strip().split()
    self.ID= str(linedict['ID'].strip())

  def print_status(self):
    print "src:", self.src
    print "dest:", self.dest
    print "cur_step:", self.cur_step
    print "time:", self.time
    print "ligand_com:", self.ligand_com
    print "receptor_com:", self.receptor_com
    print "receptor_start_com:", self.receptor_start_com

def get_md_transition_statistics(transitions):
  'parse the transition data to obtain transition counts'
  counts = {} # all the sources and their destinations
  total_counts = {} # keeps track of all counts to any destination
  total_times = {} # the total time of all counts, to be averaged later
  avg_times = {} # the times to transition out of each source
  #site_index = self.site
  for transition in transitions:
    source = transition.src
    dest = transition.dest
    time = transition.time
    src_key = '%d_%d' % (site_index, source)
    dest_key = '%d_%d' % (site_index, dest)
    if src_key in counts.keys():
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

def extract_info(struct_filename,traj_filename):
  ref = mda.Universe(struct_filename)
  mobile = mda.Universe(struct_filename,traj_filename,topology_format='PRMTOP')
  frame = mda.coordinates.DCD.DCDReader(dcd[])
    mda.analysis.align.alignto(frame, ref, select="resname BCD and name CA", mass_weighted=True)  
    first = frame(0)
    last = frame(-1)
  plane_normal = []
  lig_com = []
  rec_com = []
  for frame in trajectory:
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
  return (plane_normal, lig_com, rec_com)
  return (first, last)

def measure_angle(extract_info):
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

def compare_angles(measure_angle):
  same_secondary=[]
    for angle in compare_angles(theta):
    if (first > 90.0 && last > 90.0):
      same_secondary.append(angle)
  same_primary=[]
    for angle in compare_angles(theta):
    if (first < 90.0 && last < 90.0):
       same_primary.append(angle)
  primary_secondary_trans=[]
    for angle in compare_angles(theta):
    if (first < 90.0 && last > 90.0):
       primary_secondary_trans.append(angle)
  secondary_primary_trans=[]
    for angle in compare_angles(theta):
    if (first > 90.0 && last < 90.0):
       secondary_primary.append(angle)
  return (same_secondary, same_primary, primary_secondary_trans, secondary_primary_trans)


def main():
  parser= argparse.ArgumentParser(description='Analyzes forward reverse output for systems with locked milestones')
  parser.add_argument('anchor', metavar='ANCHOR', type=str, help="the anchor that will be analyzed")
  parser.add_argument('-m', '--milestones', dest="milestones", type=str, help="Milestones file") # This should contain most of what the user needs
  args = parser.parse_args() # parse all the arguments
  args = vars(args)

  anchor= args['anchor']
  print "Performing analysis for anchor: " , anchor
  transitions= parse_md_transitions()
  Transition.print_status(transitions[0])
  #get_md_transitions(transitions)
  struct_filename= '../building/holo.prmtop'
  for transition in transitions:
    #print "ID:", str(transition.ID)
    traj_filename= str("forward."+transition.ID+".0.dcd")
    #print traj_filename 
  



if __name__ == "__main__": main()

