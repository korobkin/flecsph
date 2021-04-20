#!/usr/bin/env python3

# Copyright (c) 2017 Triad National Security, LLC
# All rights reserved.

# Author : Alexander Kaltenborn
# Date : July.8.2020, Updated : April.13.2021
# This is an attempt to have a python script to make a wide variety of binaries
#   with generous control over characteristics

"""
Takes a(/two) single star h5part file(/s) and produces binary system at a specified orbital separation:
"""

#  IMPORTANT if there are asymmetries in initial files, orient this way:
#            M1 ->            COM             <- M2 (-dir -1)
#           ***                                  ***
#        **********                          **********
#       *************          x           *************
#        **********                          **********
#           ***                                  ***
#

# Store data in h5part format
import h5py, sys, argparse, os
import numpy as np

# IMPORTANT CONSTANTS
G_newt  = 6.67430e-8 #cm3/g/s2
M_solar = 1.988435e33 #g
R_solar = 6.957e10 #cm
c_light = 2.99792458e10 #cm/s

USE_PN = False
USE_TIDALLY_LOCKED = True

narange  = np.arange; narray = np.array; nzero = np.zeros
nsum     = np.sum;      nmin = np.min;    nmax = np.max;     nmean = np.mean
nsign    = np.sign;     ncos = np.cos;    nsin = np.sin;  narctan2 = np.arctan2
ln       = np.log;     nsqrt = np.sqrt;   nabs = np.abs

my_description = """
Takes a(/two) single star h5part file(/s) and produces binary system at a specified orbital separation:
  star.h5part (left-justified) [+ star2.h5part (right-justified)]  ==>   binary.h5part"""
my_usage = """
    %(prog)s [-f|--file <filename(s)>] [-i|--identical] [-pm|--pointmass <val>] [-a|--orbsep <val>] [-h|--help]"""

parser = argparse.ArgumentParser(description=my_description, usage=my_usage, epilog="EXAMPLE: $python %(prog)s -f wd.h5part -i -a 1.e9")

parser.add_argument("-f",     "--file",         action="store",     type=argparse.FileType('r'),                        help="input file(s) in FleCSPH hdf5 format",         nargs="+",dest="infile")
parser.add_argument("-of",    "--outfile",      action="store",     type=argparse.FileType('w'),default="binary.h5part",help="output file in FleCSPH hdf5 format",                     dest="outfile")
parser.add_argument("-i",     "--identical",    action="store_true",                            default=False,          help="set this if you want an equal mass binary",              dest="ident")
parser.add_argument("-tl",    "--tidal",        action="store_true",                            default=True,           help="set this if you want the binary to be tidally locked",   dest="tl")
parser.add_argument("-pn",    "--postnewt",     action="store_true",                            default=False,          help="set this if you want post newtonian correction",         dest="pn")
parser.add_argument("-pm",    "--pointmass",    action="store",     type=float,                 default=-1.,            help="total mass of the point star (default: -1.)",            dest="pm")
parser.add_argument("-a",     "--orbsep",       action="store",     type=float,                 default=-1.,            help="orbital separation in cm.",                              dest="orbsep")
parser.add_argument("-pmd",   "--pointmassrho", action="store",     type=float,                 default=4.e+14,         help="density of the point mass (default: 4e14)",              dest="pmd")
parser.add_argument("-pmh",   "--pointmassh",   action="store",     type=float,                 default=1.,             help="smoothing len of the point mass (default: 1)",           dest="pmh")
parser.add_argument("-pmu",   "--pointmassu",   action="store",     type=float,                 default=0.,             help="int. energy of the point mass (default: 0)",             dest="pmu")
parser.add_argument("-pmp",   "--pointmasspres",action="store",     type=float,                 default=2.e+28,         help="pressure of the point mass (default: 2e28)",             dest="pmp")
parser.add_argument("-pmabar","--pointmassabar",action="store",     type=int,                   default=52,             help="Abar of the point mass (default: 52)",                   dest="pmabar")
parser.add_argument("-pmye",  "--pointmassye",  action="store",     type=float,                 default=0.5,            help="Ye of the point mass (default: 0.5)",                    dest="pmye")
parser.add_argument("-dir",   "--direction_2",  action="store",     type=int,                   default=1,              help="orientation of the 2nd star (+/-1)*x_pos (default:+1)",  dest="dir2")
args = parser.parse_args()

def read_H5data(h5_in,string,file,dataset,size,default=0):
    # take the data from the h5part file and store it in an array of the correct size
    try:
        array = h5_in[dataset+string];
    except:
        print("No data for {} in {} setting to default = {}".format(string.split("/")[1],file,default))
        array = nzero((size))
        array[:] = default
    return array

def store_data(string,file,h5data,size,default=0):
    # take the data from the h5part file and store it in an array of the correct size
    try:
        if(size == len(h5data[()])):
            array = narray(h5data[()])
    except:
        print("No data for {} in {} setting to default = {}".format(string,file,default))
        array = nzero((size))
        array[:] = default
    return array

def get_reduced_mass(m1,m2):
    # Reduced mass mu = m1*m2/(m1+m2)
    mt = m1 + m2
    return (m1*m2)/mt

def get_center_of_mass(mtot,m,x,y,z):
  # Find the center of mass of a star
  x_com = nsum(x*m)/mtot
  y_com = nsum(y*m)/mtot
  z_com = nsum(z*m)/mtot
  return x_com,y_com,z_com

def center_star(mtot,m,x,y,z):
  # Re-centers star so it sits at the origin
  x_com,y_com,z_com = get_center_of_mass(mtot,m,x,y,z)
  xnew = x - x_com
  ynew = y - y_com
  znew = z - z_com
  return xnew,ynew,znew

def get_offsets(m1,m2,separation):
  # Returns radii of mtot1 and mtot2 in binary
  # r1,r2 are DISPLACEMENTS, not magnitudes
  mt = m1 + m2
  r2 = separation*m1/mt
  r1 = r2 - separation
  return r1,r2

def separate(r1,r2,x1,y1,z1,x2,y2,z2):
  #Separate stars by radius and return new position arrays.
  #By default we separate along x-axis.
  x1[:] = x1[:] + r1
  y1[:] = y1[:]
  z1[:] = z1[:]
  x2[:] = x2[:] + r2
  y2[:] = y2[:]
  z2[:] = z2[:]

def get_keplerian_velocity(m_mine,m_other,r_mine,separation):
  # The tangential velocity for stable circular orbit
  mt = m_mine + m_other
  angular_velocity = nsqrt(G_newt*mt/(separation**3.))
  return angular_velocity*r_mine

def get_pn_velocity_angular(m_mine,m_other,r_mine,separation):
  # The tangential velocity for initial data for PN2.5
  mt     = m_mine + m_other
  eta    = (m_mine/mt)*(m_other/mt)
  factor = G_newt*mt/(separation*c_light**2)
  omega_kepp = nsqrt(G_newt*mt/separation**3)
  omega = omega_kepp*(1 - 0.5*factor*(3. - eta)
                      + (factor**2)*(15. + 47.*eta + 3.*(eta**2))/8.)
  vtan = omega*r_mine
  return vtan

def get_pn_rdot(m_mine,m_other,separation):
  mt   = m_mine + m_other
  eta  = (m_mine/mt)*(m_other/mt)
  rdot =-(64./5.)*(G_newt**3)*(mt**2)*m_other*eta/((separation**3)*(c_light**5))
  return rdot

def get_velocities(m_mine,m_other,separation,x,y,z,rstar):
  # Get the velocities of a particle.
  # Tangential velocity is along y-axis.
  # For counterclockwise rotation, vy > 0 iff x > 0
  if USE_TIDALLY_LOCKED:
      r = nsqrt(x**2 + y**2)
      theta = narctan2(y,x)
  else:
      r = nabs(rstar)
      theta = 0 if nsign(rstar) > 0 else np.pi
  if USE_PN:
      vtan = get_pn_velocity_angular(m_mine,m_other,r,separation)
      vrad = get_pn_rdot(m_mine,m_other,separation)
  else:
      vtan = get_keplerian_velocity(m_mine,m_other,r,separation)
      vrad = 0
  vx = vrad*ncos(theta) - vtan*nsin(theta)
  vy = vrad*nsin(theta) + vtan*ncos(theta)
  vz = 0
  return vx,vy,vz

def main():
  # read the input file
  try:
    f1 = h5py.File(args.infile[0].name,'r')
  except:
    sys.exit ("ERROR: cannot read first input file %s" % args.infile)
  if(len(args.infile) == 2):
    try:
      f2 = h5py.File(args.infile[1].name,'r')
    except:
      sys.exit ("ERROR: cannot read second input file %s" % args.infile)

  try:
    os.remove("binary.h5part")
  except:
    pass

  if (args.orbsep == -1.0):
    sys.exit ("ERROR: no orbital separation supplied in arguments")

  USE_PN = args.pn
  print("Using post newtonian correction = {}".format(USE_PN))
  USE_TIDALLY_LOCKED = args.tl
  print("Assuming the binary is tidally locked = {}".format(USE_TIDALLY_LOCKED))

  # Open the first input file HDF5
  dataset1 = ("Step#%d" % (len(list(f1.keys()))-1))
  x1  = f1[dataset1+"/x"];   x1 = narray(x1[()])
  y1  = f1[dataset1+"/y"];   y1 = narray(y1[()])
  size1 = len(x1)
  z1  = read_H5data(f1,"/z",   "file1",dataset1,size1);      z1 = store_data("z",  "file1", z1,size1)
  vx1 = read_H5data(f1,"/vx",  "file1",dataset1,size1);     vx1 = store_data("vx", "file1",vx1,size1)
  vy1 = read_H5data(f1,"/vy",  "file1",dataset1,size1);     vy1 = store_data("vy", "file1",vy1,size1)
  vz1 = read_H5data(f1,"/vz",  "file1",dataset1,size1);     vz1 = store_data("vz", "file1",vz1,size1)
  P1  = read_H5data(f1,"/P",   "file1",dataset1,size1);      P1 = store_data("P",  "file1", P1,size1)
  d1  = read_H5data(f1,"/rho", "file1",dataset1,size1);      d1 = store_data("rho","file1", d1,size1)
  m1  = read_H5data(f1,"/m",   "file1",dataset1,size1);      m1 = store_data("m",  "file1", m1,size1)
  h1  = read_H5data(f1,"/h",   "file1",dataset1,size1);      h1 = store_data("h",  "file1", h1,size1)
  u1  = read_H5data(f1,"/u",   "file1",dataset1,size1);      u1 = store_data("u",  "file1", u1,size1)
  a1  = read_H5data(f1,"/Abar","file1",dataset1,size1,12);   a1 = store_data("Abar","file1",a1,size1,12)
  ye1 = read_H5data(f1,"/Ye",  "file1",dataset1,size1,0.5); ye1 = store_data("Ye", "file1",ye1,size1,0.5)
  Temp1 = read_H5data(f1,"/temp","file1",dataset1,size1,1e4); Temp1 = store_data("temp","file1",Temp1,size1,1e4)
  type1  = nzero((size1));  type1[:] = 0
  state1 = nzero((size1)); state1[:] = 1
  # Try to open the second input file HDF5
  if(len(args.infile) == 2):
    dataset2 = ("Step#%d" % (len(list(f2.keys()))-1))
    x2  = f2[dataset2+"/x"];   x2 = narray(x2[()])
    y2  = f2[dataset2+"/y"];   y2 = narray(y2[()])
    size2 = len(x2)
    z2  = read_H5data(f2,"/z",   "file2",dataset2,size2);      z2 = store_data("z",  "file2", z2,size2)
    vx2 = read_H5data(f2,"/vx",  "file2",dataset2,size2);     vx2 = store_data("vx", "file2",vx2,size2)
    vy2 = read_H5data(f2,"/vy",  "file2",dataset2,size2);     vy2 = store_data("vy", "file2",vy2,size2)
    vz2 = read_H5data(f2,"/vz",  "file2",dataset2,size2);     vz2 = store_data("vz", "file2",vz2,size2)
    P2  = read_H5data(f2,"/P",   "file2",dataset2,size2);      P2 = store_data("P",  "file2", P2,size2)
    d2  = read_H5data(f2,"/rho", "file2",dataset2,size2);      d2 = store_data("rho","file2", d2,size2)
    m2  = read_H5data(f2,"/m",   "file2",dataset2,size2);      m2 = store_data("m",  "file2", m2,size2)
    h2  = read_H5data(f2,"/h",   "file2",dataset2,size2);      h2 = store_data("h",  "file2", h2,size2)
    u2  = read_H5data(f2,"/u",   "file2",dataset2,size2);      u2 = store_data("u",  "file2", u2,size2)
    a2  = read_H5data(f2,"/Abar","file2",dataset2,size2,12);   a2 = store_data("Abar","file2",a2,size2,12)
    ye2 = read_H5data(f2,"/Ye",  "file2",dataset2,size2,0.5); ye2 = store_data("Ye", "file2",ye2,size2,0.5)
    Temp2 = read_H5data(f2,"/temp","file2",dataset2,size2,1e4); Temp2 = store_data("temp","file2",Temp2,size2,1e4)
    type2  = nzero((size2));  type2[:] = 0
    state2 = nzero((size2)); state2[:] = 2
  else:
    if(args.ident):
      x2     = nzero((size1));     x2[:] =     x1[:]
      y2     = nzero((size1));     y2[:] =     y1[:]
      z2     = nzero((size1));     z2[:] =     z1[:]
      vx2    = nzero((size1));    vx2[:] =    vx1[:]
      vy2    = nzero((size1));    vy2[:] =    vy1[:]
      vz2    = nzero((size1));    vz2[:] =    vz1[:]
      P2     = nzero((size1));     P2[:] =     P1[:]
      d2     = nzero((size1));     d2[:] =     d1[:]
      m2     = nzero((size1));     m2[:] =     m1[:]
      h2     = nzero((size1));     h2[:] =     h1[:]
      u2     = nzero((size1));     u2[:] =     u1[:]
      a2     = nzero((size1));     a2[:] =     a1[:]
      ye2    = nzero((size1));    ye2[:] =    ye1[:]
      Temp2  = nzero((size1));  Temp2[:] =  Temp1[:]
      type2  = nzero((size1));  type2[:] =  type1[:]
      state2 = nzero((size1)); state2[:] =         2
      size2  = size1
      args.dir2 = -1
    else:
      x2     = nzero((1));     x2[0] = 0.
      y2     = nzero((1));     y2[0] = 0.
      z2     = nzero((1));     z2[0] = 0.
      vx2    = nzero((1));    vx2[0] = 0.
      vy2    = nzero((1));    vy2[0] = 0.
      vz2    = nzero((1));    vz2[0] = 0.
      P2     = nzero((1));     P2[0] = args.pmp
      d2     = nzero((1));     d2[0] = args.pmd
      m2     = nzero((1));     m2[0] = args.pm
      h2     = nzero((1));     h2[0] = args.pmh
      u2     = nzero((1));     u2[0] = args.pmu
      a2     = nzero((1));     a2[0] = args.pmabar
      ye2    = nzero((1));    ye2[0] = args.pmye
      Temp2  = nzero((1));  Temp2[0] = 0.
      type2  = nzero((1));  type2[0] = 0
      state2 = nzero((1)); state2[0] = 3
      size2  = 1
  print("Calculating total mass star 1...")
  M_star1 = nsum(m1)
  print("Mass of star 1: {0:3.2f} Ms".format(M_star1/M_solar))
  print("Calculating total mass star 2...")
  M_star2 = nsum(m2)
  print("mass of star 2: {0:3.2f} Ms".format(M_star2/M_solar))
  newsize = size1 + size2
  part_id = narange(newsize)
  Mtot = M_star1 + M_star2
  sep = args.orbsep
  print("Calculating X and Y coordinates.")
  x1,y1,z1 = center_star(Mtot,m1,x1,y1,z1)
  x2,y2,z2 = center_star(Mtot,m2,x2,y2,z2)
  x2 *= args.dir2
  y2 *= args.dir2
  r1,r2 = get_offsets(M_star1,M_star2,sep)
  separate(r1,r2,x1,y1,z1,x2,y2,z2)
  # get velocities
  print("Getting Keplerian velocities.")
  vx1,vy1,vz1 = get_velocities(M_star1,M_star2,sep,x1,y1,z1,r1)
  vx2,vy2,vz2 = get_velocities(M_star2,M_star1,sep,x2,y2,z2,r2)

  temparray = nzero((newsize))
  # Open the output file HDF5
  print("Writing out to the file '{}'".format(args.outfile.name))
  out = h5py.File(args.outfile.name,'w')
  grp = out.create_group("Step#0")
  print("  setting position x, y, z...")
  temparray[:size1] = x1; temparray[size1:] = x2
  grp.create_dataset("x",data=temparray)
  temparray[:size1] = y1; temparray[size1:] = y2
  grp.create_dataset("y",data=temparray)
  temparray[:size1] = z1; temparray[size1:] = z2
  grp.create_dataset("z",data=temparray)
  print("  setting velocity vx, vy, vz...")
  temparray[:size1] = vx1; temparray[size1:] = vx2
  grp.create_dataset("vx",data=temparray)
  temparray[:size1] = vy1; temparray[size1:] = vy2
  grp.create_dataset("vy",data=temparray)
  temparray[:size1] = vz1; temparray[size1:] = vz2
  grp.create_dataset("vz",data=temparray)
  print("  setting pressure...")
  temparray[:size1] = P1; temparray[size1:] = P2
  grp.create_dataset("P",data=temparray)
  print("  setting density...")
  temparray[:size1] = d1; temparray[size1:] = d2
  grp.create_dataset("rho",data=temparray)
  print("  setting mass...")
  temparray[:size1] = m1; temparray[size1:] = m2
  grp.create_dataset("m",data=temparray)
  print("  setting smoothing length...")
  temparray[:size1] = h1; temparray[size1:] = h2
  grp.create_dataset("h",data=temparray)
  print("  setting int. energy...")
  temparray[:size1] = u1; temparray[size1:] = u2
  grp.create_dataset("u",data=temparray)
  print("  setting Abar...")
  temparray[:size1] = a1; temparray[size1:] = a2
  grp.create_dataset("Abar",data=temparray)
  print("  setting Ye...")
  temparray[:size1] = ye1; temparray[size1:] = ye2
  grp.create_dataset("ye",data=temparray)
  print("  setting temperature...")
  temparray[:size1] = Temp1; temparray[size1:] = Temp2
  grp.create_dataset("temp",data=temparray)
  print("  setting type...")
  temparray[:size1] = type1; temparray[size1:] = type2
  grp.create_dataset("type",data=temparray)
  print("  setting state...")
  temparray[:size1] = state1; temparray[size1:] = state2
  grp.create_dataset("state",data=temparray)
  print("  setting id")
  grp.create_dataset("id",data=part_id)

  print("Done creating hdf5 file.")
  out.close()
  f1.close()
  if(len(args.infile) == 2):
    f2.close()

if __name__ == "__main__":
    main()
