import simulation_class as smc
import argparse
import os

parser = argparse.ArgumentParser(description='Script for modifying, compiling lt and xyz template files and copying the LAMMPS input script output into a standalone simulation folder')
parser.add_argument('-N','--total_monomers',dest='total_monomers',type=int,default=3000, help="Total number of monomers in the system")
parser.add_argument('-fA','--fractionA',dest='fA',type=float,default=0.5, help="Fraction of monomers of species A")
parser.add_argument('--epsAA',dest='eAA',type=float,default=1.0,help="Attraction between monomers AA in kbT")
parser.add_argument('--epsBB',dest='eBB',type=float,default=1.0,help="Attraction between monomers BB in kbT")
parser.add_argument('--epsAB',dest='eAB',type=float,default=1.0,help="Attraction between monomers AB in kbT")
parser.add_argument('--angle',dest='angle',type=int,default=20,help="Angle strength of polymer chain in kbT")
parser.add_argument('-p',dest='p',type=float,default=0.9,help="Monomer conversion where simulation will end")
parser.add_argument('-f','--folder',dest='dest_folder',help='Destination folder for simulation')
parser.add_argument('--servername',dest='servername',default=None,help='If sending simulation to a remote server this is the IP/hostname of the server.')
parser.add_argument('--send_to_cluster',dest='send_to_cluster',action='store_true',help='If True send the project to the destination folder on the remote server specified by dest_folder')
parser.add_argument('--xyz',dest='xyz_folder',help='Path to folder containing xyz file for initial configuration')
parser.add_argument('--lt',dest='lt_folder',help='Path to folder containing lt files for simulation')
parser.add_argument('--slurm',dest='slurm',action='store_true',help='If true use slurm queueing system else run job from console')

args = parser.parse_args()

epsilons = (args.eAA,args.eBB,args.eAB)

project_path = os.path.abspath(os.path.expanduser(args.dest_folder))

sim = smc.Simulation(total_monomers=args.total_monomers,
                     monomer_attractions=epsilons,
                     p = args.p,monomer_A_fraction=args.fA,
                     angle_strength=args.angle,
                     lt_dir = os.path.expanduser(args.lt_folder),
                     xyz_dir = os.path.expanduser(args.xyz_folder),
                     send_to_cluster=args.send_to_cluster,
                     servername=args.servername)

sim.compile_simulation()
if not args.send_to_cluster:
    sim.move_simulation_files(project_path,slurm=args.slurm)
    sim.start_simulation(slurm  = args.slurm,singularity=os.path.expanduser("~/copoly_bondreact_kspace.sif"))
else:
    sim.move_simulation_files_remote(project_path,slurm=args.slurm)
    sim.start_simulation_remote()
