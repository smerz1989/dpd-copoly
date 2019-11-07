import os
import subprocess as sb
import glob
import shutil
import server_class as svc
import re

class Simulation(object):
    """Object encapsulating simulation used for creating and deploying LAMMPS simulation files.

    Parameters
    ----------
    total_monomers : int
        Total number of monomers in the simulation.
    
    monomer_A_fraction : float
        Fraction of monomers in the simulation that have identity A, the remaining monomers will have identity B

    p : float
        The reaction extent at which the simulation will stop

    total_polymers : int
        The total number of polymers in the simulation
    
    polymer_A_fraction : float
        Fraction of polymers with identity A
    """
    def __init__(self,sequences,total_polymers=2,polymer_A_fraction=0.5,monomer_attractions=(1,1,1),angle_strength=20,
                    dump_frequency=1000,xyz_dir = os.path.abspath('../../xyzs/'),lt_dir = os.path.abspath('../../lt_files/'),
                    send_to_cluster=False,servername=None):
        self.sequences = sequences
        self.total_polymers = total_polymers
        self.dump_frequency = dump_frequency
        self.xyz_dir = os.path.abspath(xyz_dir)
        self.lt_dir = os.path.abspath(lt_dir)
        self.eAA,self.eBB,self.eAB = monomer_attractions
        self.angle_strength = angle_strength
        self.send_to_cluster=send_to_cluster
        if self.send_to_cluster:
            self.server_connection = svc.ServerConnection()

    def create_polymer_xyz_files(self,number):
        print('placeholder')

    def create_polymer_lt_files(self,sequences,numbers):
        module_directory, filename = os.path.split(__file__)
        shutil.copy(os.path.join(module_directory,'genpoly_lt.py'),self.lt_dir)
        os.chdir(self.lt_dir)
	sequence_str=""
	for sequence in sequences:
            sequence_str += "\n".join(['ABEAD' if monomer=='A' else 'BBEAD' for monomer in sequence])+'\n'
        with open('sequence.txt','w') as seq_file:
            seq_file.write(sequence_str*number)
        with open('cuts.txt','w') as cut_file:
	    cut_str=""
	    for number, sequence in zip(numbers,sequences):
                seq_length = len(sequence)
                cut_str = '\n'.join([str(i) for i in range(seq_length,seq_length*number,seq_length)])
                cut_file.write(cut_str)
        sb.call(["python","genpoly_lt.py","-header","import copolyff.lt\nimport a_bead.lt\nimport b_bead.lt",
                                          "-inherits","COPOLYFF",
                                          "-bond","COPOLYFF/Backbone","monomer","monomer",
                                          "-sequence","sequence.txt",
                                          "-cuts","cuts.txt",
                                          "-polymer-name","DPDPoly"],stdin=open('./new_coords.raw','r'),stdout=open('dpdpoly.lt','w'))
        sb.call(["sed","-i","s/\\\n/\\n/g","dpdpoly.lt"])


    def compile_simulation(self,packmol_path='packmol'):
        #for sequence in self.sequences:
        self.create_polymer_lt_files(self.sequences,[self.total_polymers*polymer_A_fraction,self.total_polymers-self.total_polymers*polymer_A_fraction])
        self.change_monomer_attraction()
        self.change_angle_strength()
        os.chdir(self.lt_dir)
        print(os.listdir(os.path.abspath('.')))
        sb.call(["moltemplate.sh","-atomstyle","angle","system.lt"])
        sb.call(["sed","-i",'s/a\\"/a\\"\ extra\/special\/per\/atom\ 4\ extra\/bond\/per\/atom\ 2\ extra\/angle\/per\/atom\ 2/g',"system.in"])
        sb.call(["sed","-i",'s/\!\(.*\)\!/\$\{\\1\}/g',"system.in.run"])
        sb.call(["sed","-i",'s/\!(\(.*\))/\$(\\1)/g',"system.in.run"])
    
    def change_monomer_attraction(self):
        cur_path = os.path.abspath('.')
        os.chdir(self.lt_dir)
        sb.call(["sed","-i",'/\@atom\:A\ \@atom\:A/ s/twopiece\ [0-9]\?\.\?[0-9]\?/twopiece\ '+str(self.eAA)+'/g','copolyff.lt'])
        sb.call(["sed","-i",'/\@atom\:A\ \@atom\:B/ s/twopiece\ [0-9]\?\.\?[0-9]\?/twopiece\ '+str(self.eAB)+'/g','copolyff.lt'])
        sb.call(["sed","-i",'/\@atom\:B\ \@atom\:B/ s/twopiece\ [0-9]\?\.\?[0-9]\?/twopiece\ '+str(self.eBB)+'/g','copolyff.lt'])
        os.chdir(cur_path)

    def change_extent_of_reaction(self):
        cur_path = os.path.abspath('.')
        os.chdir(self.lt_dir)
        sb.call(["sed","-i",'/if\ \\"/ s/>\ \?[0-9]\?\.\?[0-9]\?[0-9]\?/>\ '+str(self.p)+'/g'])
        os.chdir(cur_path)

    def change_monomer_count(self):
        curr_dir = os.path.abspath('.')
        num_monA = int(self.total_monomers*self.monomer_A_fraction)
        num_monB = self.total_monomers-num_monA
        os.chdir(self.xyz_dir)
        sb.call(["sed","-i",'/abead/,/bbead/ s/number\ [0-9]\+/number\ '+str(num_monA)+'/g',"np.inp"])
        sb.call(["sed","-i",'/bbead/,+2 s/number\ [0-9]\+/number\ '+str(num_monB)+'/g',"np.inp"])
        os.chdir(self.lt_dir)
        sb.call(["sed",'-i','s/ABEAD\ \[[0-9]\+\]/ABEAD\ \['+str(num_monA)+'\]/g',"system.lt"])
        sb.call(["sed",'-i','s/BBEAD\ \[[0-9]\+\]/BBEAD\ \['+str(num_monB)+'\]/g',"system.lt"])
        os.chdir(curr_dir)

    def change_angle_strength(self):
        cur_path = os.path.abspath('.')
        os.chdir(self.lt_dir)
        sb.call(["sed","-i",'s/angle_coeff\(.*\)\ [0-9]\?\.\?[0-9]\?\ /angle_coeff\\1\ '+str(self.angle_strength)+'\ /g','copolyff.lt'])
        os.chdir(cur_path)

    def move_simulation_files(self,dest_dir,slurm):
        dest_folder = os.path.abspath(dest_dir+'/copoly_{}polymers_{}sequence_{}epsAA_{}epsBB_{}epsAB'.format(self.total_polymers,
                                                                                                    self.sequence,
                                                                                                    self.eAA,self.eBB,self.eAB))
        self.dest_folder = dest_folder
        if not os.path.exists(dest_folder):
            os.makedirs(dest_folder)
        for simfile in glob.glob(r''+self.lt_dir+'/system.*'):
            shutil.copy(simfile,os.path.abspath(dest_folder))
        for simfile in glob.glob(r''+self.lt_dir+'/*.txt'):
            shutil.copy(simfile,os.path.abspath(dest_folder))
        if slurm:
            shutil.copy("submit.sbatch",dest_folder)

    def move_simulation_files_remote(self,dest_folder,slurm):
        self.dest_folder = dest_folder+'/copoly_{}monomers_{}percentA_{}epsAA_{}epsBB_{}epsAB'.format(self.total_monomers,
                                                                                                 int(100*self.monomer_A_fraction),
                                                                                                    self.eAA,self.eBB,self.eAB)
        print("Moving files to directory: {}".format(self.dest_folder))
        if not self.server_connection.check_if_file_exists(self.dest_folder):
            self.server_connection.mkdir(self.dest_folder)
        for simfile in glob.glob(r''+self.lt_dir+'/system.*'):
            self.server_connection.send_file(simfile,self.dest_folder+'/'+os.path.basename(simfile))
        for simfile in glob.glob(r''+self.lt_dir+'/*.txt'):
            self.server_connection.send_file(simfile,self.dest_folder+'/'+os.path.basename(simfile))
        if slurm:
            self.server_connection.send_file(self.lt_dir+"/submit.sbatch",self.dest_folder+'/submit.sbatch')

    def analyze_simulation(self):
        print("placeholder")
    
    def start_simulation(self,slurm=True,singularity="",lmp_file="lmp"):
        os.chdir(self.dest_folder)
        singularity_path =os.path.abspath(os.path.expanduser(singularity))
        print("Singularity path is {}".format(singularity_path))
        if slurm:
            sb.call(["sbatch","submit.sbatch"])
        elif os.path.exists(singularity_path):
            sb.Popen(["singularity","run",singularity,"-i","system.in"],stdout=open('lmp_output.out','w'))
        else:
            sb.Popen([lmp_file,"-i","system.in"],stdout=open("lmp_output.out",'w'))


    def start_simulation_remote(self):
        stdin, stdout, stderr = self.server_connection.ssh_client.exec_command('cd {} \n sbatch submit.sbatch \n'.format(self.dest_folder))
        submit_status = stdout.read().decode('utf-8')
        print(submit_status)
        self.jobID = int(re.search(r'[0-9]+',submit_status).group(0))
        print(self.jobID) 






