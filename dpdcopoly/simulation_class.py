import os
import subprocess as sb
import glob
import shutil
import server_class as svc
import re

class Simulation(object):
    def __init__(self,total_monomers=3000,monomer_A_fraction=0.5,monomer_attractions=(1,1,1),p=0.9,angle_strength=20,
                    dump_frequency=1000,lt_dir = os.path.abspath('../../lt_files/'),
                    xyz_dir = os.path.abspath('../../xyzs/'),send_to_cluster=False,servername=None):
        self.total_monomers = total_monomers
        self.monomer_A_fraction = monomer_A_fraction
        self.p = p
        self.dump_frequency = dump_frequency
        self.lt_dir = os.path.abspath(lt_dir)
        self.xyz_dir = os.path.abspath(xyz_dir)
        self.eAA,self.eBB,self.eAB = monomer_attractions
        self.angle_strength = angle_strength
        self.send_to_cluster=send_to_cluster
        if self.send_to_cluster:
            self.server_connection = svc.ServerConnection()


    def compile_simulation(self,packmol_path='packmol'):
        self.change_monomer_count()
        self.change_monomer_attraction()
        self.change_angle_strength()
        os.chdir(self.xyz_dir)
        try:
            sb.call([packmol_path],stdin=open('np.inp'))
        except OSError as error:
            print(error)
            print(("\nPackmol is not found in packmol in path."
                   "  Add packmol directory to PATH environment\n"
                    "variable or pass directory to the compile_simulation function as packmol_path argument.\n"))
            raise
        os.chdir(self.lt_dir)
        sb.call(["moltemplate.sh","-xyz",self.xyz_dir+"/np.xyz","-atomstyle","angle","system.lt"])
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
        dest_folder = os.path.abspath(dest_dir+'/copoly_{}monomers_{}percentA_{}epsAA_{}epsBB_{}epsAB'.format(self.total_monomers,
                                                                                                 int(100*self.monomer_A_fraction),
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
        if slurm:
            sb.call(["sbatch","submit.sbatch"])
        elif os.path.exists(singularity):
            sb.Popen(["singularity","run",singularity,"-i","system.in"],stdout=open('lmp_output.out','w'))
        else:
            sb.call([lmp_file,"-i","system.in"],stdout=open("lmp_output.out",'w'))


    def start_simulation_remote(self):
        stdin, stdout, stderr = self.server_connection.ssh_client.exec_command('cd {} \n sbatch submit.sbatch \n'.format(self.dest_folder))
        submit_status = stdout.read().decode('utf-8')
        print(submit_status)
        self.jobID = int(re.search(r'[0-9]+',submit_status).group(0))
        print(self.jobID) 






