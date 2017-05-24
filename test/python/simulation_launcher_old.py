from mpi4py import MPI
import chaste
import sys
import os


def launch_job():
    

    rank = MPI.COMM_WORLD.Get_rank()
    new_comm = MPI.COMM_WORLD.Split(color=rank, key=rank)
    new_rank = new_comm.Get_rank()
    
    # cwd=os.getcwd()
    # os.mkdir(str(rank))
    # directory=os.path.join(cwd,str(rank))
    # print(rank,directory)
    # os.chdir(directory)
    
    print rank
    new_comm.Spawn("python", args=["test_launch_script.py", "-i "+str(rank)], maxprocs=2)


#     process_list = [1, 5, 7,  9, 13, 115, 3443, 4443, 23]
#     num_per_rank = int(len(process_list)/size)
#     
#     start_num = rank*num_per_rank
#     end_num = (1+rank)*num_per_rank-1
#     if rank == size-1:
#         end_num = len(process_list)-1
#     
#     for idx in range(start_num, end_num+1):
#         print "Launching Job: " + str(idx) + " on rank: " + str(rank)
#         
#         cmd = 'python test_launch_script.py -i '+ str(process_list[idx])
#         p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
#         out, err = p.communicate() 
#         result = out.split('\n')
#         for lin in result:
#             if not lin.startswith('#'):
#                 print lin.strip()
#                 
#         print "Completed Job: " + str(idx) + " on rank: " + str(rank)
                
if __name__=="__main__":
    
    work_dir = "Python/Cornea/TestLauncher/"
    file_handler = chaste.core.OutputFileHandler(work_dir, False)
    launch_job()