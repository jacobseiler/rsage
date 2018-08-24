import numpy as np

from mpi4py import MPI

def update_cum_stats(mean_pool, std_pool, N_pool, mean_local, std_local, N_local):
    '''
    Update the cumulative statistics (such as Stellar Mass Function, Mvir-Ngamma, fesc-z) that are saved across files.
    Pooled mean formulae taken : from https://www.ncbi.nlm.nih.gov/books/NBK56512/
    Pooled variance formulae taken from : https://en.wikipedia.org/wiki/Pooled_variance

    Parameters
    ----------
    mean_pool, std_pool, N_pool : array of floats with length equal to the number of bins (e.g. the mass bins for the Stellar Mass Function).
        The current mean, standard deviation and number of data points within in each bin.  This is the array that will be updated in this function.
    mean_local, std_local, N_local : array of floats with length equal to the number of bins.
        The mean, standard deviation and number of data points within in each bin that will be added to the pool.

    Returns
    -------
    mean_pool, std_pool, N_pool : (See above)
    The updated arrays with the local values added and accounted for within the pools.

    Units
    -----
    All units are kept the same as the input units.
    Values are in real-space (not log-space).
    '''
   
    N_times_mean_local = np.multiply(N_local, mean_local)
    N_times_var_local = np.multiply(N_local - 1, np.multiply(std_local, std_local)) # Actually N - 1 because of Bessel's Correction 
                                        # https://en.wikipedia.org/wiki/Bessel%27s_correction).  #
    N_times_mean_pool = np.add(N_times_mean_local, np.multiply(N_pool, mean_pool))
    N_times_var_pool = np.add(N_times_var_local, np.multiply(N_pool - 1, np.multiply(std_pool, std_pool)))
    N_pool = np.add(N_local, N_pool)

    '''
    print(mean_local)
    print(type(mean_local))
    print((type(mean_local).__module__ == np.__name__))
    print(isinstance(mean_local, list))
    print(isinstance(mean_local,float64))
    print(isinstance(mean_local,float32))
    '''
    if (((type(mean_local).__module__ == np.__name__) == True or
(isinstance(mean_local, list) == True)) and isinstance(mean_local, float) == False and isinstance(mean_local, int) == False and isinstance(mean_local,np.float32) == False and isinstance(mean_local, np.float64) == False): # Checks to see if we are dealing with arrays. 
        for i in range(0, len(N_pool)):
            if(N_pool[i] == 0): # This case is when we have no data points in the bin. 
                mean_pool[i] = 0.0
            else:
                mean_pool[i] = N_times_mean_pool[i]/N_pool[i]
            if(N_pool[i] < 3): # In this instance we don't have enough data points to properly calculate the standard deviation.
                std_pool[i] = 0.0
            else:
                std_pool[i] = np.sqrt(N_times_var_pool[i]/ (N_pool[i] - 2)) # We have -2 because there is two instances of N_pool contains two 'N - 1' terms. 
        
    else:
        mean_pool = N_times_mean_pool / N_pool

        if(N_pool < 3):
            std_pool = 0.0
        else:
            std_pool = np.sqrt(N_times_var_pool / (N_pool - 2))
 
    return mean_pool, std_pool


def collect_hist_across_tasks(rank, comm, hists):

    master_hists = []

    for model_number in range(len(hists)):

        if rank == 0:
            counts_total = np.zeros_like(hists[model_number])
        else:
            counts_total = None

        comm.Reduce([hists[model_number], MPI.FLOAT], [counts_total, MPI.FLOAT], op = MPI.SUM, root = 0) # Sum all the stellar mass and pass to Rank 0.

        master_hists.append(counts_total)

    return master_hists


def collect_across_tasks(rank, comm, mean_per_task, std_per_task, N_per_task, SnapList,
                         BinSnapList=[], binned=False, m_bin_low=0.0, 
                         m_bin_high=0.0, my_bin_width=0.1):
                
    """
    Reduces arrays that are unique to each task onto the master task.

    The dimensions of the input arrays will change slightly if we are collecting a statistics
    that is binned across e.g., halo mass or galaxy stellar mass.

    Parameters
    ----------

    rank: Integer
        Rank of the task.

    mean_per_task, std_per_task, N_per_task: Nested 2D (or 3D if binned == True) arrays of floats.  
                                             Outer length is equal to the number of models.
                                             Inner length is equal to the number of snapshots the data has been calculated for.
                                             Most inner length is equal to the number of bins.
        Contains the mean/standard deviation/number of objects unique for each task.

    SnapList: Nested 2D arrays of integers.  Outer length is equal to the number of models.
        Contains the snapshot numbers the data has been calculated for each model. 

    BinSnapList: Nested 2D arrays of integers. Outer length is equal to the number of models.
        Often statistics are calculated for ALL snapshots but we only wish to plot for a subset of snapshots.
        This variable allows the binned data to be collected for only a subset of the snapshots.

    binned: Boolean.
        Dictates whether the collected data is a 2D or 3D array with the inner-most array being binned across e.g., halo mass.

    Returns
    ----------

    master_mean, master_std, master_N: Nested 2D (or 3D if binned == True) arrays of floats.
                                       Shape is identical to the input mean_per_task etc.
        If rank == 0 these contain the collected statistics.
        Otherwise these will be none.

    master_bin_middle: Array of floats.
        Contains the location of the middle of the bins for the data.         
    """


    master_mean = []
    master_std = []
    master_N = []

    master_bin_middle = []

    for model_number in range(0, len(SnapList)): 

        master_mean.append([])
        master_std.append([])
        master_N.append([])

        master_bin_middle.append([])

        # If we're collecting a binned statistic (e.g., binned across halo mass), then we need to perform the collecting per snapshot.
        if binned:
            count = 0 
            for snapshot_idx in range(len(SnapList[model_number])):
                if SnapList[model_number][snapshot_idx] == BinSnapList[model_number][count]:
                    master_mean[model_number], master_std[model_number], master_N[model_number] = calculate_pooled_stats(rank, comm, master_mean[model_number], master_std[model_number], master_N[model_number], mean_per_task[model_number][snapshot_idx], std_per_task[model_number][snapshot_idx], N_per_task[model_number][snapshot_idx])
                    master_bin_middle[model_number].append(np.arange(m_bin_low,
                                                                     m_bin_high+my_bin_width, 
                                                                     my_bin_width)[:-1] 
                                                           + my_bin_width* 0.5)

                    count += 1

                    if count == len(BinSnapList[model_number]):
                        break

        else:
            master_mean[model_number], master_std[model_number], master_N[model_number] = calculate_pooled_stats(rank, comm, master_mean[model_number], master_std[model_number], master_N[model_number], 
                                                                                                                 mean_per_task[model_number], std_per_task[model_number], 
                                                                                                                 N_per_task[model_number])

            if rank == 0:
                master_mean[model_number] = master_mean[model_number][0]
                master_std[model_number] = master_std[model_number][0]
                master_N[model_number] = master_N[model_number][0]
    return master_mean, master_std, master_N, master_bin_middle


def calculate_pooled_stats(rank, comm, mean_pool, std_pool, N_pool, mean_local, std_local, N_local):
    '''
    Calculates the pooled mean and standard deviation from multiple processors and appends it to an input array.
    Formulae taken from https://en.wikipedia.org/wiki/Pooled_variance
    As we only care about these stats on the rank 0 process, we make use of junk inputs/outputs for other ranks.

    NOTE: Since the input data may be an array (e.g. pooling the mean/std for a stellar mass function).

    Parameters
    ----------

    rank: Integer
        Rank of the task.

    mean_pool, std_pool, N_pool : array of floats.
        Arrays that contain the current pooled means/standard deviation/number of data points (for rank 0) or just a junk input (for other ranks).
    mean_local, mean_std : float or array of floats.
        The non-pooled mean and standard deviation unique for each process.
    N_local : floating point number or array of floating point numbers. 
        Number of data points used to calculate the mean/standard deviation that is going to be added to the pool.
        NOTE: Use floating point here so we can use MPI.DOUBLE for all MPI functions.

    Returns
    -------
    mean_pool, std_pool : array of floats.
        Original array with the new pooled mean/standard deviation appended (for rank 0) or the new pooled mean/standard deviation only (for other ranks).

    Units
    -----
    All units are the same as the input.
    All inputs MUST BE real-space (not log-space).
    '''

    if isinstance(mean_local, list) == True:    
        if len(mean_local) != len(std_local):
            print("len(mean_local) = {0} \t len(std_local) = {1}".format(len(mean_local), len(std_local)))
            raise ValueError("Lengths of mean_local and std_local should be equal")
   
    if ((type(mean_local).__module__ == np.__name__) == True or (isinstance(mean_local, list) == True)): # Checks to see if we are dealing with arrays. 
    
        N_times_mean_local = np.multiply(N_local, mean_local)
        N_times_var_local = np.multiply(N_local, np.multiply(std_local, std_local))
        
        N_local = np.array(N_local).astype(float)
        N_times_mean_local = np.array(N_times_mean_local).astype(np.float32)

        if rank == 0: # Only rank 0 holds the final arrays so only it requires proper definitions.
            N_times_mean_pool = np.zeros_like(N_times_mean_local) 
            N_pool_function = np.zeros_like(N_local)
            N_times_var_pool = np.zeros_like(N_times_var_local)

            N_times_mean_pool = N_times_mean_pool.astype(np.float64) # Recast everything to double precision then use MPI.DOUBLE.
            N_pool_function = N_pool_function.astype(np.float64)
            N_times_var_pool = N_times_var_pool.astype(np.float64)
        else:
            N_times_mean_pool = None
            N_pool_function = None
            N_times_var_pool = None


        N_times_mean_local = N_times_mean_local.astype(np.float64)
        N_local = N_local.astype(np.float64)
        N_times_var_local = N_times_var_local.astype(np.float64)

        comm.Reduce([N_times_mean_local, MPI.DOUBLE], [N_times_mean_pool, MPI.DOUBLE], op = MPI.SUM, root = 0) # Sum the arrays across processors.
        comm.Reduce([N_local, MPI.DOUBLE],[N_pool_function, MPI.DOUBLE], op = MPI.SUM, root = 0)   
        comm.Reduce([N_times_var_local, MPI.DOUBLE], [N_times_var_pool, MPI.DOUBLE], op = MPI.SUM, root = 0)
        
    else:
    
        N_times_mean_local = N_local * mean_local
        N_times_var_local = N_local * std_local * std_local

        N_times_mean_pool = comm.reduce(N_times_mean_local, op = MPI.SUM, root = 0)
        N_pool_function = comm.reduce(N_local, op = MPI.SUM, root = 0)
        N_times_var_pool = comm.reduce(N_times_var_local, op = MPI.SUM, root = 0)
    
    if rank == 0:

        mean_pool_function = np.zeros((len(N_pool_function)))
        std_pool_function = np.zeros((len(N_pool_function)))

        for i in range(0, len(N_pool_function)):
            if N_pool_function[i] == 0:
                mean_pool_function[i] = 0.0
            else:
                mean_pool_function[i] = np.divide(N_times_mean_pool[i], N_pool_function[i])
            if N_pool_function[i] < 3:
                std_pool_function[i] = 0.0
            else:
                std_pool_function[i] = np.sqrt(np.divide(N_times_var_pool[i], N_pool_function[i]))
       
        mean_pool.append(mean_pool_function)
        std_pool.append(std_pool_function)
        N_pool.append(N_pool_function)

        return mean_pool, std_pool, N_pool
    else:
    
        return mean_pool, std_pool, N_pool_function # Junk return because non-rank 0 doesn't care.

