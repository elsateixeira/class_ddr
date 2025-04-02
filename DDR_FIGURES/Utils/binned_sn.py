import numpy as np

def read_file(fn, header='#', sep=' ', skip=None):
    with open(fn, 'r') as file:
        start = True
        for iline, line in enumerate(file.readlines()):
            if skip is not None:
                if isinstance(skip, str):
                    if line.strip().startswith(skip):
                        continue
                elif iline <= skip:
                    continue
            if start:
                names = [name.strip() for name in line[len(header):].split(sep)]
                values = {name: [] for name in names}
                start = False
                continue
            line = [el for el in line.split(sep) if el]
            for name, value in zip(names, line):
                try: value = float(value)
                except ValueError: pass  # str
                values[name].append(value)
    return {name: np.array(value) for name, value in values.items()}

def build_covariance(fn):
    """Run once at the start to build the covariance matrix for the data"""
    print("Loading covariance from {}".format(fn))
    
    with open(fn, 'r') as file:
        size = int(file.readline())
    return np.loadtxt(fn, skiprows=1).reshape(size, size)

def binned_SN(df,covmat,N_bins=20,spacing='linear'):
        
    redshifts=df['zHD']
    distance_modulus=df['MU_SH0ES']
    
    # Define the number of bins and create bin edges
    num_bins = N_bins
    if spacing=='log':
        bin_edges = np.geomspace(df['zHD'].min(), df['zHD'].max(), num_bins + 1)
    else:
        bin_edges = np.linspace(df['zHD'].min(), df['zHD'].max(), num_bins + 1)
        
    # Digitize the redshifts into bins
    bin_indices = np.digitize(redshifts, bin_edges)

    # Initialize arrays to hold the mean and standard deviation for each bin
    mean_distance_modulus = np.zeros(num_bins)
    std_distance_modulus = np.zeros(num_bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Calculate mean and standard deviation for each bin
    for i in range(1, num_bins + 1):
        # Select data points and their corresponding covariance submatrix that fall into the current bin
        bin_data_indices = np.where(bin_indices == i)[0]
        bin_data = distance_modulus[bin_data_indices]
        
        if len(bin_data) > 0:
            bin_covariance_matrix = covmat[np.ix_(bin_data_indices, bin_data_indices)]
            bin_weights = np.linalg.inv(bin_covariance_matrix)
            bin_weights_sum = np.sum(bin_weights)
            
            mean_distance_modulus[i - 1] = np.sum(bin_weights @ bin_data) / bin_weights_sum
            std_distance_modulus[i - 1] = np.sqrt(1 / bin_weights_sum)
        else:
            mean_distance_modulus[i - 1] = np.nan  # If no data points in bin, set mean to NaN
            std_distance_modulus[i - 1] = np.nan  # If no data points in bin, set std to NaN

    # Return the binned data with error bars for the means and standard deviations
    return (bin_centers, mean_distance_modulus, std_distance_modulus)

if __name__=='__main__':
    
    # Load the full PantheonPlus dataset
    df=read_file('data/Pantheon+SH0ES.dat')
    
    # Trim those SN with z>0.01 used for Cosmology
    zmask=df['zHD']>0.01
    PanP={key: value[zmask] for key,value in df.items()}
    covmat=build_covariance('data/Pantheon+SH0ES_STAT+SYS.cov')[np.ix_(zmask,zmask)]
    
    z_bin,mean_mu,std_mu=binned_SN(PanP,covmat,N_bins=30,spacing='log')
    
    #Store the results
    fn='data/binned_PanP.npz'
    np.savez_compressed(fn,z_bin=z_bin,mean_mu=mean_mu,std_mu=std_mu)
    print(f'Binned dataset stored in {fn} !')
    pass