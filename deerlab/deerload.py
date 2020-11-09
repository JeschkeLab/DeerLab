from warnings import warn
import re
import numpy as np
import os
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------
def deerload(fullbasename, plot=False, full_output=False, *args,**kwargs):
    r"""
    Load file in BES3T format (Bruker EPR Standard for Spectrum Storage and Transfer)
    
    * .DSC: description file
    * .DTA: data file
    
    used on Bruker ELEXSYS and EMX machines. This function can handle experiments acquired 
    in a multi-dimensional setup. 

    Parameters
    ----------
    fullbasename : string
        Full name of data file.
    
    full_output : boolean, optional
        Return the parameter file entries as a third output. Disabled by default.
    
    plot : boolean, optional
        Display a plot of the real and imaginary parts of the loaded data. Disabled by default.

    Returns
    -------
    t : ndarray
        Time axis in microseconds.Its structure depends on the dimensionality of the experimental datasets:

        * 1D-datasets: the time axis of ``N`` points is returned as a one-dimensional ndarray of shape ``(N,)``
        * 2D-datasets: the ``M`` time-axes of ``N`` points are returned as a two-dimensional ndarray of shape ``(N,M)``. The i-th axis can be accessed via ``t[:,i]``.
        
    V : ndarray
        Experimental signal(s). Its structure depends on the dimensionality of the experimental datasets:

        * 1D-datasets: the signal of ``N`` points is returned as a one-dimensional ndarray of shape ``(N,)``
        * 2D-datasets: the ``M`` signals of ``N`` points are returned as a two-dimensional ndarray of shape ``(N,M)``. The i-th  signal can be accessed via ``V[:,i]``.
        
    pars : dict
        Parameter file entries, returned if ``full_output`` is ``True``.

    Notes
    -----
    Code based on BES3T version 1.2 (Xepr >=2.1).
    """
    filename = fullbasename[:-4]
    fileextension = fullbasename[-4:].upper() # case insensitive extension
    if fileextension in ['.DSC','.DTA']:
        filename_dsc = filename + '.DSC'
        filename_dta = filename + '.DTA'
    else:
        raise ValueError("Only Bruker BES3T files with extensions .DSC or .DTA are supported.")
    
    if not os.path.exists(filename_dta):
        filename_dta = filename_dta[:-4] + filename_dta[-4:].lower()
        filename_dsc = filename_dsc[:-4] + filename_dsc[-4:].lower()
    
    # Read descriptor file (contains key-value pairs)
    parameters = read_description_file(filename_dsc)
    parDESC = parameters["DESC"]
    parSPL = parameters["SPL"]
    
    # XPTS, YPTS, ZPTS specify the number of data points along x, y and z.
    if 'XPTS' in parDESC:
        nx = int(parDESC['XPTS'])
    else:
        raise ValueError('No XPTS in DSC file.')
    
    if 'YPTS' in parDESC:
        ny = int(parDESC['YPTS'])
    else:
        ny = 1
    if 'ZPTS' in parDESC:
        nz = int(parDESC['ZPTS'])
    else:
        nz = 1
    
    # BSEQ: Byte Sequence
    # BSEQ describes the byte order of the data, big-endian (BIG, encoding = '>') or little-endian (LIT, encoding = '<').
    if 'BSEQ' in parDESC:
        if 'BIG' == parDESC['BSEQ']:
            byteorder = '>' 
        elif 'LIT' == parDESC['BSEQ']:
            byteorder = '<'
        else:
            raise ValueError('Unknown value for keyword BSEQ in .DSC file!')
    else:
        warn('Keyword BSEQ not found in .DSC file! Assuming BSEQ=BIG.')
        byteorder = '>'
    
    # IRFMT: Item Real Format
    # IIFMT: Item Imaginary Format
    # Data format tag of BES3T is IRFMT for the real part and IIFMT for the imaginary part.
    if 'IRFMT' in parDESC:
        IRFTM = parDESC["IRFMT"]
        if 'C' == IRFTM:
            dt_spc = np.dtype('int8')
        elif 'S' == IRFTM:
            dt_spc = np.dtype('int16')
        elif 'I' == IRFTM:
            dt_spc = np.dtype('int32')
        elif 'F' == IRFTM:
            dt_spc = np.dtype('float32')
        elif 'D' == IRFTM:
            dt_spc = np.dtype('float64')
        elif 'A' == IRFTM:
            raise TypeError('Cannot read BES3T data in ASCII format!')
        elif ('0' or 'N') == IRFTM:
            raise ValueError('No BES3T data!')
        else:
            raise ValueError('Unknown value for keyword IRFMT in .DSC file!')
    else:
        raise ValueError('Keyword IRFMT not found in .DSC file!')
    
    # IRFMT and IIFMT must be identical.
    if "IIFMT" in parDESC:
        if  parDESC["IIFMT"] != parDESC["IRFMT"]:
            raise ValueError("IRFMT and IIFMT in DSC file must be identical.")
    
    # Preallocation of the abscissa
    maxlen = max(nx,ny,nz)
    abscissa = np.full((maxlen,3),np.nan)
    # Construct abscissa vectors
    AxisNames = ['X','Y','Z']
    Dimensions = [nx,ny,nz]
    for a in AxisNames:
        index = AxisNames.index(a)
        axisname = a+'TYP'
        axistype = parDESC[axisname]
        if Dimensions[index] > 1:
            pass
        else:
            if 'IGD'== axistype:
                # Nonlinear axis -> Try to read companion file (.XGF, .YGF, .ZGF)
                companionfilename=str(filename+'.'+a+'GF')
                if 'D' == parDESC[str(a+'FMT')]:
                    dt_axis = np.dtype('float64')
                elif 'F' == parDESC[str(a+'FMT')]:
                    dt_axis = np.dtype('float32')
                elif 'I' == parDESC[str(a+'FMT')]:
                    dt_axis = np.dtype('int32')
                elif 'S' == parDESC[str(a+'FMT')]:
                    dt_axis = np.dtype('int16')
                else:
                    raise ValueError('Cannot read data format {0} for companion file {1}'.format(str(a+'FMT'),companionfilename))

                dt_axis = dt_axis.newbyteorder(byteorder)
                # Open and read companion file
                with open(companionfilename,'rb') as fp:
                    if fp > 0:
                        abscissa[:Dimensions[index],index] = np.frombuffer(fp.read(),dtype=dt_axis)
                    else:
                        warn('Could not read companion file {0} for nonlinear axis. Assuming linear axis.'.format(companionfilename))
                axistype='IDX'
        if axistype == 'IDX':
            minimum = float(parDESC[str(a+'MIN')])
            width = float(parDESC[str(a+'WID')])
            npts = int(parDESC[str(a+'PTS')])
            if width == 0:
                warn('Warning: {0} range has zero width.\n'.format(a))
                minimum = 1.0
                width = len(a) - 1.0
            abscissa[:Dimensions[index],index] = np.linspace(minimum,minimum+width,npts)
        if axistype == 'NTUP':
            raise ValueError('Cannot read data with NTUP axes.')
    
    # In case of column filled with NaN, erase the column in the array
    abscissa = abscissa[:,~np.isnan(abscissa).all(axis=0)]
    dt_data = dt_spc
    dt_spc = dt_spc.newbyteorder(byteorder)

    # Read data matrix and separate complex case from real case.
    data = np.full((nx,ny,nz),np.nan)
    # reorganize the data in a "complex" way as the real part and the imaginary part are separated
    # IKKF: Complex-data Flag
    # CPLX indicates complex data, REAL indicates real data.
    if 'IKKF' in parDESC:
        if parDESC['IKKF'] == 'REAL':
            data = np.full((nx,ny,nz),np.nan) 
            with open(filename_dta,'rb') as fp:
                 data = np.frombuffer(fp.read(),dtype=dt_spc).reshape(nx,ny,nz)
            data = np.copy(data)
        elif parDESC['IKKF'] == 'CPLX':
            dt_new = np.dtype('complex')
            data = np.full((2*nx*ny*nz,1),np.nan)
            with open(filename_dta,'rb') as fp:
                data = np.frombuffer(fp.read(),dtype=dt_spc)
                # Check if there is multiple harmonics (High field ESR quadrature detection)
                harmonics = np.array([[False] * 5]*2) # outer dimension for the 90 degree phase
                for j,jval in enumerate(['1st','2nd','3rd','4th','5th']):
                    for k,kval in enumerate(['','90']):
                        thiskey = 'Enable'+jval+'Harm'+kval
                        if thiskey in parameters.keys() and parameters[thiskey]:
                            harmonics[k,j] = True
                n_harmonics = sum(harmonics)[0]
                if n_harmonics != 0:
                    ny = int(len(data)/nx/n_harmonics)
            
            # copy the data to a writable numpy array
            data = np.copy(data.astype(dtype=dt_data).view(dtype=dt_new).reshape(nx,ny,nz))
        else:
            raise ValueError("Unknown value for keyword IKKF in .DSC file!")
    else:
        warn("Keyword IKKF not found in .DSC file! Assuming IKKF=REAL.")
    
    if nz == 1:
        data = data.reshape(nx,ny)
            
     # ns -> us converesion
    abscissa /= 1e3

    # Ensue proper numpy formatting
    abscissa,data = np.atleast_1d(abscissa,data)
    abscissa = np.squeeze(abscissa)
    data = np.squeeze(data)

    if plot:
        plt.plot(abscissa,np.real(data),abscissa,np.imag(data))
        plt.xlabel("Time (μs)")
        plt.ylabel("Intensity [a.u.]")
        plt.grid(alpha=0.3)
        plt.show()
    
    if full_output:
        return abscissa, data, parameters
    else:
        return abscissa, data


def read_description_file(DSCFileName):
    """
    Parameters = readDSCfile(DSCFileName)
    Reads a Bruker BES3T .DSC file DSCFileName and returns a dictionary in Parameters.
    """
    with open(DSCFileName,'r',encoding='utf-8',errors='ignore') as f:
        allLines = f.readlines()
    
    # Remove lines with comments
    allLines = [l for l in allLines if not l.startswith("*")]

    # Remove newlines
    allLines = [l.rstrip("\n") for l in allLines]
    
    # Remove empty lines
    allLines = [l for l in allLines if l]
    
    # Merge any line ending in \n\ with the subsequent one
    allLines2 = []
    val = ""
    for line in allLines:
        val = "".join([val, line])    
        if val.endswith("\\"):
            val = val.strip("\\")
        else:
            allLines2.append(val)
            val = ""
    allLines = allLines2
    
    Parameters = {}
    SectionName = None
    DeviceName = None
    
    # Regular expressions to match layer/section headers, device block headers, and key-value lines
    reSectionHeader = re.compile(r"#(\w+)\W+(\d+.\d+)")
    reDeviceHeader = re.compile(r"\.DVC\W+(\w+),\W+(\d+\.\d+)")
    reKeyValue = re.compile(r"(\w+)\W+(.*)")
    
    for line in allLines:
        
        # Layer/section header (possible values: #DESC, #SPL, #DSL, #MHL)
        mo1 = reSectionHeader.search(line) 
        if mo1:
            SectionName = mo1.group(1)
            SectionVersion = mo1.group(2)
            if SectionName not in {"DESC","SPL","DSL","MHL"}:
                raise ValueError("Found unrecognized section " + SectionName + ".")
            Parameters[SectionName] = {"_version": SectionVersion}
            DeviceName = None
            continue
        
        # Device block header (starts with .DVC)
        mo2 = reDeviceHeader.search(line)
        if mo2:
            DeviceName = mo2.group(1)
            DeviceVersion = mo2.group(2)
            Parameters[SectionName][DeviceName] = {"_version": DeviceVersion}
            continue
        
        # Key/value entry
        mo3 = reKeyValue.search(line)
        if not mo3:
            raise ValueError("Key/value pair expected.")        
        if not SectionName:
            raise ValueError("Found a line with key/value pair outside any layer.")
        if SectionName=="DSL" and not DeviceName:
            raise ValueError("Found a line with key-value pair outside .DVC in #DSL layer.")
        
        Key = mo3.group(1)
        Value = mo3.group(2)
        if DeviceName:
            Parameters[SectionName][DeviceName][Key] = Value
        else:
            Parameters[SectionName][Key] = Value
    
        # Assert DESC section is present
        if "DESC" not in Parameters:
            raise ValueError("Missing DESC section in .DSC file.")
    
    return Parameters
