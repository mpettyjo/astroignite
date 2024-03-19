#!/usr/bin/env python
# coding: utf-8

# In[68]:

import re
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------------

# readiqmol function: Function to parse data from IQMol output files
def readiqmol(filename, sflow = 0.970, sfmid = 0.953, sfhigh = 0.9475):
    iq = {"Frequency":[], "Intensity":[], "nomodes":0, "units":'Wavenumber'}
    with open(filename, 'r') as iqmol_out:
        for line in iqmol_out:
            if "Mode:" in line:
                pass
                #print(line.split())
            if "Frequency:" in line:
                #print(line.split())
                iq["Frequency"].append(float(line.split()[1]))
                try:
                    iq["Frequency"].append(float(line.split()[2]))
                    iq["Frequency"].append(float(line.split()[3]))
                except:
                    pass
            if "IR Intens:" in line:
                #print(line.split())
                iq["Intensity"].append(float(line.split()[2]))
                try:
                    iq["Intensity"].append(float(line.split()[3]))
                    iq["Intensity"].append(float(line.split()[4]))
                except:
                    pass
                    
    freqs_scaled = []
    for frequency in iq["Frequency"]:
        if frequency < 1000:
            freqs_scaled.append(frequency*sflow)
        if 1000 <= frequency < 2000:
            freqs_scaled.append(frequency*sfmid)
        if frequency >= 2000:
            freqs_scaled.append(frequency*sfhigh)
            
    iq["Freq_Scaled"] = freqs_scaled
    iq["nomodes"] += len(iq["Frequency"])
    iq["Intensity_norm"] = normintens(iq["Intensity"]).tolist()
    
    iq["Frequency_um"] = np.zeros(len(iq["Freq_Scaled"])).tolist()
    for i in range(iq["nomodes"]):
        iq["Frequency_um"][i] = 10000/iq["Freq_Scaled"][i]
    

    return iq
    
# --------------------------------------------------------------------------------

# normintens function: Function to normalised the intensities from the IQMol output file
def normintens(intensities):
    norm = np.max(intensities)
    normintensities = np.zeros(len(intensities))
    for i in range(len(intensities)):
        normintensities[i] = intensities[i]/norm
    return normintensities

# --------------------------------------------------------------------------------

# makespectra function: Function to visualise the spectra from IQMol output files
def makespectra(mol, units="CM",  npts = 100, peakwidth = 20):

    if units=="CM":
        freq = mol['Freq_Scaled']
        intens = mol['Intensity']
    elif units =="MICRO":
        freq = mol['Frequency_um']
        intens = mol['Intensity'] 
        peakwidth = peakwidth/200.0
            
    minf = min(freq)-(min(freq)*0.1)
    maxf = max(freq)+(max(freq)*0.1)
    x = np.linspace(minf,maxf,npts)
    y = np.zeros(npts)
    
    for i in range(len(freq)):
        #This add each peak sequentially based on the intensity of the peak and the 
        y = y - np.exp(-(2.77/(2*peakwidth)**2)*(x-freq[i])**2)*intens[i]
    
    ymin = min(y)
    ynorm= []
    for i in range(len(y)):
        ynorm.append(y[i]/ymin)
    mol['Frequency_width'] = x.tolist()
    mol['Intensity_width'] = np.array(ynorm).tolist()
    return mol
    
    
def givewidth(freq, intens, units="CM",  npts = 100, peakwidth = 20):
    if units=="CM":
        pass
    elif units =="MICRO":
        peakwidth = peakwidth/200.0
    
    minf = min(freq)-(min(freq)*0.1)
    maxf = max(freq)+(max(freq)*0.1)
    x = np.linspace(minf,maxf,npts)
    y = np.zeros(npts)
    
    for i in range(len(freq)):
        #This add each peak sequentially based on the intensity of the peak and the 
        y = y - np.exp(-(2.77/(2*peakwidth)**2)*(x-freq[i])**2)*intens[i]
    
    ymin = min(y)
    ynorm= []
    for i in range(len(y)):
        ynorm.append(y[i]/ymin)
        
    x = x.tolist()
    y = (-y).tolist()
    return x, y

# --------------------------------------------------------------------------------

# readNIST function: function to parse data from the NIST txt files
def readNIST(NISTfile, **kwargs):
    NIST = {"Frequency":[], "Intensity":[], "Molecule":'', "units":''}
    with open(NISTfile,'r')as NISTf:
        cnt = 0
        for line in NISTf: 
            cnt +=1
            if "TITLE" in line:
                NIST["Molecule"] = str(line.split('=')[-1])
                #print(molecule)
            if "XUNITS" in line:
                NIST["units"] = line.split('=')[-1]
                print ("NIST units are in:", NIST["units"])
            if cnt >10 and "#" not in line:
                NIST["Intensity"].append(1-float(line.split()[1]))
                NIST["Frequency"].append(float(line.split()[0]))

    Nintens = list(range(0, len(NIST["Intensity"])))
    minfreq = min(NIST["Frequency"])
    maxfreq = max(NIST["Frequency"])
    npts = len(NIST["Frequency"])
    maxintens = max(NIST["Intensity"])
    minintens = min(NIST["Intensity"])
#     print(maxintens)
    
    #Normailse
    if abs(minintens)> abs(maxintens):
        for i in range(len(NIST["Intensity"])):
            Nintens[i]= NIST["Intensity"][i]/minintens
    else:
        for i in range(len(NIST["Intensity"])):
            Nintens[i]= NIST["Intensity"][i]/maxintens
    
    NIST["Intensity_norm"] = np.array(Nintens).tolist()
    
    
    if "MICROMETERS" in NIST["units"]:
        print(f"Units in micrometers - use peakwidth of {1/200} and change units of comparison stick spectra")
#         for i in range(calnomodes):
#             calfreq[i] = 10000/calfreq[i]
#         #change pw 
#         peakwidth = peakwidth/200
#         #print(calfreq)
    return NIST

# --------------------------------------------------------------------------------

# figsettings function: Function providing some basic settings for figure creation in python
def figsettings(ax1, units = "CM"):
#     def cmtoum(x):
#         return (x**-1)/(10**-4)
#     def umtocm(x):
#         return (x**-1)/(10**-4)
    def tick_function(X):
        V = 1.0000000e4/(X)
        return ["%.2f" % z for z in V]
    sns.set_context('talk', font_scale=1.1)
    
    
    ax2 = ax1.twiny()
    if "CM" in units:
    	ax1.set_xlabel('Wavenumber (cm$^{-1}$) \n', labelpad=15)
    	ax2.set_xlabel("\n Wavelength ($\mu m$)", labelpad=15)
    elif "MICRO" in units:
        ax1.set_xlabel(r"Wavelength / $\mu m$", labelpad=15)
        ax2.set_xlabel(r'Wavenumber / cm$^{-1}$', labelpad=15)
    ax1.set_ylabel('Absorbance (arb. units)', labelpad=15)
    ax1.set_xlim(0,)
#     ax1.set_ylim(0,)
    ax1Ticks = ax1.get_xticks()
    ax2.set_xticks(ax1Ticks)
    ax2.set_xbound(ax1.get_xbound())
    ax2.set_xticklabels(tick_function(ax1Ticks))
    
    return ax1, ax2

# --------------------------------------------------------------------------------

# GWP function: Function to calculate global warming potentials using IQMol data
def GWP_calculator(mol_data, molecule, lifetime, th):
    
    from scipy import integrate
    
    periodic_table = {
    'H': 1.00794 ,'He': 4.002602,'Li': 6.941,'Be': 9.012182,'B': 10.811,'C': 12.0107,'N': 14.0067,'O': 15.9994,
    'F': 18.9984032,'Ne': 20.1797,'Na': 22.98977,'Mg': 24.305,'Al': 26.981538,'Si': 28.0855,'P': 30.973761,'S': 32.065,
    'Cl': 35.453,'Ar': 39.948,'K': 39.0983,'Ca': 40.078,'Sc': 44.95591,'Ti': 47.867,'V': 50.9415,'Cr': 51.9961,'Mn': 54.938049,
    'Fe': 55.845,'Co': 58.9332,'Ni': 58.6934,'Cu': 63.546,'Zn': 65.409,'Ga': 69.723,'Ge': 72.64,'As': 74.9216,
    'Se': 78.96,'Br': 79.904,'Kr': 83.798,'Rb': 85.4678,'Sr': 87.62,'Y': 88.90585,'Zr': 91.224,'Nb': 92.90638,
    'Mo': 95.94,'Tc': 98,'Ru': 101.07,'Rh': 102.9055,'Pd': 106.42,'Ag': 107.8682,'Cd': 112.411,'In': 114.818,
    'Sn': 118.71,'Sb': 121.76,'Te': 127.6,'I': 126.90447,'Xe': 131.293,'Cs': 132.90545,'Ba': 137.327,'La': 138.9055,
    'Ce': 140.116,'Pr': 140.90765,'Nd': 144.24,'Pm': 145,'Sm': 150.36,'Eu': 151.964,'Gd': 157.25,'Tb': 158.92534,
    'Dy': 162.5,'Ho': 164.93032,'Er': 167.259,'Tm': 168.93421,'Yb': 173.04,'Lu': 174.967,'Hf': 178.49,'Ta': 180.9479,
    'W': 183.84,'Re': 186.207,'Os': 190.23,'Ir': 192.217,'Pt': 195.078,'Au': 196.96655,'Hg': 200.59,'Tl': 204.3833,
    'Pb': 207.2,'Bi': 208.98038,'Po': 209,'At': 210,'Rn': 222,'Fr': 223,'Ra': 226,'Ac': 227,
    'Th': 232.0381,'Pa': 231.03588,'U': 238.02891,'Np': 237,'Pu': 244,'Am': 243,'Cm': 247,'Bk': 247,
    'Cf': 251,'Es': 252,'Fm': 257,'Md': 258,'No': 259,'Lr': 262,'Rf': 261,'Db': 262,
    'Sg': 266,'Bh': 264,'Hs': 277,'Mt': 268,'Ds': 281,'Rg': 272,'Cn': 285,'Nh': 286,
    'Fl': 289,'Mc': 289,'Lv': 293,'Ts': 294,'Og': 294,
    }
    
    freq_wd = [
    5,15,25,35,45,55,65,75,85,95,105,115,125,135,145,155,165,175,185,195,
    205,215,225,235,245,255,265,275,285,295,305,315,325,335,345,355,365,375,385,395,
    405,415,425,435,445,455,465,475,485,495,505,515,525,535,545,555,565,575,585,595,
    605,615,625,635,645,655,665,675,685,695,705,715,725,735,745,755,765,775,785,795,
    805,815,825,835,845,855,865,875,885,895,905,915,925,935,945,955,965,975,985,995,
    1005,1015,1025,1035,1045,1055,1065,1075,1085,1095,1105,1115,1125,1135,1145,1155,1165,1175,1185,1195,
    1205,1215,1225,1235,1245,1255,1265,1275,1285,1295,1305,1315,1325,1335,1345,1355,1365,1375,1385,1395,
    1405,1415,1425,1435,1445,1455,1465,1475,1485,1495,1505,1515,1525,1535,1545,1555,1565,1575,1585,1595,
    1605,1615,1625,1635,1645,1655,1665,1675,1685,1695,1705,1715,1725,1735,1745,1755,1765,1775,1785,1795,
    1805,1815,1825,1835,1845,1855,1865,1875,1885,1895,1905,1915,1925,1935,1945,1955,1965,1975,1985,1995,
    2005,2015,2025,2035,2045,2055,2065,2075,2085,2095,2105,2115,2125,2135,2145,2155,2165,2175,2185,2195,
    2205,2215,2225,2235,2245,2255,2265,2275,2285,2295,2305,2315,2325,2335,2345,2355,2365,2375,2385,2395,
    2405,2415,2425,2435,2445,2455,2465,2475,2485,2495]
    
    rf_pfw = [
    0.001,0.0056,0.0145,0.0213,0.0363,0.0421,0.0676,0.0724,0.0978,0.112,0.131,0.191,0.174,0.202,0.255,0.233,0.348,0.287,0.388,0.439,
    0.337,0.484,0.42,0.584,0.514,0.536,0.665,0.611,0.62,0.776,0.678,0.825,0.732,0.856,0.862,0.836,1.02,0.99,1.06,1,
    1.28,1.08,1.33,1.41,1.34,1.26,1.65,1.46,1.68,1.88,1.67,1.67,1.97,2.2,1.96,2.21,2.01,2.13,1.61,1.19,
    1.36,0.54,0.864,0.365,0.0777,0.0573,0.0003,0.0496,0.0777,0.222,0.51,0.806,0.607,1.32,1.44,2.37,2.95,3.03,3.15,2.74,
    3.09,3.2,3.19,3.23,3.16,3.05,3.17,3.1,3.04,3.1,3.01,3.03,2.93,2.95,2.81,2.81,2.81,2.67,2.72,2.61,
    2.31,1.82,1.47,1.37,1.51,1.11,1.78,2.18,2.22,2.18,2.1,2.07,1.99,1.79,1.88,1.79,1.72,1.52,1.52,1.52,
    1.67,1.23,1.15,1.26,0.906,0.815,0.418,0.455,0.411,0.359,0.207,0.296,0.371,0.237,0.231,0.28,0.173,0.18,0.181,0.122,
    0.182,0.098,0.12,0.0907,0.122,0.0569,0.0917,0.0613,0.0658,0.0528,0.0372,0.0462,0.0341,0.0335,0.0325,0.0281,0.0337,0.042,0.0739,0.0635,0.0518,
    0.0401,0.0346,0.0274,0.0261,0.019,0.0242,0.0233,0.0194,0.0167,0.0182,0.0179,0.0353,0.0156,0.0226,0.0265,0.0272,0.02,0.0308,0.0199,0.0402,
    0.04,0.03,0.0315,0.0228,0.0569,0.0255,0.066,0.0439,0.0561,0.0452,0.0404,0.0496,0.0735,0.036,0.077,0.0581,0.106,0.081,0.0466,0.131,
    0.0708,0.115,0.133,0.0767,0.106,0.0582,0.0543,0.0953,0.0752,0.101,0.0883,0.0912,0.086,0.0985,0.098,0.0863,0.0789,0.0494,0.0309,0.0202,
    0.0176,0.0204,0.0078,0.0053,0.003,0.0016,0.0007,0.0006,0.0001,0.0001,0,0,0,0,0,0,0,0.0027,0.0451,0.0473,0.0455,0.0436,0.0412,0.0385,0.0373,0.0367,0.0338,0.0339,0.0333]

    calc_rf = []
    for i in np.arange(0,len(freq_wd)):
        calc_rf.append([0])
    
    # Calculate radiative transfer for species i
    # -------------------------------------------

    freqs_considered = []
    ints_considered = []
    for freq,ints in zip(mol_data['Freq_Scaled'],mol_data['Intensity']):
        if freq <= 2500:
            freqs_considered.append(freq)
            ints_considered.append(ints)

    for location in np.arange(0,len(freq_wd)):
        for freq,ints in zip(freqs_considered,ints_considered):
            if freq_wd[location]-5 < freq < 5+freq_wd[location]:
                calc_rf[location].append((rf_pfw[location]/(10**-15))*ints*(1/(6.023*10**23))*1000*100)
                
    ai = 0
    for i in calc_rf:
        value = sum(i)
        ai = ai + value
    
    # Calculate molecular mass for species i
    # ---------------------------------------
    mass = 0
    for elements in re.findall('[A-Z][^A-Z]*', molecule.split('_')[0].rstrip()):
        elements_list = re.split('(\d+)', elements.split()[0])

        if len(elements_list) == 1:
            mass = mass + periodic_table[elements_list[0]]

        else:
            elements_list.pop(-1)
            mass = mass + float(elements_list[-1])*(periodic_table[elements_list[0]])
    
    # Calculate GWP for species i
    # -----------------------------
    AGWP_CO2 = {20:0.192, 50:0.676, 100:2.223}
    
    func = lambda t: np.exp(-t/lifetime)
    GWP = ai*(integrate.quad(func,0,th)[0]/AGWP_CO2[th])*(44.0095/mass)*1000
    
    return print('Computed radiative forcing (ai) for',molecule,':',ai,'W/(m^2*ppb)','\nGWP',GWP,'kgCO2/kg'+molecule)

# --------------------------------------------------------------------------------

# GWP function with k: Function to calculate global warming potentials using IQMol data for molecules with no lifetimes available
def GWP_kconstant(mol_data, molecule, th, A, ER_factor):

    from scipy import integrate
    
    periodic_table = {
    'H': 1.00794 ,'He': 4.002602,'Li': 6.941,'Be': 9.012182,'B': 10.811,'C': 12.0107,'N': 14.0067,'O': 15.9994,
    'F': 18.9984032,'Ne': 20.1797,'Na': 22.98977,'Mg': 24.305,'Al': 26.981538,'Si': 28.0855,'P': 30.973761,'S': 32.065,
    'Cl': 35.453,'Ar': 39.948,'K': 39.0983,'Ca': 40.078,'Sc': 44.95591,'Ti': 47.867,'V': 50.9415,'Cr': 51.9961,'Mn': 54.938049,
    'Fe': 55.845,'Co': 58.9332,'Ni': 58.6934,'Cu': 63.546,'Zn': 65.409,'Ga': 69.723,'Ge': 72.64,'As': 74.9216,
    'Se': 78.96,'Br': 79.904,'Kr': 83.798,'Rb': 85.4678,'Sr': 87.62,'Y': 88.90585,'Zr': 91.224,'Nb': 92.90638,
    'Mo': 95.94,'Tc': 98,'Ru': 101.07,'Rh': 102.9055,'Pd': 106.42,'Ag': 107.8682,'Cd': 112.411,'In': 114.818,
    'Sn': 118.71,'Sb': 121.76,'Te': 127.6,'I': 126.90447,'Xe': 131.293,'Cs': 132.90545,'Ba': 137.327,'La': 138.9055,
    'Ce': 140.116,'Pr': 140.90765,'Nd': 144.24,'Pm': 145,'Sm': 150.36,'Eu': 151.964,'Gd': 157.25,'Tb': 158.92534,
    'Dy': 162.5,'Ho': 164.93032,'Er': 167.259,'Tm': 168.93421,'Yb': 173.04,'Lu': 174.967,'Hf': 178.49,'Ta': 180.9479,
    'W': 183.84,'Re': 186.207,'Os': 190.23,'Ir': 192.217,'Pt': 195.078,'Au': 196.96655,'Hg': 200.59,'Tl': 204.3833,
    'Pb': 207.2,'Bi': 208.98038,'Po': 209,'At': 210,'Rn': 222,'Fr': 223,'Ra': 226,'Ac': 227,
    'Th': 232.0381,'Pa': 231.03588,'U': 238.02891,'Np': 237,'Pu': 244,'Am': 243,'Cm': 247,'Bk': 247,
    'Cf': 251,'Es': 252,'Fm': 257,'Md': 258,'No': 259,'Lr': 262,'Rf': 261,'Db': 262,
    'Sg': 266,'Bh': 264,'Hs': 277,'Mt': 268,'Ds': 281,'Rg': 272,'Cn': 285,'Nh': 286,
    'Fl': 289,'Mc': 289,'Lv': 293,'Ts': 294,'Og': 294,
    }

    freq_wd = [
    5,15,25,35,45,55,65,75,85,95,105,115,125,135,145,155,165,175,185,195,
    205,215,225,235,245,255,265,275,285,295,305,315,325,335,345,355,365,375,385,395,
    405,415,425,435,445,455,465,475,485,495,505,515,525,535,545,555,565,575,585,595,
    605,615,625,635,645,655,665,675,685,695,705,715,725,735,745,755,765,775,785,795,
    805,815,825,835,845,855,865,875,885,895,905,915,925,935,945,955,965,975,985,995,
    1005,1015,1025,1035,1045,1055,1065,1075,1085,1095,1105,1115,1125,1135,1145,1155,1165,1175,1185,1195,
    1205,1215,1225,1235,1245,1255,1265,1275,1285,1295,1305,1315,1325,1335,1345,1355,1365,1375,1385,1395,
    1405,1415,1425,1435,1445,1455,1465,1475,1485,1495,1505,1515,1525,1535,1545,1555,1565,1575,1585,1595,
    1605,1615,1625,1635,1645,1655,1665,1675,1685,1695,1705,1715,1725,1735,1745,1755,1765,1775,1785,1795,
    1805,1815,1825,1835,1845,1855,1865,1875,1885,1895,1905,1915,1925,1935,1945,1955,1965,1975,1985,1995,
    2005,2015,2025,2035,2045,2055,2065,2075,2085,2095,2105,2115,2125,2135,2145,2155,2165,2175,2185,2195,
    2205,2215,2225,2235,2245,2255,2265,2275,2285,2295,2305,2315,2325,2335,2345,2355,2365,2375,2385,2395,
    2405,2415,2425,2435,2445,2455,2465,2475,2485,2495]

    rf_pfw = [
    0.001,0.0056,0.0145,0.0213,0.0363,0.0421,0.0676,0.0724,0.0978,0.112,0.131,0.191,0.174,0.202,0.255,0.233,0.348,0.287,0.388,0.439,
    0.337,0.484,0.42,0.584,0.514,0.536,0.665,0.611,0.62,0.776,0.678,0.825,0.732,0.856,0.862,0.836,1.02,0.99,1.06,1,
    1.28,1.08,1.33,1.41,1.34,1.26,1.65,1.46,1.68,1.88,1.67,1.67,1.97,2.2,1.96,2.21,2.01,2.13,1.61,1.19,
    1.36,0.54,0.864,0.365,0.0777,0.0573,0.0003,0.0496,0.0777,0.222,0.51,0.806,0.607,1.32,1.44,2.37,2.95,3.03,3.15,2.74,
    3.09,3.2,3.19,3.23,3.16,3.05,3.17,3.1,3.04,3.1,3.01,3.03,2.93,2.95,2.81,2.81,2.81,2.67,2.72,2.61,
    2.31,1.82,1.47,1.37,1.51,1.11,1.78,2.18,2.22,2.18,2.1,2.07,1.99,1.79,1.88,1.79,1.72,1.52,1.52,1.52,
    1.67,1.23,1.15,1.26,0.906,0.815,0.418,0.455,0.411,0.359,0.207,0.296,0.371,0.237,0.231,0.28,0.173,0.18,0.181,0.122,
    0.182,0.098,0.12,0.0907,0.122,0.0569,0.0917,0.0613,0.0658,0.0528,0.0372,0.0462,0.0341,0.0335,0.0325,0.0281,0.0337,0.042,0.0739,0.0635,0.0518,
    0.0401,0.0346,0.0274,0.0261,0.019,0.0242,0.0233,0.0194,0.0167,0.0182,0.0179,0.0353,0.0156,0.0226,0.0265,0.0272,0.02,0.0308,0.0199,0.0402,
    0.04,0.03,0.0315,0.0228,0.0569,0.0255,0.066,0.0439,0.0561,0.0452,0.0404,0.0496,0.0735,0.036,0.077,0.0581,0.106,0.081,0.0466,0.131,
    0.0708,0.115,0.133,0.0767,0.106,0.0582,0.0543,0.0953,0.0752,0.101,0.0883,0.0912,0.086,0.0985,0.098,0.0863,0.0789,0.0494,0.0309,0.0202,
    0.0176,0.0204,0.0078,0.0053,0.003,0.0016,0.0007,0.0006,0.0001,0.0001,0,0,0,0,0,0,0,0.0027,0.0451,0.0473,0.0455,0.0436,0.0412,0.0385,0.0373,0.0367,0.0338,0.0339,0.0333]

    calc_rf = []
    for i in np.arange(0,len(freq_wd)):
        calc_rf.append([0])

    # Calculate radiative transfer for species i
    # -------------------------------------------

    freqs_considered = []
    ints_considered = []
    for freq,ints in zip(mol_data['Freq_Scaled'],mol_data['Intensity']):
        if freq <= 2500:
            freqs_considered.append(freq)
            ints_considered.append(ints)

    for location in np.arange(0,len(freq_wd)):
        for freq,ints in zip(freqs_considered,ints_considered):
            if freq_wd[location]-5 < freq < 5+freq_wd[location]:
                calc_rf[location].append((rf_pfw[location]/(10**-15))*ints*(1/(6.023*10**23))*1000*100)

    ai = 0
    for i in calc_rf:
        value = sum(i)
        ai = ai + value
    
    # Calculate molecular mass for species i
    # ---------------------------------------
    mass = 0
    for elements in re.findall('[A-Z][^A-Z]*', molecule.split('_')[0].rstrip()):
        elements_list = re.split('(\d+)', elements.split()[0])

        if len(elements_list) == 1:
            mass = mass + periodic_table[elements_list[0]]

        else:
            elements_list.pop(-1)
            mass = mass + float(elements_list[-1])*(periodic_table[elements_list[0]])
        
    # Calculate lifetime for species i
    # ---------------------------------

    k = A*np.exp(-ER_factor*(1/277))
    lifetime = 5.7*(6.62*10**-15/k)

    # Calculate GWP for species i
    # -----------------------------
    AGWP_CO2 = {20:0.192, 50:0.676, 100:2.223}

    func = lambda t: np.exp(-t/lifetime)
    GWP = ai*(integrate.quad(func,0,th)[0]/AGWP_CO2[th])*(44.0095/mass)*1000
    
    return print('Computed radiative forcing (ai) for',molecule,':',ai,'W/(m^2*ppb)','\nGWP',GWP,'kgCO2/kg'+molecule+'\n\nCalculated k: ',k,'(cm^3/(mol*sec))\nCalculated lifetime: ',lifetime,' years')
    
    # --------------------------------------------------------------------------------
