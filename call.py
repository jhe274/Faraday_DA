import os, Plot

cwd = os.getcwd()
dir_path = os.path.join(cwd, 'Research', 'PhD Project', 'Faraday Rotation Measurements', 'K vapor cell')
Bristol = os.path.join(dir_path, 'Bristol data')
TC300 = os.path.join(dir_path, 'TC300 data')
TriMod = os.path.join(dir_path, 'Triple modulation data')
Analysis = os.path.join(dir_path, 'Data Analysis')

Plot.theta_mea_vs_nu(Bristol, TriMod, 12, 5)