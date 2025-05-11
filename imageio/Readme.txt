You can use python package tifffile to speed up tiff file loading and saving.

Configuration Pipeline:
[For Windows OS User]:
(1) Install Anaconda (with adaptable python version to MATLAB)
(2) Make sure '<install folder>\anaconda3\' and '<install folder>\anaconda3\Scripts\' in the environment variables
(3) Create new environment named 'reg3d' or other with tifffile package under Anaconda
(4) Reg3D App -> Advanced -> Python Interpreter, replace python.exe with the instance in 'reg3d'

[For Linux/Unix OS User]:
Try your best! You can do it.

================================================================================
You can use Nikon SDK for fast reading *.nd2 file. Note that the support is enabled as default.
If you open .nd2 file with MATLAB warning in computing server(Ubuntu)@silab: 
'Warning: There was an error loading the library "/home/<youe account>/MATLAB
Add-Ons/Apps/Reg3D<Version>/imgio/ImageReader/NikonReader/bin/linux/libNd2ReadSdk.so"
<library name>.so.52: cannot open shared object file: No such file or directory'
Try to connect the server administrator and run the following commands:
(1) put the following files under '/usr/lib/', > sudo cp ./<library name>.so.52 /usr/lib/
(2) set the mod as '755', > sudo chmod 755 /usr/lib/<library name>.so.52
(3) refresh dynamic library buffer, > sudo ldconfig
(4) check the libraries status, files will be located in memory, > ldd ./<library name>.so.52
(5) restart your MATLAB and Reg3D App
[library]
libicuuc.so.52
libicudata.so.52
libicui18n.so.52
libicutu.so.52

Note that these files will be found in folder '/home/<youe account>/MATLAB
Add-Ons/Apps/Reg3D<Version>/imgio/ImageReader/NikonReader/extern/linux/'
