Before using, we need to setup the parallel computing environment:
0. Make sure your MATLAB Version is not lower than R2021b.
1. If you are using this tool on your local machine, please use parallel 
profile 'local' as default. If your MATLAB version is higher than R2022a,
use profile 'processes' as recommanded.
2. If you are using this tool on computing server(IW4220@SiLab), please 
set the parallel profile 'Reg3D_Server.mlsettings' in 'bin\' for stable running.
NOTE: For advanced user, you can modify the parallel setting in manual if 
there are some parallel bugs (lost connection, out of memory, ... etc).
If bugs come, just reduce the number of workers and run registration again.
NOTE: For advanced user, you can modify the threads number per worker, such
as setting between 1 and 2. Because of the 'hyper-thread' technology, you can 
get fully performance, but parpool may be not very stable. (modify the file
register3D_tc_auto.m -> line 1035, register3D_oc_auto.m -> line 
or Reg3D_Server.mlsettings)

[Two Channels Registration Pipeline]
1. select the template as fix frame, for better performance, we recommand
you select the local stable best
2. using 'translate' as transformation for object location estimating
3. set the initial iteration step bigger than default(6.25e-2), which could 
be choosed from 0.1 to 0.15
4. set the iteration factor bigger than default(0.5), which could be choosed
from 0.6 to 0.8
5. set the robustness factor bigger than default(0.5), which could be choosed
from 0.6 to 0.9
6. check the 'Enable' checkbox to confirm registration(global) setting
7. make sure you are working on source 'raw', and click the 'Run' button
Until the app is stopped, you can go to next step#

8. if 'translate' is not accurate enough, and the object has some rotation 
motion or stretch motion, using 'rotation' or 'similarity' as transformation
9. you may decrease the minimal iteration step for good accuracy, such as 
1e-6 or smaller, where 1e-5 as default
10. set the initial iteration step as default(6.25e-2)
11. set the iteration factor as default(0.5) or smaller
12. set the robustness factor smaller than default(0.5), which could be choosed
from 0.01 to 0.2
13. you can crop the image in XY plane for better registration speed 
14. make sure you are working on source 'registrated', and click the 'Run' button
Until the app is stopped, you can go to next step#
[Warning] Do not run on 'registrated' too many times(more than 3 times), 
because of interpolation may decrease the spatial and temporal resolution.

15. if there is also some local motion, you can choose registration(local)
setting panel, using the template which came from 'global' is better
16. modify the cumulative field smooth value, the value is smaller, the 
local motion registration is stronger
17. you can skip this step, if marker channel with low photobleaching, or 
you can check the 'compensate photobleaching' checkbox for better performance
18. check the 'Enable' checkbox to confirm registration(local) setting
19. make sure you are working on source 'registrated', and click the 'Run' button
Until the app is stopped, you can go to next step#

20. if there are some bad registration frames, you can cut them or use 
registration(manual) for manual resiatration

[One Channel Registration Pipeline]
1. 

Go and have fun!


Copy Rights
Weihan Li, 2022-2024