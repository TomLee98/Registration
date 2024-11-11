此文件夹将每个文件的 3D 原始体积数据存储为 *.tif。

您可以按照以下流程生成训练体积：
[1] 通过 Fiji 读取显微镜原始图像
[2] 通过“Image”->“Duplicate”或“Ctrl+Shift+D”选择一个体积（如果维度大于 3）
[3] 选择 ROI，然后通过“Image”->“Crop”或“Ctrl+Shift+X”裁剪
[4] 通过“File”->“Save As”->“Tiff...”将体积保存为当前文件夹中的 tiff 文件


This folder stores the 3D raw volumes data each file as *.tif. 

You can follow the next pipeline to generate your training volume:
[1] Read your microscope raw image by Fiji
[2] select one volume (if dimension more than 3) by 'Image'->'Duplicate' or 'Ctrl+Shift+D'
[3] select ROI, then crop by 'Image'->'Crop' or 'Ctrl+Shift+X'
[4] save the volume as tiff file by 'File'->'Save As'->'Tiff...' in current folder