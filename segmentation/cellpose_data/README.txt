您可以按照以下流程生成自定义的 3D 数据训练集：

[1] 生成原始体积数据，详情请参阅 ./raw_files/README.txt
[2] 使用脚本“raw2train.m”生成训练数据和测试数据
[3] 使用“labelset.m”手动生成训练数据和测试数据的掩码

如果您至少有一张高性能显卡（>=Nvidia RTX2070，显存>=8GB），
您可以在 ReTiNA 应用中基于“cyto”训练自定义模型。

请注意，目前仅支持“ORN”和“KC”细胞/核训练。

玩得开心！:)



You can follow the next pipeline to generate custom training set of 3D data:

[1] generate raw volume data, details in ./raw_files/README.txt
[2] generate the train data and test data by using script 'raw2train.m'
[3] generate mask of train data and test data in manual by using 'labelset.m'

If you have one high performance graphics card (>=Nvidia RTX2070, GPU Memory>=8GB) at least, 
you can train a custom model based on 'cyto' in ReTiNA app.
Note that only "ORN" and "KC" cell/nuclear training are supported currently.

Have fun! :)
