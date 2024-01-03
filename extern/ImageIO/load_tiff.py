from tifffile import imread

def loadtiff(file_name):
    try:
        vol = imread(file_name)
    except IOError:
        status = -1
    else:
        status = 0
    
    return status, vol

~, vol = loadtiff(file)