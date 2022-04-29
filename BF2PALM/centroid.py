def centroid(vertexes):
    x_list = [vertex [0] for vertex in vertexes]
    y_list = [vertex [1] for vertex in vertexes]
    len = len(vertexes)
    x = sum(x_list)/_len
    y = sum(y_list)/_len
    
    return(x, y)