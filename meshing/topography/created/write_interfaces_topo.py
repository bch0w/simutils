"""
Copy-pasted from prepare_meshfem.py and removed from there too. 
Hard coding topography lines in the interfaces.dat file is not required if this
script is to be mroe general so I am just moving that here for now
"""
def write_interfaces(template, dir_name, layers, interfaces, lat_min, lon_min,   
                     suppress_utm_proj, fids=[], topo="default"):                
    """                                                                          
    Write the interfaces.dat file as well as the corresponding flat interface    
    layers. Topo will need to be written manually                                
                                                                                 
    :type template: str                                                          
    :param template: fid of the interfaces.dat template, preformatted            
    :type dir_name: str                                                          
    :param dir_name: directory to store output files                             
    :type layers: list                                                           
    :param layers: number of elements in each layer                              
    :type interfaces: list                                                       
    :param interfaces: depths in km of each interface                            
    :param lat_min: minimum latitude of the mesh                                 
    :param lon_min: minimum longitude of the mesh                                
    :type suppress_utm_proj: bool                                                
    :param suppress_utm_proj: if True, will use Lat/Lon coordinates, else        
        uses a UTM projection conversion                                         
    :type fids: list                                                             
    :param fids: for custom interface names, starting from the top going down,   
        excluding topography and the bottom of the mesh. If left blank, defaults 
        to 1, 2 ... until the number of layers                                   
    """                                                                          
    # Ensure counting layers from bottom                                         
    layers_from_bottom = layers   

    # Template for setting a flat layer                                          
    flat_layer = (f".{str(suppress_utm_proj).lower()}. 2 2 {lon_min:.1f}d0 "     
                  f"{lat_min:.1f}d0 180.d0 180.d0")                              
                                                                                 
    # Hardcoded topo layers define the structure of the underlying               
    # single-column topography file that must be generated externally            
    topo_fid = "interface_topo.dat"                                              
    logger.info(f"\tsetting topography to '{topo}', points to '{topo_fid}'")     
    if topo == "nznorth":                                                        
        topo = ".false. 720 720 173.d0 -43.d0 0.00833d0 0.00833d0"               
    elif topo == "nzsouth":                                                      
        topo = ".true. 899 859 38192d0 -5288202d0 1000.00d0 1000.00d0"           
    elif topo == "nznorth_ext":                                                  
        topo = ".true. 763 850 115822.d0 5358185.d0 1000.00d0 1000.00d0"         
    elif topo == "c2s_nznorth_ext":                                              
        topo = ".true. 560 818 -361383.5d0 -405861.5d0 1290.7d0 992.3d0"         
    elif topo == "c2s_nalaska":                                                  
        topo = ".true. 288 293 104449d0 6903153d0 5000.00d0 5000.00d0"           
    else:                                                                        
        topo = flat_layer                                                        
                                                                                 
    # Write to a new file                                                        
    with open(os.path.join(dir_name, "interfaces.dat"), "w") as f:               
        f.write("# number of interfaces\n")                                      
        f.write(" {ninterfaces}\n".format(ninterfaces=len(interfaces) + 1))      
                                                                                 
        # Write flat layers up to topo                                           
        for i, (interface, fid) in enumerate(zip(interfaces, reversed(fids))):   
            f.write("# interface number {}\n".format(i+1))                       
            f.write(f" {flat_layer}\n")                                          
            f.write(f" {fid}\n")                                                 
                                                                                 
        # Write topo interface                                                   
        f.write("# interface number {} (topo)\n".format(i+2))                    
        f.write(f" {topo}\n")                                                    
        f.write(f" {topo_fid}\n")                                                
                                                                                 
        # Write number of elements in each layer                                 
        f.write("# number of spectral elements per layer (from bottom)\n")       
        for j, layer in enumerate(layers_from_bottom):                           
            f.write(f" {layer}\n")                                               
                                                                                 
    # Write the individual interface files                                       
    for fid, interface in zip(fids, interfaces):                                 
        # Skip top interface (topography)                                        
        logger.info(f"WRITING INTERACE {fid}")                                   
        with open(os.path.join(dir_name, fid), "w") as f:                        
            for i in range(4):                                                   
                f.write("-{}\n".format(abs(int(interface * 1E3))))   
