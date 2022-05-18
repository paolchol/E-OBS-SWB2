update_targpest = function(file, outpath, ext = '.dat'){
  
  #ext: extension
  
  outfile = paste0(outpath, 'targpest_updated', ext)
  
  f = read.fwf(file)
  write(f, file = '')
  
}

