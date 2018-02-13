get_impute_data = function(index, data,two=T,segs) {
  if(two) {
    if(index == 1) {return(data[2:3])}
    if(index == segs) {return(data[(index-1):(index-2)])}
    return(c(data[index-1] ,data[index+1]))
  } else {
    if(index == 1) {return(data[2:5])}
    if(index == 2) {return(c(data[1], data[3:5]))}
    if(index == (segs-1)) {return(c(data[index+1], data[(index-1):(index-3)]))}
    if(index == segs) {return(data[(index-1):(index-4)])}
    return(c(data[index-1], data[index+1], data[index-2], data[index+2]))
  }
}
