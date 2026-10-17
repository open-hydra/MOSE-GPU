from ORION import read_TEC, write_TEC

[x,y,z,var,names] = read_TEC("field.tec")

var[0][0] *= 10

write_TEC('out.tec',x,y,z,var)