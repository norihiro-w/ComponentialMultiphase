GeoSys-ST: Source Terms ------------------------------------------------
#SOURCE_TERM
 $PCS_TYPE
  CMP_CapPresFormular
 $PRIMARY_VARIABLE
  NONWET_PRESSURE
 $GEO_TYPE
  POLYLINE PLY_LEFT
 $DIS_TYPE   
  CONSTANT_NEUMANN 5.57e-6 
#STOP


 $DIS_TYPE  
  CONSTANT_NEUMANN 0.0002

 $DIS_TYPE  
  RIVER 2
  4  3 1.00e-6 10.0 1.3 1.0
  5  3 1.00e-6 10.0 1.3 1.0