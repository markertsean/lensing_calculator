#ifndef INPUT_FUNCTIONS
#define INPUT_FUNCTIONS

#include <cosmo.h>
#include <iostream>

void ReadNbodyHalo( 
		   PosType           **xpos,  //2D position of particles we will be filling in   
		   IndexType        Npoints,  //Number of particles/lines to read
		   std::string  inpFileName); //Name of the file to read


void   ReadInpFile(
		   userInfo        &inpInfo,  //Object we write to, contains parameters governing options
		   std::string userFileName); //Name of the file to read


void setCosmoParameters( 
			InputParams    params,  //Parameters used by GLAMER code
			COSMOLOGY   &inpCosmo); //Cosmology, default defined by user



void  setHaloParameters( 
			InputParams    params, //Parameters used by GLAMER code    
			COSMOLOGY   &inpCosmo, //Cosmology, default defined by user
			haloInfo     &inpInfo, //Object we are filling from input
 
			std::string   objType="lens"); //Default lens, can be source



#endif
