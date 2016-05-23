#ifndef INPUT_FUNCTIONS
#define INPUT_FUNCTIONS

#include <cosmo.h>
#include <iostream>


void readFitsHeader( const std::string inputFile ,  // Name of the FITS file
                           haloInfo      &myHalo ,  // Halo to put into into
                           userInfo   &userInput ); // User info to take from header


/*
void   ReadInpFile( 	      userInfo        &inpInfo ,  // Object we write to, contains parameters governing options
                   const std::string    userFileName ); // Name of the file to read
*/


/*
void ReadNbodyHalo(
		   double         xpos[][3],  //2D position of particles we will be filling in
		   int              Npoints,  //Number of particles/lines to read
		   std::string  inpFileName); //Name of the file to read
*/
/*
void setCosmoParameters(
			InputParams    params,  //Parameters used by GLAMER code
			COSMOLOGY   &inpCosmo); //Cosmology, default defined by user



void  setHaloParameters(
			InputParams    params, //Parameters used by GLAMER code
			COSMOLOGY   &inpCosmo, //Cosmology, default defined by user
			haloInfo     &inpInfo, //Object we are filling from input

			std::string   objType="lens"); //Default lens, can be source
*/


#endif
