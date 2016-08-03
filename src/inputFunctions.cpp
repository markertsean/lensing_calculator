#include <cstring>
#include <CCfits/CCfits>
#include <slsimlib.h>
#include <lensing_classes.h>




// Read the FITS header of our file
void readFitsHeader( const std::string inputFile ,  // Name of the FITS file
                           haloInfo      &myHalo ,  // Halo to put info into
                           userInfo   &userInput ){ // User info to take from header


  // Opens fits file to read header
  std::auto_ptr<CCfits::FITS> ff(new CCfits::FITS(inputFile, CCfits::Read));
	CCfits::PHDU* h0 = &ff->pHDU();

  // Values to read in
  std::string catalog;
  int N_pixels_v, N_pixels_h, id;
  float fov, redshift, mass, r_vir, c, ba, ca, phi, theta, integ, physicalsize;

  // Attempt to read header items, if fail abort
  try {
    h0->readKey( "CATALOG"      , catalog      );
    h0->readKey( "FOV"          , fov          );
    h0->readKey( "N_pixels_v"   , N_pixels_h   );
    h0->readKey( "N_pixels_h"   , N_pixels_v   );
    h0->readKey( "PHYSICALSIZE" , physicalsize );
    h0->readKey( "REDSHIFT"     , redshift     );
    h0->readKey( "ID"           , id           );
    h0->readKey( "MASS"         , mass         );
    h0->readKey( "RVIR"         , r_vir        );
    h0->readKey( "C"            , c            );
    h0->readKey( "b/a"          , ba           );
    h0->readKey( "c/a"          , ca           );
    h0->readKey( "PHI"          , phi          );
    h0->readKey( "THETA"        , theta        );
    h0->readKey( "INTEG"        , integ        );
  }
  catch (...) {

    std::string aborting = std::string("FITS header must contain: " )+
                           std::string("CATALOG, "      )+
                           std::string("FOV, "          )+
                           std::string("N_pixels_v, "   )+
                           std::string("N_pixels_h, "   )+
                           std::string("PHYSICALSIZE, " )+
                           std::string("REDSHIFT, "     )+
                           std::string("ID, "           )+
                           std::string("MASS, "         )+
                           std::string("RVIR, "         )+
                           std::string("C, "            )+
                           std::string("b/a, "          )+
                           std::string("c/a, "          )+
                           std::string("PHI, "          )+
                           std::string("THETA, "        )+
                           std::string("INTEG, "        );

    std::cout << aborting << std::endl << "Aborting." << std::endl;

    logMessage(  aborting   );
    logMessage( "Aborting." );
    exit(1);
  }

  logMessage(  std::string("FITS header contains: " )+
               std::string("CATALOG, "      )+                              catalog        + std::string("\n") +
               std::string("FOV, "          )+ std::to_string((long double) fov          ) + std::string("\n") +
               std::string("N_pixels_v, "   )+ std::to_string((long long  ) N_pixels_v   ) + std::string("\n") +
               std::string("N_pixels_h, "   )+ std::to_string((long long  ) N_pixels_h   ) + std::string("\n") +
               std::string("PHYSICALSIZE, " )+ std::to_string((long double) physicalsize ) + std::string("\n") +
               std::string("REDSHIFT, "     )+ std::to_string((long double) redshift     ) + std::string("\n") +
               std::string("ID, "           )+ std::to_string((long long  ) id           ) + std::string("\n") +
               std::string("MASS, "         )+ std::to_string((long double) mass         ) + std::string("\n") +
               std::string("RVIR, "         )+ std::to_string((long double) r_vir        ) + std::string("\n") +
               std::string("C, "            )+ std::to_string((long double) c            ) + std::string("\n") +
               std::string("b/a, "          )+ std::to_string((long double) ba           ) + std::string("\n") +
               std::string("c/a, "          )+ std::to_string((long double) ca           ) + std::string("\n") +
               std::string("PHI, "          )+ std::to_string((long double) phi          ) + std::string("\n") +
               std::string("THETA, "        )+ std::to_string((long double) theta        ) + std::string("\n") +
               std::string("INTEG, "        )+ std::to_string((long double) integ        ) );


  // Put halo values into halo data type
  myHalo.setID    (       id );
  myHalo.setM     (     mass );
  myHalo.setZ     ( redshift );
  myHalo.setRmax  (    r_vir / 1000);
  myHalo.setC     (        c );
  myHalo.setBA    (       ba );
  myHalo.setCA    (       ca );
  myHalo.setPhi   (      phi );
  myHalo.setTheta (    theta );

  // Put values of image into the user input data type
  userInput.setCatType    ( catalog      );
  userInput.setIntegLength( integ        );
  userInput.setPhysFOV    ( fov          );
  userInput.setAngFOV     ( physicalsize );
  userInput.setNpixV      ( N_pixels_v   );
  userInput.setNpixH      ( N_pixels_h   );
  userInput.setNpix       ( N_pixels_h   *
                            N_pixels_v   );

}



// Read parameters not included in paramfile (Nbins, Nthreads, etc.)
void readInpFile(          userInfo  &inpInfo  ,   // Info needed for the rest of the code
                  const std::string  inputFile ){  // Input file name
  FILE *pFile;
  char   inpC1[35],inpC2[35];

  // Attempt to open file, if successful go line by line
  //  finding the variable name and value

  pFile = fopen(inputFile.c_str(),"r");

  if (pFile!=NULL){

    logMessage( std::string( "Reading file: ") + inputFile );

    // Scan variables
    while ( fscanf(pFile,"%s%s",inpC1,inpC2) != EOF ){
      std::string inpS = std::string(inpC1);
           if ( inpS=="N_bins"      ){        inpInfo.setNbins    (        atoi(inpC2) );      }
      else if ( inpS=="N_sources"   ){        inpInfo.setNsrc     (        atoi(inpC2) );      }
      else if ( inpS=="N_threads"   ){        inpInfo.setNthreads (        atoi(inpC2) );      }
      else if ( inpS=="sourceRadius"){        inpInfo.setSourceRadius(     atof(inpC2) );      }
      else if ( inpS=="cosmo"       ){        inpInfo.setCosmology( std::string(inpC2) );      }
      else if ( inpS=="fox2012F"    ){        inpInfo.setFoxH2012F( std::string(inpC2) );      }
      else if ( inpS=="fox2123F"    ){        inpInfo.setFoxH2123F( std::string(inpC2) );      }
      else{

          // Abort if unrecognized variables

          std::cout << " Couldn't recognize input from " << inputFile <<
                     ": " << inpS << std::endl << std::endl;
          logMessage( std::string("Unrecognized input: ") + inpS );
        exit(1);
      }
    }
    fclose(pFile);
  }
  // Abort if couldn't open file
  else{
                std::cout << "Couldn't find file: " << inputFile << std::endl;
    logMessage( std::string( "Couldn't find file: ") + inputFile );
    logMessage( std::string( "Aborting." ) );
    exit(1);
  }

  // Required parameters, abort if missing
  if ( inpInfo.getCosmology() == " " ||
       inpInfo.getNbins    () == -1  ||
       inpInfo.getNsrc     () == -1  ){

    std::cout << inputFile <<             " must contain cosmo     = WMAP or PLANCK" << std::endl <<
                                          "              N_bins    = #"              << std::endl <<
                                          "              N_sources = #"              << std::endl ;

    logMessage(  inputFile + std::string( " must contain cosmo=WMAP or PLANCK" ) );
    logMessage(  inputFile + std::string( " must contain N_bins=#" ) );
    logMessage(  inputFile + std::string( " must contain N_sources=#" ) );
    logMessage(              std::string( "Aborting." ) );
    exit(1);
  }

  logMessage( std::string("N_bins    = ") + std::to_string((long long) inpInfo.getNbins()     ) +
              std::string("N_sources = ") + std::to_string((long long) inpInfo.getNsrc()      ) +
              std::string("N_threads = ") + std::to_string((long long) inpInfo.getNthreads()  ) +
              std::string("cosmo     = ") + std::string   (            inpInfo.getCosmology() ) +
              std::string("FoxH2012F = ") + std::string   (            inpInfo.getFoxH2012F() ) +
              std::string("FoxH2123F = ") + std::string   (            inpInfo.getFoxH2123F() ) );

}



// Reads the fox H tables, saved in log
einTable readFoxH( userInfo &u, const int fileType ){

  double minX, maxX;
  double minA, maxA;

  int x_bins, a_bins;

  einTable einKappa;

//  std::string          myFile = u.getFoxH2012F();//"src/foxH2012.dat";
//  if ( fileType == 2 ) myFile = u.getFoxH2123F();//"src/foxH2123.dat";
  std::string          myFile = "src/foxH2012.dat";
  if ( fileType == 2 ) myFile = "src/foxH2123.dat";

  FILE *pFile;

  pFile = fopen(myFile.c_str(),"r");

  if (pFile!=NULL){

//    logMessage( std::string( "Reading file: ") + myFile );

    fscanf( pFile, "%16lf%16lf%4i",&minX,&maxX,&x_bins);
    fscanf( pFile, "%16lf%16lf%4i",&minA,&maxA,&a_bins);

    // Allocate the files
    einKappa.setX_min( minX );
    einKappa.setX_max( maxX );
    einKappa.setA_min( minA );
    einKappa.setA_max( maxA );

    einKappa   .setBins( a_bins, x_bins );


//    logMessage( std::string( "Number of alpha bins: ") + std::to_string( (long long) a_bins ) );
//    logMessage( std::string( "Number of x     bins: ") + std::to_string( (long long) x_bins ) );


    for ( int i = 0; i < a_bins; ++i ){ // Each row is a new alpha
    for ( int j = 0; j < x_bins; ++j ){ // Each column is a different x

      double inpVal;

      fscanf( pFile, "%16lf", &inpVal); // Goes across rows, then down columns

      einKappa.setVal( i, j, inpVal );

    }
    }

  } else {

//    logMessage( std::string( "Cannot open FoxH file: ") + myFile );

    std::cout << "Couldn't open FoxH file: " << myFile << std:: endl;
    exit(0);

  }

//  logMessage( std::string( "FoxH Read in complete" ) );

  return einKappa;
}
