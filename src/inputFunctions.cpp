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
  myHalo.setRmax  (    r_vir );
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
              std::string("cosmo     = ") + std::string   (            inpInfo.getCosmology() ) );

}

/*Reads Nbody particle input, only takes position "f10f10f10"
//
//input:
//      PosType 2D object, accepts the xyz coordinates from dark matter particles
//
//output:
//      void
//      Positions of halos in xpos, PosType 2D object [xyz][index]
//
*//*
void ReadNbodyHalo( double xpos[][3], int Npoints, std::string inpFileName){

  int  iNlines(Npoints);

  FILE * pFile;
  double x,y,z;
  pFile = fopen(inpFileName.c_str(),"r");

  for (int ii=0;ii<iNlines;ii++){
    fscanf(pFile,"%10lf%10lf%10lf",&x,&y,&z);
    xpos[ii][0]=x;
    xpos[ii][1]=y;
    xpos[ii][2]=z;
  }
  fclose(pFile);
}*/


/*gets cosmo parameters from paramfile, puts in cosmology object
//
//input:
//     param object, cosmology object to write to
//
//output:
//     void, inpCosmo to write to
//
*/
/*
void setCosmoParameters( InputParams params, COSMOLOGY &inpCosmo ){
  double tempVal;
  std::string tempStr;
  bool tempBol;

  //Need to specify the cosmology in input file
  if (params.get("Omega_baryon",tempVal))
    inpCosmo.setOmega_baryon(   tempVal);
  else
    std::cerr << " No param in paramfile: Omega_baryon\n" ;
  if (params.get("Omega_matter",tempVal))
    inpCosmo.setOmega_matter(   tempVal);
  else
    std::cerr << " No param in paramfile: Omega_matter\n" ;
  if (params.get("Omega_lambda",tempVal))
    inpCosmo.setOmega_lambda(   tempVal);
  else
    std::cerr << " No param in paramfile: Omega_lambda\n" ;
  if (params.get("Omega_neutrino",tempVal))
    inpCosmo.setOmega_neutrino(   tempVal);
  else
    std::cerr << " No param in paramfile: Omega_neutrino\n" ;
  if (params.get("hubble",tempVal))
    inpCosmo.sethubble(   tempVal);
  else
    std::cerr << " No param in paramfile: h\n" ;
  if (params.get("sigma_8",tempVal))
    inpCosmo.setSigma8(    tempVal);
  else
    std::cerr << " No param in paramfile: sigma_8\n" ;

  //By attempting to access each variable from paramfile, supresses
  // VERY annoying output of all variables used and unused
  params.get("MOKA_analyze"                     ,tempVal);
  params.get("MOKA_background_field"            ,tempVal);
  params.get("MOKA_input_file"                  ,tempStr);
  params.get("MOKA_input_params"                ,tempVal);
  params.get("Omega_baryon"                     ,tempVal);
  params.get("Omega_lambda"                     ,tempVal);
  params.get("Omega_matter"                     ,tempVal);
  params.get("Omega_neutrino"                   ,tempVal);
  params.get("QSO_colors_file"                  ,tempStr);
  params.get("QSO_kcorrection_file"             ,tempStr);
  params.get("SourceSBType"                     ,tempStr);
  params.get("deflection_off"                   ,tempVal);
  params.get("gamma_uniform_1"                  ,tempStr);
  params.get("gamma_uniform_2"                  ,tempStr);
  params.get("gauss_r2"                         ,tempVal);
  params.get("hubble"                           ,tempVal);
  params.get("kappa_uniform"                    ,tempStr);
  params.get("lensing_off"                      ,tempVal);
  params.get("main_DM_halo_type"                ,tempStr);
  params.get("main_NDistortionModes"            ,tempVal);
  params.get("main_Rmax"                        ,tempVal);
  params.get("main_axis_ratio"                  ,tempVal);
  params.get("main_concentration"               ,tempVal);
  params.get("main_core"                        ,tempVal);
  params.get("main_ellip_method"                ,tempStr);
  params.get("main_galaxy_halo_type"            ,tempVal);
  params.get("main_halo_on"                     ,tempVal);
  params.get("main_mass"                        ,tempVal);
  params.get("main_perturb_beta"                ,tempVal);
  params.get("main_perturb_gamma"               ,tempVal);
  params.get("main_perturb_hexopole"            ,tempVal);
  params.get("main_perturb_kappa"               ,tempVal);
  params.get("main_perturb_monopole"            ,tempVal);
  params.get("main_perturb_octopole"            ,tempVal);
  params.get("main_perturb_quadrapole"          ,tempVal);
  params.get("main_pos_angle"                   ,tempVal);
  params.get("main_rscale"                      ,tempVal);
  params.get("main_sigma"                       ,tempVal);
  params.get("main_slope"                       ,tempVal);
  params.get("main_stars_N"                     ,tempVal);
  params.get("main_stars_bending_point"         ,tempVal);
  params.get("main_stars_fraction"              ,tempVal);
  params.get("main_stars_hi_mass_slope"         ,tempVal);
  params.get("main_stars_lo_mass_slope"         ,tempVal);
  params.get("main_stars_mass"                  ,tempVal);
  params.get("main_stars_mass_function"         ,tempVal);
  params.get("main_stars_max_mass"              ,tempVal);
  params.get("main_stars_min_mass"              ,tempVal);
  params.get("main_sub_Ndensity"                ,tempVal);
  params.get("main_sub_Rmax"                    ,tempVal);
  params.get("main_sub_alpha"                   ,tempVal);
  params.get("main_sub_beta"                    ,tempVal);
  params.get("main_sub_mass_max"                ,tempVal);
  params.get("main_sub_mass_min"                ,tempVal);
  params.get("main_sub_type"                    ,tempVal);
  params.get("main_zlens"                       ,tempVal);
  params.get("outputfile"                       ,tempStr);
  params.get("pixelmap_padding_factor"          ,tempVal);
  params.get("pixelmap_zeromean"                ,tempBol);
  params.get("pixelmaps_input_file"             ,tempStr);
  params.get("pixelmaps_on"                     ,tempVal);
  params.get("read_redshift_planes"             ,tempVal);
  params.get("redshift_planes_file"             ,tempStr);
  params.get("shapelets_band"                   ,tempStr);
  params.get("shapelets_folder"                 ,tempStr);
  params.get("sigma_8"                          ,tempVal);
  params.get("source_BHmass"                    ,tempVal);
  params.get("source_band"                      ,tempVal);
  params.get("source_fK"                        ,tempVal);
  params.get("source_gamma"                     ,tempVal);
  params.get("source_inclin"                    ,tempVal);
  params.get("source_input_galaxy_file"         ,tempStr);
  params.get("source_mag_limit"                 ,tempVal);
  params.get("source_nuo"                       ,tempVal);
  params.get("source_opening_ang"               ,tempVal);
  params.get("source_r_in"                      ,tempVal);
  params.get("source_r_out"                     ,tempVal);
  params.get("source_sb_limit"                  ,tempVal);
  params.get("z_lens"                           ,tempVal);
  params.get("z_source"                         ,tempVal);
  params.get("zsource_reference"                ,tempVal);
  params.get("field_Nplanes"                    ,tempVal);
  params.get("field_buffer"                     ,tempVal);
  params.get("field_fov"                        ,tempVal);
  params.get("field_input_simulation_center_DEC",tempVal);
  params.get("field_input_simulation_center_RA" ,tempVal);
  params.get("field_input_simulation_format"    ,tempStr);
  params.get("field_input_simulation_path"      ,tempStr);
  params.get("field_input_simulation_radius"    ,tempVal);
  params.get("field_internal_profile"           ,tempStr);
  params.get("field_internal_profile_galaxy"    ,tempVal);
  params.get("field_mass_func_alpha"            ,tempVal);
  params.get("field_mass_func_type"             ,tempVal);
  params.get("field_min_mass"                   ,tempVal);
  params.get("field_off"                        ,tempVal);
  params.get("field_prof_internal_slope_pl"     ,tempVal);
  params.get("field_prof_internal_slope_pnfw"   ,tempVal);

}*/

/*Read in from paramfile halo info, redshift mass etc
//
//input:
//     paramfile var, haloinfo var, string with either lens or source
//
//output:
//     void, inpInfo overwritten
//
*//*
void setHaloParameters(InputParams params, COSMOLOGY &inpCosmo, haloInfo &inpInfo, std::string objType){
  double tempVal;
  std::string tempStr;
  bool tempBol;
  inpInfo.setCosmo(inpCosmo);
  if ( objType=="lens"){
    if (params.get("field_internal_profile",tempStr))
      inpInfo.setProfile( tempStr );
    if (params.get("main_zlens"            ,tempVal))
      inpInfo.setZ(       tempVal );
    if (params.get("main_Rmax"             ,tempVal))
      inpInfo.setRmax(    tempVal );
    if (params.get("main_mass"             ,tempVal))
      inpInfo.setM(       tempVal );
    if (params.get("main_concentration"    ,tempVal))
      inpInfo.setC(       tempVal );
  }
  else if ( objType=="source"){
    if (params.get("z_source"              ,tempVal))
      inpInfo.setZ(       tempVal );
    if (params.get("main_zlens"            ,tempVal))
      inpInfo.setLensZ(   tempVal );
  }
  else
    std::cerr << "Couldn't recognize inp halo type\n";

 //By attempting to access each variable from paramfile, supresses
  // VERY annoying output of all variables used and unused
  params.get("MOKA_analyze"                     ,tempVal);
  params.get("MOKA_background_field"            ,tempVal);
  params.get("MOKA_input_file"                  ,tempStr);
  params.get("MOKA_input_params"                ,tempVal);
  params.get("Omega_baryon"                     ,tempVal);
  params.get("Omega_lambda"                     ,tempVal);
  params.get("Omega_matter"                     ,tempVal);
  params.get("Omega_neutrino"                   ,tempVal);
  params.get("QSO_colors_file"                  ,tempStr);
  params.get("QSO_kcorrection_file"             ,tempStr);
  params.get("SourceSBType"                     ,tempStr);
  params.get("deflection_off"                   ,tempVal);
  params.get("gamma_uniform_1"                  ,tempStr);
  params.get("gamma_uniform_2"                  ,tempStr);
  params.get("gauss_r2"                         ,tempVal);
  params.get("hubble"                           ,tempVal);
  params.get("kappa_uniform"                    ,tempStr);
  params.get("lensing_off"                      ,tempVal);
  params.get("main_DM_halo_type"                ,tempStr);
  params.get("main_NDistortionModes"            ,tempVal);
  params.get("main_Rmax"                        ,tempVal);
  params.get("main_axis_ratio"                  ,tempVal);
  params.get("main_concentration"               ,tempVal);
  params.get("main_core"                        ,tempVal);
  params.get("main_ellip_method"                ,tempStr);
  params.get("main_galaxy_halo_type"            ,tempVal);
  params.get("main_halo_on"                     ,tempVal);
  params.get("main_mass"                        ,tempVal);
  params.get("main_perturb_beta"                ,tempVal);
  params.get("main_perturb_gamma"               ,tempVal);
  params.get("main_perturb_hexopole"            ,tempVal);
  params.get("main_perturb_kappa"               ,tempVal);
  params.get("main_perturb_monopole"            ,tempVal);
  params.get("main_perturb_octopole"            ,tempVal);
  params.get("main_perturb_quadrapole"          ,tempVal);
  params.get("main_pos_angle"                   ,tempVal);
  params.get("main_rscale"                      ,tempVal);
  params.get("main_sigma"                       ,tempVal);
  params.get("main_slope"                       ,tempVal);
  params.get("main_stars_N"                     ,tempVal);
  params.get("main_stars_bending_point"         ,tempVal);
  params.get("main_stars_fraction"              ,tempVal);
  params.get("main_stars_hi_mass_slope"         ,tempVal);
  params.get("main_stars_lo_mass_slope"         ,tempVal);
  params.get("main_stars_mass"                  ,tempVal);
  params.get("main_stars_mass_function"         ,tempVal);
  params.get("main_stars_max_mass"              ,tempVal);
  params.get("main_stars_min_mass"              ,tempVal);
  params.get("main_sub_Ndensity"                ,tempVal);
  params.get("main_sub_Rmax"                    ,tempVal);
  params.get("main_sub_alpha"                   ,tempVal);
  params.get("main_sub_beta"                    ,tempVal);
  params.get("main_sub_mass_max"                ,tempVal);
  params.get("main_sub_mass_min"                ,tempVal);
  params.get("main_sub_type"                    ,tempVal);
  params.get("main_zlens"                       ,tempVal);
  params.get("outputfile"                       ,tempStr);
  params.get("pixelmap_padding_factor"          ,tempVal);
  params.get("pixelmap_zeromean"                ,tempBol);
  params.get("pixelmaps_input_file"             ,tempStr);
  params.get("pixelmaps_on"                     ,tempVal);
  params.get("read_redshift_planes"             ,tempVal);
  params.get("redshift_planes_file"             ,tempStr);
  params.get("shapelets_band"                   ,tempStr);
  params.get("shapelets_folder"                 ,tempStr);
  params.get("sigma_8"                          ,tempVal);
  params.get("source_BHmass"                    ,tempVal);
  params.get("source_band"                      ,tempVal);
  params.get("source_fK"                        ,tempVal);
  params.get("source_gamma"                     ,tempVal);
  params.get("source_inclin"                    ,tempVal);
  params.get("source_input_galaxy_file"         ,tempStr);
  params.get("source_mag_limit"                 ,tempVal);
  params.get("source_nuo"                       ,tempVal);
  params.get("source_opening_ang"               ,tempVal);
  params.get("source_r_in"                      ,tempVal);
  params.get("source_r_out"                     ,tempVal);
  params.get("source_sb_limit"                  ,tempVal);
  params.get("z_lens"                           ,tempVal);
  params.get("z_source"                         ,tempVal);
  params.get("zsource_reference"                ,tempVal);
  params.get("field_Nplanes"                    ,tempVal);
  params.get("field_buffer"                     ,tempVal);
  params.get("field_fov"                        ,tempVal);
  params.get("field_input_simulation_center_DEC",tempVal);
  params.get("field_input_simulation_center_RA" ,tempVal);
  params.get("field_input_simulation_format"    ,tempStr);
  params.get("field_input_simulation_path"      ,tempStr);
  params.get("field_input_simulation_radius"    ,tempVal);
  params.get("field_internal_profile"           ,tempStr);
  params.get("field_internal_profile_galaxy"    ,tempVal);
  params.get("field_mass_func_alpha"            ,tempVal);
  params.get("field_mass_func_type"             ,tempVal);
  params.get("field_min_mass"                   ,tempVal);
  params.get("field_off"                        ,tempVal);
  params.get("field_prof_internal_slope_pl"     ,tempVal);
  params.get("field_prof_internal_slope_pnfw"   ,tempVal);
}
*/
