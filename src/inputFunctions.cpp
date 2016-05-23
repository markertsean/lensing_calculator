#include <cstring>
#include <slsimlib.h>
#include <lensing_classes.h>


/*Read parameters not included in paramfile (Nbins, inpfile, etc.)
//
//input:
//     userInfo class object to overwrite
//
//output:
//     void, overwrite inpInfo
//
*/
void ReadInpFile( userInfo &inpInfo, std::string inputFile ){
  FILE *pFile;
  char   inpC1[35],inpC2[35];

  // Attempt to open file, if successful go line by line
  //  finding the variable name and value

  pFile = fopen(inputFile.c_str(),"r");

  if (pFile!=NULL){
    while ( fscanf(pFile,"%s%s",inpC1,inpC2) != EOF ){
      std::string inpS = std::string(inpC1);
      if      ( inpS=="angFOV"      ){        inpInfo.angFOV      = atof(inpC2);      }
      else if ( inpS=="N_pixels"    ){        inpInfo.N_pixels    = atoi(inpC2);      }
      else if ( inpS=="N_bins"      ){        inpInfo.N_bins      = atoi(inpC2);      }
      else if ( inpS=="N_sources"   ){        inpInfo.N_sources   = atoi(inpC2);      }
      else if ( inpS=="N_particles" ){        inpInfo.N_particles = atoi(inpC2);      }
      else if ( inpS=="N_partSmooth"){        inpInfo.N_partSmooth= atoi(inpC2);      }
      else if ( inpS=="num_threads" ){        inpInfo.num_threads = atoi(inpC2);      }
      else if ( inpS=="model"       ){        inpInfo.model       = atoi(inpC2);      }
      else if ( inpS=="nbody"       ){        inpInfo.nbody       = atoi(inpC2);      }
      else if ( inpS=="massMap"     ){        inpInfo.nbody       = atoi(inpC2);      }
      else if ( inpS=="R_max"       ){        inpInfo.R_max       = atof(inpC2);      }
      else if ( inpS=="readFile"    ){        inpS                = std::string(inpC2);        inpInfo.readFile    =      inpS  ;      }
      else{
          std::cerr << " Couldn't recognize input from " << inputFile <<
                     ": " << inpS << std::endl << std::endl;
        exit(0);
      }
    }
    fclose(pFile);
  }
  else{
    std::cerr << "Couldn't find file: lensuserParams.dat" << std::endl;
    exit(0);
  }

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
