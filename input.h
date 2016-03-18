//
//  input.h
//  
//
//  Created by Markus Scheucher on 17/02/16.
//
//

#ifndef __input_h__
#define __input_h__

//-- END time --
#define tmax (double) 300       //endtime

//-- simulation box --
#define nx (int) 240            //gridsize_X
#define ny (int) 800            //gridsize_Y
#define boxx (double) 12        //length_X
#define boxy (double) 40        //length_Y
#define na (int) 1              //number_of_boundary_layers
#define thick (double) 20.      //boundary_layer_thickness

//-- boundary consitions --
#define x_boundaries (int) 1	//0:transmissive_1:periodic
#define y_boundaries (int) 0	//0:transmissive_1:periodic

//-- simulation methods --
#define cleaning (int) 1        //0:no_cleaning_1:fieldCD_2:fluxCT_3:fluxCD
#define cmethod (int) 0         //0:Barmin_et_al_1:max_cmax_over_grid
#define nout (int) 60           //number_of_output_files_per_parameter
#define limiter (int) 0         //0:Woodward_1:MinMod

//-- plasma parameters --
#define cfl_input (double) 0.8	//courant_number
#define kappa (double) 1.666667	//kappa
#define mu (double) 1.0         //mu_0

//-- sinusoidal disturbance --
#define kx (double) 0.52        //Wavenumber_of_Disturbance_in_Boundary_Layer
#define dist (double) 0.01      //Amplitude_of_Disturbance_in_Boundary_Layer

//-- solar wind parameters --
#define iR0 (double) 0.1        //START_Density-LOWspeed
#define iVx0 (double) 0.5       //START_X-Velocity-LOWspeed
#define iVy0 (double) 0.        //START_Y-Velocity-LOWspeed
#define iVz0 (double) 0.        //START_Z-Velocity-LOWspeed
#define iP0 (double) 0.0268     //START_Pressure-LOWspeed
#define iBx0 (double) 0.        //START_X-Magnetic-Field-LOWspeed
#define iBy0 (double) 0.        //START_Y-Magnetic-Field-LOWspeed
#define iBz0 (double) 0.49376	//START_Z-Magnetic-Field-LOWspeed

//-- filament parameters --
#define iR1 (double) 1.         //START_Density-HIGHspeed
#define iVx1 (double) 1.        //START_X-Velocity-HIGHspeed
#define iVy1 (double) 0.        //START_Y-Velocity-HIGHspeed
#define iVz1 (double) 0.        //START_Z-Velocity-HIGHspeed
#define iP1 (double) 0.0268     //START_Pressure-HIGHspeed
#define iBx1 (double) 0.        //START_X-Magnetic-Field-HIGHspeed
#define iBy1 (double) 0.        //START_Y-Magnetic-Field-HIGHspeed
#define iBz1 (double) 0.49376	//START_Z-Magnetic-Field-HIGHspeed

//-- low pressure handling --
#define plim (double) 0.         //minimum allowed local thermal pressure
#define pset (double) 0.01       //SET pressure value if p < plim

#endif
