/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

// declare cell definitions here 

Cell_Definition custom_cell;
double hif_base_concentration;

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	hif_base_concentration = parameters.doubles("hif_base_concentration");

	initialize_default_cell_definition(); 
	
	create_default_cell_definition();

	create_custom_cell_definition();

	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
/* now this is in XML 
	default_microenvironment_options.X_range = {-1000, 1000}; 
	default_microenvironment_options.Y_range = {-1000, 1000}; 
	default_microenvironment_options.simulate_2D = true; 
*/
	
	// make sure to override and go back to 2D 
	// if( default_microenvironment_options.simulate_2D == false )
	// {
	// 	std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
	// 	default_microenvironment_options.simulate_2D = true; 
	// }
	
/* now this is in XML 	
	// no gradients need for this example 

	default_microenvironment_options.calculate_gradients = false; 
	
	// set Dirichlet conditions 

	default_microenvironment_options.outer_Dirichlet_conditions = true;
	
	// if there are more substrates, resize accordingly 
	std::vector<double> bc_vector( 1 , 38.0 ); // 5% o2
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
	
	// set initial conditions 
	default_microenvironment_options.initial_condition_vector = { 38.0 }; 
*/
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	// create some cells near the origin
	
	double cell_radius = cell_defaults.phenotype.geometry.radius;
	double spacing = 2.0*cell_radius;
	double tissue_radius = parameters.doubles("tissue_radius");

	// create_circular_tissue(-500.0, 0.0, cell_defaults, tissue_radius, spacing);
	// create_circular_tissue(500.0, 0.0, custom_cell, tissue_radius, spacing);
	create_circular_tissue(0.0, 0.0, custom_cell, tissue_radius, spacing);

	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring 
	
	std::vector<std::string> output = false_cell_coloring_cytometry(pCell); 
	int hif_index = pCell->custom_data.find_variable_index("hif_concentration");
		
	if( pCell->phenotype.death.dead == false && pCell->type == 1 )
	{
		 output[0] = "black"; 
		 std::string color = "rgb(255, 0, 0)";
		 if ( pCell->custom_data[hif_index] > hif_base_concentration) {
			 color = "rgb(0,0,255)";
		 }
		 output[2] = color; 
	}
	
	return output; 
}

void create_default_cell_definition(void)
{
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// Name the default cell type 
	
	cell_defaults.type = 0; 
	cell_defaults.name = "tumor cell"; 
	
	// set default cell cycle model 

	cell_defaults.functions.cycle_model = flow_cytometry_separated_cycle_model; 
	
	// set default_cell_functions; 
	
	cell_defaults.functions.update_phenotype = update_cell_and_death_parameters_O2_based; 
	
	// needed for a 2-D simulation: 
	
	/* grab code from heterogeneity */ 
	
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment );	
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	// set the rate terms in the default phenotype 

	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 

	int G0G1_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );
	int S_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::S_phase );

	// initially no necrosis 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 

	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 10; 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 38;

	cell_defaults.custom_data.add_variable("hif_concentration", "dimensionless", hif_base_concentration);
	cell_defaults.custom_data.add_variable("ldh_level", "dimensionless", 0.0);
	cell_defaults.custom_data.add_variable("pdh_level", "dimensionless", 0.0);
	cell_defaults.custom_data.add_variable("pdk_level", "dimensionless", 0.0);
	
}

void create_custom_cell_definition(void)
{
	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 

	// Now, let's define another cell type. 
	// It's best to just copy the default and modify it. 
	
	// make this cell type randomly motile, less adhesive, greater survival, 
	// and less proliferative 
	
	custom_cell = cell_defaults; 
	custom_cell.type = 1; 
	custom_cell.name = "test cell";

	// make sure the new cell type has its own reference phenotype
	custom_cell.parameters.pReference_live_phenotype = &( custom_cell.phenotype );  

	custom_cell.functions.update_phenotype = my_custom_phenotype_update;
	custom_cell.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 10;
	custom_cell.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0;
	custom_cell.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 38;
	custom_cell.phenotype.death.rates[necrosis_model_index] = 0.0;
	// test_cell.phenotype.death.rates[apoptosis_model_index] = 0.0;

	// custom_cell.parameters.o2_hypoxic_threshold = parameters.doubles("o2_hypoxic_threshold");
	// custom_cell.parameters.o2_hypoxic_response = parameters.doubles("o2_hypoxic_response");
	// custom_cell.parameters.o2_hypoxic_saturation = parameters.doubles("o2_hypoxic_saturation");
	
	//custom_cell.functions.custom_cell_rule = simulate_metabolism;

	//Custom data and funtions	
	// custom_cell.custom_data.add_variable("hif_concentration", "dimensionless", hif_base_concentration);
	// custom_cell.custom_data.add_variable("ldh_level", "dimensionless", 0.0);
	// custom_cell.custom_data.add_variable("pdh_level", "dimensionless", 0.0);
	// custom_cell.custom_data.add_variable("pdk_level", "dimensionless", 0.0);
}

void create_circular_tissue (double center_x, double center_y, Cell_Definition& cell_definition, double tissue_radius, double spacing)
{
	Cell* pC;

	int i = 0;
	int j = 0;
	double z = 0.0;
	double dist = 0.0;
	
	while (dist < tissue_radius)
	{
		while ( dist < tissue_radius)
		{
			pC = create_cell(cell_definition);
			pC -> assign_position( (center_x+i*spacing), (center_y+j*spacing), z);
			
			pC = create_cell(cell_definition);
			pC -> assign_position( (center_x+i*spacing), (center_y-j*spacing), z);
			
			
			pC = create_cell(cell_definition);
			pC -> assign_position( (center_x-i*spacing), (center_y+j*spacing), z);
			
			pC = create_cell(cell_definition);
			pC -> assign_position( (center_x-i*spacing), (center_y-j*spacing), z);

			j += 1;
			dist = hypot( center_x-(center_x+i*spacing), center_y-(center_y+j*spacing) );
		}
		i += 1;
		j = 0;
		dist = hypot( center_x-(center_x+i*spacing), center_y-(center_y+j*spacing) );
	}
}

void my_custom_phenotype_update( Cell* pCell, Phenotype& phenotype, double dt )
{
	if (pCell->phenotype.death.dead == true) { return; }
	update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);
	simulate_metabolism(pCell, phenotype, dt);
}

void my_custom_uptake_rates_update(Cell* pCell)
{

}

void my_custom_secretion_rates_update(Cell* pCell)
{

}

void simulate_metabolism(Cell* pCell, Phenotype& phenotype, double dt)
{
	compute_hif_concentration(pCell, phenotype, dt);
	compute_ldh_concentration(pCell, phenotype, dt);
	compute_pdh_concentration(pCell, phenotype, dt);
	my_custom_uptake_rates_update(pCell);
	my_custom_secretion_rates_update(pCell);
}

void compute_hif_concentration(Cell* pCell, Phenotype& phenotype, double dt)
{
	int o2_index = microenvironment.find_density_index("oxygen");
	double pO2 = pCell->nearest_density_vector()[o2_index];
	
	int hif_index = pCell->custom_data.find_variable_index("hif_concentration");

	if ( pO2 >= pCell->parameters.o2_hypoxic_response ) {
		pCell->custom_data[hif_index] = hif_base_concentration;
	} 
	else {
		pCell->custom_data[hif_index] = hif_base_concentration*exp( 2.5*(1-(pO2/760.0) ) );
	}
}

void compute_ldh_concentration(Cell* pCell, Phenotype& phenotype, double dt)
{

}

void compute_pdh_concentration(Cell* pCell, Phenotype& phenotype, double dt)
{

}
