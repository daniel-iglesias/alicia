/***************************************************************************
 *   Copyright (C) 2007-2017 by Daniel Iglesias                            *
 *   Copyright (c) 2013-2018, European Atomic Energy Community (EURATOM)   *
 *   https://github.com/daniel-iglesias/alicia                             *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 2 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU Lesser General Public      *
 *   License along with this program; if not, write to the                 *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
//#include <unistd.h>

#include "simulation.h"
#include "node.h"
#include "common.h"

using namespace std;
using namespace mknix;

int main(int argc, char *argv[])
{
    lmx::setMatrixType( 3 );
    lmx::setLinSolverType( 2 );

    if (argc >= 3 && argc < 5) {
        //     cout << "EXECUTING: $ alicia " << argv[1] << endl;
        //     system("echo -n '1. Current Directory is '; pwd");
        mknix::Simulation mySimulation;
        mySimulation.setOutputFilesDetail( 1 ); // 0 - nothing, 1 - timing, 2 - All
        mySimulation.inputFromFile(argv[1]);

        std::string systemName = "tiles";
        std::ifstream ifile(argv[2]);

        boxFIR box(1); //If this is 1, then no filtering happens, use bigger ints for more smoothing

        std::ofstream interface_file ("interface_coordinates.txt");
        int numberOfInterfaceNodes = mySimulation.getInterfaceNumberOfNodes(systemName);
        vector<double> interface_coordinates;
        // Fill x coordinates for constrained nodes...
        for(int i=0; i<numberOfInterfaceNodes; ++i){
            interface_coordinates.push_back(mySimulation.getInterfaceNode(systemName, i)->getX());
            interface_file << interface_coordinates.back() << "\t";
        }
        interface_file << endl;
        cout << "NODES: " << numberOfInterfaceNodes << endl;
        cout << interface_coordinates[0] << " ... " << interface_coordinates.back() << endl;

        const std::vector<std::vector<double> > vectors = read_lines(ifile) ;

        double t_0 = vectors[2][0];
        double time = t_0;
        double t_n = vectors.back()[0];
        int iterations = (vectors.size()-3);
        double delta_t = (t_n - t_0)/iterations;
        int infoEveryNthIteration = iterations/100+1; // Every 1 %
        cout << "t_0 = " << t_0 << endl;
        cout << "t_n = " << t_n << endl;
        cout << "iterations = " << iterations << endl;
        cout << "delta_t = " << delta_t << endl;
        if (argc >= 4){
            iterations =2; //doing dry runs to check time parameters
            if (argv[3] == "dry") return EXIT_SUCCESS;
        }

        // initializing coordinate map:
        map<double, double> input_load_step;
        int input_points_num = vectors[0].size();
        for(int i=0; i<input_points_num; ++i){
            input_load_step[ vectors[0][i] ]= 0.;
            interface_file << vectors[0][i] << "\t";
        }
        interface_file << endl;
//         mySimulation.setInitialTemperatures(0); // Validation & tortures
//        mySimulation.setInitialTemperatures(200); // 3D tortures (3 & 4)
        mySimulation.setInitialTemperatures(150); // Reference in-vessel value
//         mySimulation.setInitialTemperatures(23); // Reference MAST-U
        mySimulation.init(2); // Output Vervosity=0

        std::ofstream temps_filtered ("temps_filtered.txt");
        std::ofstream temps_interface ("temps_interface.txt");
        std::ofstream flux_interface ("flux_interface.txt");
        std::ofstream flux_max_min ("flux_interface_max_min.txt");

        // Output 2D files m-type matrix preparation: size of time vector in first element, x-coordinates in the rest of the 1st row.
        temps_interface << iterations << "\t";
        flux_interface << iterations << "\t";
        for(int i=0; i<numberOfInterfaceNodes; ++i){
            temps_interface << mySimulation.getInterfaceNode(systemName, i)->getX() << "\t";
            flux_interface << mySimulation.getInterfaceNode(systemName, i)->getX() << "\t";
        }
        temps_interface << endl;
        flux_interface << endl;

        //         double factor(1./0.00387); //width of the element
        //         factor*=0.836695; // power equilibrium between flat and 3D shadowed as measured in CAD; nominal position
        //         factor*=0.85; // Reduction of thermal capacity due to the gap (less material to heat)
        //         factor*=0.9; // Effect of the change in section of the lamella; the mid part of the lamella is thinner.

        for(int i=0; i<iterations; ++i){
            vector<double> filtered_vector ( vectors[i+2] );
            filtered_vector.erase(filtered_vector.begin()); // Remove time information
            vector<double> interface_fluxes( interface_coordinates.size() );

             box.filter(filtered_vector);
            //             cout << "input_points_num = " << input_points_num << endl;
            //             cout << "input_load_step.size() = " << input_load_step.size() << endl;
            //             cout << "filtered_vector.size() = " << filtered_vector.size() << endl;
            for(int j=0; j<input_points_num; ++j){
                input_load_step.at( vectors[0][j] ) = filtered_vector[j];
                temps_filtered << input_load_step.at( vectors[0][j] ) << '\t';
            }
            temps_filtered << '\n';

            temps_interface << time << '\t';
            flux_interface << time << '\t';

            for(int j=0; j<numberOfInterfaceNodes; ++j){
                mySimulation.getInterfaceNode(systemName,j)->setqt( interpolate1D( interface_coordinates[j], input_load_step ) );
                temps_interface << mySimulation.getInterfaceNode(systemName,j)->getqt() << '\t';
            }
            temps_interface << '\n';
            mySimulation.solveStep();
            // First node:
            {
                std::string constraintName( "INT_" + std::to_string(0) );
                interface_fluxes.push_back( mySimulation.getConstraintOutput(constraintName, systemName, 0)
                / ( std::abs(interface_coordinates[1]-interface_coordinates[0]) / 2. ) ) ; // scaled with the global Qperp factor and border condition
                flux_interface << interface_fluxes.back() << '\t';
            }

            for(int j=1; j<numberOfInterfaceNodes-1; ++j){
                std::string constraintName( "INT_" + std::to_string(j) );
                //                 cout << "Getting constraint data from " << systemName << "." << constraintName << endl;
                interface_fluxes.push_back( mySimulation.getConstraintOutput(constraintName, systemName, 0)
                / ( std::abs(interface_coordinates[j+1]-interface_coordinates[j-1]) / 2.) ) ;
                flux_interface << interface_fluxes.back() << '\t';
            }
            // Last node:
            {
                std::string constraintName( "INT_" + std::to_string(numberOfInterfaceNodes-1) );
                interface_fluxes.push_back( mySimulation.getConstraintOutput(constraintName, systemName, 0)
                / ( std::abs(interface_coordinates.back()-interface_coordinates[numberOfInterfaceNodes-2]) / 2. ) ) ; // scaled with the global Qperp factor and border condition
                flux_interface << interface_fluxes.back() << '\t';
            }

            flux_interface << '\n';
            flux_max_min << time << '\t'
            << *max_element(std::begin(interface_fluxes), std::end(interface_fluxes))  << '\t'
            << *min_element(std::begin(interface_fluxes), std::end(interface_fluxes))
            << '\n';
            time += delta_t;
            if( i % infoEveryNthIteration == 0 ){
                cout << int(double(i)/iterations*100) << " % Complete," << "\t TIME = " << time << ",\t STEP = " << i << endl;
            }
        }
        mySimulation.endSimulation();

    } else if(argc != 2) {
        cout<< "Need two parameters: 1. mknix input file, and 2. Temperature map."<<endl
        << "Usage: $./alicia \"mknix_file_name\" \"IR_temps_file\""<<endl;
    }
    return EXIT_SUCCESS;
}
