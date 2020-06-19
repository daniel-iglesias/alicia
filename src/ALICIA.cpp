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
#include <unistd.h>
#include <stdlib.h>
#include <math.h>

// Public domain JSON library.
#include "json/json.h"
#include "json/json-forwards.h"
#include "jsoncpp.cpp"

extern "C"{
    #include "ppf.h"
}

#include "simulation.h"
#include "node.h"
#include "common.h"

using namespace std;
using namespace mknix;

int writePPF(int shot){
    /* PPF Docs
    https://users.euro-fusion.org/pages/data-ppf-jpf/ppfuserguide/quickstart.shtml
    */    

    // As we don't know what form we should save the data in, I'm going to just save an arbitrary DataType.
    cout << '\n' << "Saving arbitrary shot " << shot << " data." << '\n';
    
    /* Read
    int NTMAX = 1024;
    int NXMAX = 1;
    int NDMAX = NTMAX * NXMAX;
    
    int ISHOT = 0;
    int ISEQ = 0;
    int ILUN = 0;
    int IERR = 0;
    
    int IRDAT[13];
    int IWDAT[13];
    int INDAT[10];
    
    float DATA[NDMAX];
    float X[NXMAX];
    float T[NTMAX];
    
    char DDA[5];
    char DTYPE[5];
    char IHDAT[60];
    
    // Initialise the PPF system.
    PPFGO(&ISHOT, &ISEQ, &ILUN, &IERR);
    if (IERR != 0){
        cout << "PPFGO error: " << IERR << '\n';
        return EXIT_FAILURE;
	}
	
    // Set up the PPF system to read to a public record.
	PPFUID("JETPPF\0", "R\0", 8, 1);
	
	// Specify the shot we want to read.
	ISHOT = shot;
	PPFGO(&ISHOT, &ISEQ, &ILUN, &IERR);
    if (IERR != 0){
        cout << "PPFGO error: " << IERR << '\n';
        return EXIT_FAILURE;
	}
	
	strcpy(DDA, "MAGN\0");
    strcpy(DTYPE, "IPLA\0");
       
    // First and last X point requested.
    IRDAT[0] = 1;
    IRDAT[1] = NXMAX;

    // First and last T point requested.
    IRDAT[2] = 1;
    IRDAT[3] = NTMAX;

    // Dimensions of DATA X, and T arrays.
    IRDAT[4] = NDMAX;
    IRDAT[5] = NXMAX;
    IRDAT[6] = NTMAX;
       
    // Set flags to return x- and t-vectors.
    IRDAT[7] = 1;
    IRDAT[8] = 1;
	
	// Grab the data.
    PPFGET(&ISHOT, DDA, DTYPE, IRDAT, IHDAT, IWDAT, DATA, X, T, &IERR, 5, 5, 60);
    if (IERR != 0){
	    cout << "PPFGET error: " << IERR << '\n';
        return EXIT_FAILURE;
	}
	
    // Output whatever we desire!
    for (int i = 0; i < 13; i++){
	    cout << i << " " << IRDAT[i] << '\n';
    }
	//*/
	
    //* Write
    int ISHOT = 0;
    int ISEQ = 0;
    int ILUN = 0;
    int IERR = 0;
    
    int ITIME = 0;
    int IDATE = 0;
    
    int IRDAT[13];
    int IWDAT[13];
    int INDAT[10];
    
    float DATA[1];
    
    char DDA[5];
    char DTYPE[5];
    char IHDAT[60];

    // PPF header comment. Empty, but required.
    char ICOM[25];
    
    // Initialise the PPF system.
    PPFGO(&ISHOT, &ISEQ, &ILUN, &IERR);
    if (IERR != 0){
        cout << '\n' << "PPFGO error: " << IERR << '\n';
        return EXIT_FAILURE;
    }
	
    // Set up the PPF system to write to a private record.
    PPFUID("tclayton\0", "W\0", 9, 1);
    
    // Specify the shot we want to write.
    ISHOT = shot;
	PPFGO(&ISHOT, &ISEQ, &ILUN, &IERR);
    if (IERR != 0){
        cout << '\n' << "PPFGO error: " << IERR << '\n';
        return EXIT_FAILURE;
	}
	
	// Create temporary PPF file.
	INDAT[0] = shot;
	INDAT[1] = 123000;  // 12:30:00
	INDAT[2] = 311018;  // 2018/10/31
	INDAT[3] = 0;
	INDAT[9] = 0;
	PPFOPN(&ISHOT, 0, 0, 0, INDAT, ICOM, &IERR, 24);
    if (IERR != 0){
        cout << '\n' << "PPFOPN error: " << IERR << '\n';
	    PPFABO(&IERR);
	    return EXIT_FAILURE;
    }
    
    // Write test data to the temporary PPF file.
    // Document this better when we know what we're actually supposed to be writing.
    IRDAT[5] = 1;
    IRDAT[6] = 1;
    IRDAT[7] = 0;
	IRDAT[8] = 0;
	IRDAT[9] = 0;
	IRDAT[10] = 0;
	IRDAT[11] = 0;
	IRDAT[12] = 0;
	
	char *padding = "     ";
	char *floatPadded = "F   ";
	strcpy(&IHDAT[0], "Tlsa");
	strcpy(&IHDAT[4], padding);
	strcpy(&IHDAT[8], padding);
	strcpy(&IHDAT[12], padding);
	strcpy(&IHDAT[16], padding);
	strcpy(&IHDAT[20], padding);
	strcpy(&IHDAT[24], floatPadded);
	strcpy(&IHDAT[28], padding);
	strcpy(&IHDAT[32], floatPadded);
	strcpy(&IHDAT[36], "                        ");
	
	DATA[0] = 3.141592654;
	
    strcpy(DDA, "TEST\0");
    strcpy(DTYPE, "Tlsa\0");
    PPFWRI(&ISHOT, DDA, DTYPE, IRDAT, IHDAT, IWDAT, DATA, NULL, NULL, &IERR, 5, 5, 60);
    
    if (IERR != 0){
        cout << '\n' << "PPFWRI error: " << IERR << '\n';
	    PPFABO(&IERR);
	    return EXIT_FAILURE;
    }
    
    // Scrap the temporary PPF file. To save, replace with PPFCLO.
    PPFABO(&IERR);
    //*/
    
    return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
    lmx::setMatrixType( 3 );
    lmx::setLinSolverType( 2 );

    if (argc >= 3) {
        
        int dryRuns = -1; // -1 - not a dry run, 0 - dry run, no iterations, >0 - dry run, this many iterations
        int outputFilesDetail = 1; // 0 - nothing, 1 - timing, 2 - All
        int signalFiltering = 1; // If this is 1, then no filtering happens, use bigger ints for more smoothing
        int infoEveryPercentageIteration = 1; // Possible values: 0 to 100
        int initialTemperature = 150; // 0 - Validation & tortures, 200 - 3D tortures (3 & 4), 150 - Reference in-vessel value, 23 - Reference MAST-U
        int verbosity = 2; // Possible values: 0, 1, 2
        int shot = -1; // Possible values: any valid pulse number. If not provided the results will not be saved as a PPF.
        
        // Process optional command line flags
        bool overriddenByConfig = false;
        for(int i=1; i<argc-2; i++){
            if (overriddenByConfig)
                break;
            
            if (argv[i][0] == '-'){
                char flag = argv[i][1];
                try {    
                    switch (flag){
                        case 'c':
                        case 'C': {
                        
                            // Read in config file.
                            Json::Reader reader;
                            Json::Value obj;

                            ifstream ifs(argv[i+1]);
                            if (ifs.fail()){
                                cout << "There was a problem reading the config file. Does it exist, and do you have permisison to read it?" << '\n';
                                return EXIT_FAILURE;
                            }
                            
                            if (!(reader.parse(ifs, obj))){
                                cout << "There was a problem parsing the config file. Is your syntax correct?" << '\n';
                                cout << reader.getFormattedErrorMessages() << '\n';
                                return EXIT_FAILURE;
                            }
                            
                            // Break out of external loop.
                            overriddenByConfig = true;

                            if (!obj["dryRuns"].isNull()){
                                dryRuns = obj["dryRuns"].asUInt();
                            }
                            if (!obj["outputFilesDetail"].isNull())
                                outputFilesDetail = obj["outputFilesDetail"].asUInt();
                                
                            if (!obj["signalFiltering"].isNull())
                                signalFiltering = obj["signalFiltering"].asUInt();
                                
                            if (!obj["infoEveryPercentageIteration"].isNull())
                                infoEveryPercentageIteration = obj["infoEveryPercentageIteration"].asUInt();
                                
                            if (!obj["initialTemperature"].isNull())
                                initialTemperature = obj["initialTemperature"].asUInt();
                                
                            if (!obj["verbosity"].isNull())
                                verbosity = obj["verbosity"].asUInt();
                                
                            if (!obj["shot"].isNull())
                                shot = obj["shot"].asUInt();

                            break;
                        }
                        case 'd':
                        case 'D':
                            dryRuns = std::stoi(argv[i+1]);
                            break;
                        case 'o':
                        case 'O':
                            outputFilesDetail = std::stoi(argv[i+1]);
                            break;
                        case 'f':
                        case 'F':
                            signalFiltering = std::stoi(argv[i+1]);
                            break;
                        case 'i':
                        case 'I':
                            infoEveryPercentageIteration = std::stoi(argv[i+1]);
                            break;
                        case 't':
                        case 'T':
                            initialTemperature = std::stoi(argv[i+1]);
                            break;
                        case 'v':
                        case 'V':
                            verbosity = std::stoi(argv[i+1]);
                            break;
                        case 's':
                        case 'S':
                            shot = std::stoi(argv[i+1]);
                            break;
                    }
                }
                catch (const std::invalid_argument& exception) {
                  cout << "\nError at flag: -" << flag << ": the provided argument is invalid!" << '\n';
                  return EXIT_FAILURE;
                }
                catch (const std::out_of_range& exception) {
                  cout << "\nError at flag: -" << flag << ": the provided argument is too big!" << '\n';
                  return EXIT_FAILURE;
                }
            }
        }
        
        //     cout << "EXECUTING: $ alicia " << argv[1] << endl;
        //     system("echo -n '1. Current Directory is '; pwd");
        mknix::Simulation mySimulation;
        mySimulation.setOutputFilesDetail(outputFilesDetail); 
        
        // Only output the simulation building at high verbosity.
        // As we don't want to alter the MkniX library, we do this by making cout silently discard any output.
        if (verbosity < 2)
            cout.setstate(ios_base::failbit);
        mySimulation.inputFromFile(argv[argc-2]);
        // Stop suppressing cout.
        if (verbosity < 2)
            cout.clear();

        std::string systemName = "tiles";
        std::ifstream ifile(argv[argc-1]);
        //         ifile.open("../Pulse_data/89162.9ATP.TPS/89162.9ATP.TPSO.49-40000_steps.txt");
        //         ifile.open("../Pulse_data/89445.9ATP.TPS/89445.9ATP.TPSO.51-40000_steps.txt");

        boxFIR box(signalFiltering); 

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

        //         cout << "t0 = " << vectors[2][0] << endl;
        //         cout << "t1000 = " << vectors[1002][0] << endl;
        //         cout << "size of t = " << vectors.size() << endl;

        double t_0 = vectors[2][0];
        double time = t_0;
        double t_n = vectors.back()[0];
        int iterations = (vectors.size()-3);
        double delta_t = (t_n - t_0)/iterations;
        int infoEveryNthIteration = iterations + 1;
        if (infoEveryPercentageIteration > 0){
            infoEveryNthIteration = fmax(1, iterations/ (100/infoEveryPercentageIteration));
        }
        cout << "t_0 = " << t_0 << endl;
        cout << "t_n = " << t_n << endl;
        cout << "iterations = " << iterations << endl;
        cout << "delta_t = " << delta_t << endl;
        
        if (dryRuns > -1)
            iterations = dryRuns;
            if (iterations < 1)
                return EXIT_SUCCESS;

        // initializing coordinate map:
        map<double, double> input_load_step;
        int input_points_num = vectors[0].size();
        for(int i=0; i<input_points_num; ++i){
            input_load_step[ vectors[0][i] ]= 0.;
            interface_file << vectors[0][i] << "\t";
        }
        interface_file << endl;
        mySimulation.setInitialTemperatures(initialTemperature);
        mySimulation.init(verbosity);

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
            if( infoEveryPercentageIteration > 0 && i % infoEveryNthIteration == 0 ){
                cout << int(double(i)/iterations*100) << " % Complete," << "\t TIME = " << time << ",\t STEP = " << i << endl;
            }
        }
        mySimulation.endSimulation();
        
        
        // Output the settings used during this run.
        Json::StyledStreamWriter writer;
        Json::Value obj;

        ofstream ofs("run_settings.json");
        
        obj["files"]["MkniX"] = argv[argc-2];
        obj["files"]["temperatureMap"] = argv[argc-1];
        
        obj["external"]["dryRuns"] = dryRuns;
        obj["external"]["outputFilesDetail"] = outputFilesDetail;
        obj["external"]["signalFiltering"] = signalFiltering;
        obj["external"]["infoEveryPercentageIteration"] = infoEveryPercentageIteration;
        obj["external"]["initialTemperature"] = initialTemperature;
        obj["external"]["verbosity"] = verbosity;
        obj["external"]["shot"] = shot;
        
        obj["internal"]["numberOfInterfaceNodes"] = numberOfInterfaceNodes;
        obj["internal"]["t_0"] = t_0;
        obj["internal"]["t_n"] = t_n;
        obj["internal"]["iterations"] = iterations;
        obj["internal"]["delta_t"] = delta_t;
        
        writer.write(ofs, obj);
        
        if (dryRuns < 0 && shot > 0){
            int savingFailed = writePPF(shot);
            if (savingFailed)
                cout << '\n' << "Saving failed." << '\n';
            else
                cout << '\n' << "Saving success." << '\n';
            return savingFailed;   
        }


    } else if(argc != 2) {
        cout<< "Need two parameters: 1. mknix input file, and 2. Temperature map."<<endl
        << "Usage: $./alicia \"mknix_file_name\" \"IR_temps_file\""<<endl;
    }
    return EXIT_SUCCESS;
}

