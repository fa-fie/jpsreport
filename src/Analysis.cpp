/**
 * \file        Analysis.cpp
 * \date        Oct 10, 2014
 * \version     v0.7
 * \copyright   <2009-2016> Forschungszentrum Juelich GmbH. All rights reserved.
 *
 * \section License
 * This file is part of JuPedSim.
 *
 * JuPedSim is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * JuPedSim is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with JuPedSim. If not, see <http://www.gnu.org/licenses/>.
 *
 * \section Description
 * The Analysis class represents a process of analyzing groups of pedestrian
 * trajectories from experiment or simulation. Different measurement methods
 * can be used and are defined by various parameters and functions.
 *
 *
 **/

#include "Analysis.h"

#include "general/Logger.h"
#include "methods/Method_A.h"
#include "methods/Method_B.h"
#include "methods/Method_C.h"
#include "methods/Method_D.h"
#include "methods/Method_E.h"
#include "methods/Method_F.h"
#include "methods/Method_G.h"
#include "methods/Method_H.h"
#include "methods/PedData.h"
#include "methods/VoronoiDiagram.h"

#include <algorithm> // std::min_element, std::max_element
#include <cfloat>
#include <fmt/format.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <vector>

#ifdef __linux__
#include <dirent.h>
#include <sys/stat.h>
#elif __APPLE__
#include <dirent.h>
#include <sys/stat.h>
#else
#include <direct.h>
#endif


using boost::geometry::dsv;
using namespace std;

/************************************************
 // Konstruktoren
 ************************************************/

Analysis::Analysis()
{
    _projectRootDir = "";
    _deltaF = 5; // half of the time interval that used to calculate instantaneous velocity of ped
                 // i. Here v_i = (X(t+deltaF) - X(t+deltaF))/(2*deltaF).   X is location.
    _DoesUseMethodA = false; // Method A (Zhang2011a)
    _DoesUseMethodB = false; // Method B (Zhang2011a)
    _DoesUseMethodC = false; // Method C //calculate and save results of classic in separate file
    _DoesUseMethodD = false; // Method D--Voronoi method
    _DoesUseMethodE = false; // Method E
    _DoesUseMethodF = false; // Method F
    _DoesUseMethodG = false; // Method G
    _DoesUseMethodH = false; // Method H

    _vComponent =
        "B"; // to mark whether x, y or x and y coordinate are used when calculating the velocity
    _IgnoreBackwardMovement = false;
    _lowVertexX             = 0;  // LOWest vertex of the geometry (x coordinate)
    _lowVertexY             = 0;  //  LOWest vertex of the geometry (y coordinate)
    _highVertexX            = 10; // Highest vertex of the geometry
    _highVertexY            = 10;
    _trajFormat             = FileFormat::FORMAT_PLAIN;
}

// file.txt ---> file
std::string Analysis::GetBasename(const std::string & str)
{
    unsigned found = str.find_last_of(".");
    return str.substr(0, found);
}
// c:\\windows\\winhelp.exe ---> winhelp.exe
std::string Analysis::GetFilename(const std::string & str)
{
    unsigned found = str.find_last_of("/\\");
    return str.substr(found + 1);
}

void Analysis::InitArgs(ArgumentParser * args)
{
    _geometry = args->GetGeometry();
    auto box  = GetBoundingBox(_geometry, 10);
    boost::geometry::convert(box, _boundingBox);
    LOG_INFO(
        "Bounding box: \n \t\tminX = {:.2f}\n \t\tmaxX = {:.2f} \n \t\tminY = {:.2f} \n\t\tmaxY = "
        "{:.2f}",
        box.min_corner().x() * CMtoM,
        box.max_corner().x() * CMtoM,
        box.min_corner().y() * CMtoM,
        box.max_corner().y() * CMtoM);

    // TODO Same default MA can be used for profiles and globalIFD, currently only the polygon is
    // reused
    if(args->GetIsMethodA()) {
        _DoesUseMethodA                  = true;
        vector<int> Measurement_Area_IDs = args->GetAreaIDforMethodA();
        for(unsigned int i = 0; i < Measurement_Area_IDs.size(); i++) {
            _areasForMethodA.push_back(dynamic_cast<MeasurementArea_L *>(
                args->GetMeasurementArea(Measurement_Area_IDs[i])));
        }
        _deltaT = args->GetTimeIntervalA();
    }

    if(args->GetIsMethodB()) {
        _DoesUseMethodB                  = true;
        vector<int> Measurement_Area_IDs = args->GetAreaIDforMethodB();
        for(unsigned int i = 0; i < Measurement_Area_IDs.size(); i++) {
            MeasurementArea_B * area = dynamic_cast<MeasurementArea_B *>(
                args->GetMeasurementArea(Measurement_Area_IDs[i]));
            if(area->_poly.outer().empty()) {
                LOG_WARNING("Measurement {} has 0 points, will be skipped.", area->_id);
            } else {
                _areasForMethodB.push_back(area);
            }
        }
    }

    if(args->GetIsMethodC()) {
        _DoesUseMethodC                  = true;
        vector<int> Measurement_Area_IDs = args->GetAreaIDforMethodC();
        for(unsigned int i = 0; i < Measurement_Area_IDs.size(); i++) {
            MeasurementArea_B * area = dynamic_cast<MeasurementArea_B *>(
                args->GetMeasurementArea(Measurement_Area_IDs[i]));
            if(area->_poly.outer().empty()) {
                LOG_WARNING("Measurement {} has 0 points, will be skipped.", area->_id);
            } else {
                _areasForMethodC.push_back(area);
            }
        }
    }

    if(args->GetIsMethodD()) {
        _DoesUseMethodD = true;
        // TODO[DH]: modernize loops
        for(unsigned int i = 0; i < args->_configDataD.areaIDs.size(); i++) {
            MeasurementArea_B * area = dynamic_cast<MeasurementArea_B *>(
                args->GetMeasurementArea(args->_configDataD.areaIDs[i]));
            if(area->_poly.outer().empty()) {
                area->_poly = _boundingBox;
            }
            _areasForMethodD.push_back(area);
        }
        _geoPolyMethodD = GetRoomForMeasurementArea(_areasForMethodD);
    }

    if(args->GetIsMethodE()) {
        _DoesUseMethodE                  = true;
        vector<int> Measurement_Area_IDs = args->GetAreaIDforMethodE();
        vector<int> Line_IDs             = args->GetLineIDforMethodE();
        // GetAreaIDforMethodE() and GetLineIDforMethodE() should have the same size
        for(unsigned int i = 0; i < Measurement_Area_IDs.size(); i++) {
            _areasForMethodE.push_back(dynamic_cast<MeasurementArea_B *>(
                args->GetMeasurementArea(Measurement_Area_IDs[i])));
            _linesForMethodE.push_back(
                dynamic_cast<MeasurementArea_L *>(args->GetMeasurementArea(Line_IDs[i])));
        }
        _deltaTMethodE = args->GetTimeIntervalE();
    }

    if(args->GetIsMethodF()) {
        _DoesUseMethodF                  = true;
        vector<int> Measurement_Area_IDs = args->GetAreaIDforMethodF();
        vector<int> Line_IDs             = args->GetLineIDforMethodF();
        // GetAreaIDforMethodF() and GetLineIDforMethodF() should have the same size
        for(unsigned int i = 0; i < Measurement_Area_IDs.size(); i++) {
            _areasForMethodF.push_back(dynamic_cast<MeasurementArea_B *>(
                args->GetMeasurementArea(Measurement_Area_IDs[i])));
            _linesForMethodF.push_back(
                dynamic_cast<MeasurementArea_L *>(args->GetMeasurementArea(Line_IDs[i])));
        }
        _deltaTMethodF = args->GetTimeIntervalF();
    }

    if(args->GetIsMethodG()) {
        _DoesUseMethodG                  = true;
        vector<int> Measurement_Area_IDs = args->GetAreaIDforMethodG();
        for(unsigned int i = 0; i < Measurement_Area_IDs.size(); i++) {
            _areasForMethodG.push_back(dynamic_cast<MeasurementArea_B *>(
                args->GetMeasurementArea(Measurement_Area_IDs[i])));
        }
        _deltaTMethodG         = args->GetTimeIntervalG();
        _dtMethodG             = args->GetDtMethodG();
        _numberPolygonsMethodG = args->GetNumPolyMethodG();
        _pointsMethodG         = args->GetPointsMethodG();
    }

    if(args->GetIsMethodH()) {
        _DoesUseMethodH                  = true;
        vector<int> Measurement_Area_IDs = args->GetAreaIDforMethodH();
        for(unsigned int i = 0; i < Measurement_Area_IDs.size(); i++) {
            _areasForMethodH.push_back(dynamic_cast<MeasurementArea_B *>(
                args->GetMeasurementArea(Measurement_Area_IDs[i])));
        }
        _deltaTMethodH = args->GetTimeIntervalH();
    }

    _deltaF                 = args->GetDelatT_Vins();
    _vComponent             = args->GetVComponent();
    _IgnoreBackwardMovement = args->GetIgnoreBackwardMovement();
    _geometryFileName       = args->GetGeometryFilename();
    _projectRootDir         = args->GetProjectRootDir();
    _trajFormat             = args->GetFileFormat();
    _outputLocation         = args->GetOutputLocation();

    configData_D = args->_configDataD;
}


std::map<int, polygon_2d>
Analysis::GetRoomForMeasurementArea(const std::vector<MeasurementArea_B *> & areas)
{
    std::map<int, polygon_2d> geoPoly;

    // loop over all areas
    for(auto && area : areas) {
        // search for the subroom that contains that area
        for(auto && room : _geometry) {
            point_2d point(0, 0);
            boost::geometry::centroid(area->_poly, point);

            if(boost::geometry::within(point, room)) {
                geoPoly.emplace(area->_id, room);
            }
        }

        if(geoPoly.count(area->_id) == 0) {
            throw std::runtime_error(
                fmt::format(FMT_STRING("No polygon containing the measurement id {}."), area->_id));
        }
    }

    // Get min/max values of all used rooms containing a measurement area
    std::vector<polygon_2d> rooms;
    std::for_each(
        std::begin(geoPoly), std::end(geoPoly), [&](const std::pair<int, polygon_2d> & element) {
            rooms.push_back(element.second);
        });

    auto box = GetBoundingBox(rooms);

    // These values are used for the grid when computing profiles
    _highVertexX = box.max_corner().x();
    _highVertexY = box.max_corner().y();
    _lowVertexX  = box.min_corner().x();
    _lowVertexY  = box.min_corner().y();

    return geoPoly;
}

boost::geometry::model::box<point_2d>
Analysis::GetBoundingBox(const std::vector<polygon_2d> & polygons, double extension)
{
    // Union all rooms
    boost::geometry::model::multi_polygon<polygon_2d> border;

    for(polygon_2d poly : polygons) {
        boost::geometry::model::multi_polygon<polygon_2d> tmp;
        union_(border, poly, tmp);
        border = tmp;
    }

    // Compute bounding box
    boost::geometry::model::box<point_2d> box;
    boost::geometry::envelope(border, box);

    // Slightly increase bounding box
    box.max_corner().x(box.max_corner().x() + extension * M2CM);
    box.max_corner().y(box.max_corner().y() + extension * M2CM);

    box.min_corner().x(box.min_corner().x() - extension * M2CM);
    box.min_corner().y(box.min_corner().y() - extension * M2CM);

    return box;
}

int Analysis::RunAnalysis(const fs::path & filename, const fs::path & path)
{
    PedData data;
    if(data.ReadData(
           _projectRootDir,
           _outputLocation,
           path,
           filename,
           _trajFormat,
           _deltaF,
           _vComponent,
           _IgnoreBackwardMovement) == false) {
        LOG_ERROR("Could not parse the file {}", filename);
        return EXIT_FAILURE;
    }

    //-----------------------------check whether there is pedestrian outside the whole
    // geometry--------------------------------------------
    std::map<int, std::vector<int>> _peds_t = data.GetPedIDsByFrameNr();
    for(int frameNr = 0; frameNr < data.GetNumFrames(); frameNr++) {
        vector<int> ids         = _peds_t[frameNr];
        vector<int> IdInFrame   = data.GetIdInFrame(frameNr, ids);
        vector<double> XInFrame = data.GetXInFrame(frameNr, ids);
        vector<double> YInFrame = data.GetYInFrame(frameNr, ids);
        for(unsigned int i = 0; i < IdInFrame.size(); i++) {
            point_2d p{XInFrame[i] * CMtoM, YInFrame[i] * CMtoM};
            bool isInBuilding = std::any_of(
                std::begin(_geometry), std::end(_geometry), [&p](const polygon_2d & polygon) {
                    return within(p, polygon);
                });
            if(!isInBuilding) {
                LOG_WARNING(
                    "Warning:\tAt {:d}th frame pedestrian at <x={}, y={}> is not in geometry!",
                    frameNr + data.GetMinFrame(),
                    XInFrame[i] * CMtoM,
                    YInFrame[i] * CMtoM);
            }
        }
    }
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------

    if(_DoesUseMethodA) // Method A
    {
        if(_areasForMethodA.empty()) {
            LOG_ERROR("Method A selected with no measurement area!");
            exit(EXIT_FAILURE);
        }
#pragma omp parallel for
        for(int i = 0; i < int(_areasForMethodA.size()); i++) {
            Method_A method_A;
            method_A.SetMeasurementArea(_areasForMethodA[i]);
            method_A.SetTimeInterval(_deltaT[i]);
            bool result_A = method_A.Process(data, _scriptsLocation, _areasForMethodA[i]->_zPos);
            if(result_A) {
                LOG_INFO(
                    "Success with Method A using measurement area id {}!\n",
                    _areasForMethodA[i]->_id);
            } else {
                LOG_ERROR(
                    "Failed with Method A using measurement area id {}!\n",
                    _areasForMethodA[i]->_id);
            }
        }
    }

    if(_DoesUseMethodB) // Method_B
    {
        if(_areasForMethodB.empty()) {
            LOG_ERROR("Method B selected with no measurement area!");
            exit(EXIT_FAILURE);
        }

#pragma omp parallel for
        for(int i = 0; i < int(_areasForMethodB.size()); i++) {
            Method_B method_B;
            method_B.SetMeasurementArea(_areasForMethodB[i]);
            bool result_B = method_B.Process(data);
            if(result_B) {
                LOG_INFO(
                    "Success with Method B using measurement area id {}!\n",
                    _areasForMethodB[i]->_id);
            } else {
                LOG_ERROR(
                    "Failed with Method B using measurement area id {}!\n",
                    _areasForMethodB[i]->_id);
            }
        }
    }

    if(_DoesUseMethodC) // Method C
    {
        if(_areasForMethodC.empty()) {
            LOG_ERROR("Method C selected with no measurement area!");
            exit(EXIT_FAILURE);
        }
#pragma omp parallel for
        for(int i = 0; i < int(_areasForMethodC.size()); i++) {
            Method_C method_C;
            method_C.SetMeasurementArea(_areasForMethodC[i]);
            bool result_C = method_C.Process(data, _areasForMethodC[i]->_zPos);
            if(result_C) {
                LOG_INFO(
                    "Success with Method C using measurement area id {}!\n",
                    _areasForMethodC[i]->_id);
            } else {
                LOG_ERROR(
                    "Failed with Method C using measurement area id {}!\n",
                    _areasForMethodC[i]->_id);
            }
        }
    }

    if(_DoesUseMethodD) // method_D
    {
        if(_areasForMethodD.empty()) {
            LOG_ERROR("Method D selected with no measurement area!");
            exit(EXIT_FAILURE);
        }

#pragma omp parallel for
        for(int i = 0; i < int(_areasForMethodD.size()); i++) {
            Method_D method_D;
            // TODO: setting and processing should be restructured. constructor?
            method_D.SetGeometryPolygon(_geoPolyMethodD[_areasForMethodD[i]->_id]);
            method_D.SetGeometryBoundaries(_lowVertexX, _lowVertexY, _highVertexX, _highVertexY);
            method_D.SetMeasurementArea(_areasForMethodD[i]);
            bool result_D;
            result_D = method_D.Process(configData_D, i, data, _areasForMethodD[i]->_zPos);
            if(result_D) {
                LOG_INFO(
                    "Success with Method D using measurement area id {}!\n",
                    _areasForMethodD[i]->_id);
            } else {
                LOG_ERROR(
                    "Failed with Method D using measurement area id {}!\n",
                    _areasForMethodD[i]->_id);
            }
        }
    }

    if(_DoesUseMethodE) // method_E
    {
        if(_areasForMethodE.empty()) {
            LOG_ERROR("Method E selected with no measurement area!");
            exit(EXIT_FAILURE);
        }
        for(int i = 0; i < int(_areasForMethodE.size()); i++) {
            Method_E method_E;
            method_E.SetMeasurementArea(_areasForMethodE[i]);
            method_E.SetLine(_linesForMethodE[i]);
            method_E.SetTimeInterval(_deltaTMethodE[i]);
            bool result_E = method_E.Process(data, _areasForMethodE[i]->_zPos);
            if(result_E) {
                LOG_INFO(
                    "Success with Method E using measurement area id {} and "
                    "line id {}!\n",
                    _areasForMethodE[i]->_id,
                    _linesForMethodE[i]->_id);
            } else {
                LOG_ERROR(
                    "Failed with Method E using measurement area id {} and "
                    "line id {}!\n",
                    _areasForMethodE[i]->_id,
                    _linesForMethodE[i]->_id);
            }
        }
    }

    if(_DoesUseMethodF) // method_F
    {
        if(_areasForMethodF.empty()) {
            LOG_ERROR("Method F selected with no measurement area!");
            exit(EXIT_FAILURE);
        }
        for(int i = 0; i < int(_areasForMethodF.size()); i++) {
            Method_F method_F;
            method_F.SetMeasurementArea(_areasForMethodF[i]);
            method_F.SetLine(_linesForMethodF[i]);
            method_F.SetTimeInterval(_deltaTMethodF[i]);
            bool result_F = method_F.Process(data, _areasForMethodF[i]->_zPos);
            if(result_F) {
                LOG_INFO(
                    "Success with Method F using measurement area id {} and "
                    "line id {}!\n",
                    _areasForMethodF[i]->_id,
                    _linesForMethodF[i]->_id);
            } else {
                LOG_ERROR(
                    "Failed with Method F using measurement area id {} and "
                    "line id {}!\n",
                    _areasForMethodF[i]->_id,
                    _linesForMethodF[i]->_id);
            }
        }
    }

    if(_DoesUseMethodG) // method_G
    {
        if(_areasForMethodG.empty()) {
            LOG_ERROR("Method G selected with no measurement area!");
            exit(EXIT_FAILURE);
        }
        for(int i = 0; i < int(_areasForMethodG.size()); i++) {
            Method_G method_G;
            method_G.SetMeasurementArea(_areasForMethodG[i]);
            method_G.SetTimeInterval(_deltaTMethodG[i]);
            method_G.SetDt(_dtMethodG[i]);
            method_G.SetNumberPolygons(_numberPolygonsMethodG[i]);
            method_G.SetPoints(_pointsMethodG[i]);
            bool result_G = method_G.Process(data);
            if(result_G) {
                LOG_INFO(
                    "Success with Method G using measurement area id {}!\n",
                    _areasForMethodG[i]->_id);
            } else {
                LOG_ERROR(
                    "Failed with Method G using measurement area id {}!\n",
                    _areasForMethodG[i]->_id);
            }
        }
    }

    if(_DoesUseMethodH) // method_H
    {
        if(_areasForMethodH.empty()) {
            LOG_ERROR("Method H selected with no measurement area!");
            exit(EXIT_FAILURE);
        }
        for(int i = 0; i < int(_areasForMethodH.size()); i++) {
            Method_H method_H;
            method_H.SetMeasurementArea(_areasForMethodH[i]);
            method_H.SetTimeInterval(_deltaTMethodH[i]);
            bool result_H = method_H.Process(data);
            if(result_H) {
                LOG_INFO(
                    "Success with Method H using measurement area id {}!\n",
                    _areasForMethodH[i]->_id);
            } else {
                LOG_ERROR(
                    "Failed with Method H using measurement area id {}!\n",
                    _areasForMethodH[i]->_id);
            }
        }
    }

    return 0;
}

FILE * Analysis::CreateFile(const string & filename)
{
    // create the directory for the file
    fs::path filepath = fs::path(filename.c_str()).parent_path();
    if(fs::is_directory(filepath) == false) {
        if(fs::create_directories(filepath) == false && fs::is_directory(filepath) == false) {
            LOG_ERROR("cannot create the directory <{}>", filepath.string());
            return NULL;
        }
        LOG_INFO("create the directory <{}>", filepath.string());
    }

    // create the file
    FILE * fHandle = fopen(filename.c_str(), "w");
    return fHandle;
}
