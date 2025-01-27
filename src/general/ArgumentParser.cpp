/**
 * \file        ArgumentParser.cpp
 * \date        Oct 10, 2014
 * \version     v0.8.3
 * \copyright   <2009-2018> Forschungszentrum Jülich GmbH. All rights reserved.
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
 * The ArgumentParser class define functions reading the input parameters from initial files.
 *
 *
 **/
#include "ArgumentParser.h"

#include "../Analysis.h"
#include "Compiler.h"
#include "Logger.h"
#include "tinyxml.h"

#include <boost/range/iterator_range.hpp>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <numeric>
#include <optional>
#include <sstream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif


using namespace std;

void Logs()
{
    time_t now = chrono::system_clock::to_time_t(chrono::system_clock::now());
    std::ostringstream oss;
    char foo[100];
    if(0 < std::strftime(foo, sizeof(foo), "%a %b %d %X %Y", std::localtime(&now)))
        oss << foo;
    else
        oss << "No time!";
    // else  // hack for g++ < 5
    //      oss << std::put_time(std::localtime(&now), "%a %b %d %X %Y");
    auto currentTime = oss.str();

    // first logs will go to stdout
    LOG_INFO("Starting JuPedSim - JPSreport");
    LOG_INFO("Version {}", JPSREPORT_VERSION);
    LOG_INFO("Commit id {}", GIT_COMMIT_HASH);
    LOG_INFO("Commit date {}", GIT_COMMIT_DATE);
    LOG_INFO("Build from branch {}", GIT_BRANCH);
    LOG_INFO("Build with {}({})", compiler_id, compiler_version);
}

void ArgumentParser::Usage(const std::string file)
{
    LOG_INFO("Usage: \n");
    LOG_INFO("{} inifile.xml\n", file);
    exit(EXIT_SUCCESS);
}

ArgumentParser::ArgumentParser()
{
    // Default parameter values
    _geometryFileName = "geo.xml";

    _vComponent             = "B";
    _IgnoreBackwardMovement = false;
    _isMethodA              = false;
    _delatTVInst            = 5;
    _isMethodB              = false;
    _isMethodC              = false;
    _isMethodD              = false;
    _isMethodE              = false;
    _isMethodF              = false;
    _isMethodG              = false;
    _isMethodH              = false;
    _steadyStart            = 100;
    _steadyEnd              = 1000;
    _trajectoriesLocation   = "./";
    _trajectoriesFilename   = "";
    _projectRootDir         = "./";
    _fileFormat             = FORMAT_XML_PLAIN;
}


bool ArgumentParser::ParseArgs(int argc, char ** argv)
{
    // special case of the default configuration ini.xml
    if(argc == 1) {
        LOG_INFO("Trying to load the default configuration from the file <ini.xml>");
        if(!ParseInputFiles("ini.xml")) {
            Usage(argv[0]);
        }
        return true;
    }

    string argument = argv[1];
    if(argument == "-h" || argument == "--help") {
        Usage(argv[0]);
    }
    if(argument == "-v" || argument == "--version") {
        Logs();
        exit(EXIT_SUCCESS);
    }

    // other special case where a single configuration file is submitted
    // check if inifile options are given
    if(argc == 2) {
        string prefix1 = "--ini=";
        string prefix2 = "--inifile=";

        if(!argument.compare(0, prefix2.size(), prefix2)) {
            argument.erase(0, prefix2.size());
        } else if(!argument.compare(0, prefix1.size(), prefix1)) {
            argument.erase(0, prefix1.size());
        }
        return ParseInputFiles(argument);
    }

    // more than one argument was supplied
    Usage(argv[0]);
    return false;
}

const vector<fs::path> & ArgumentParser::GetTrajectoriesFiles() const
{
    return _trajectoriesFiles;
}

const fs::path & ArgumentParser::GetProjectRootDir() const
{
    return _projectRootDir;
}

bool ArgumentParser::ParseInifile(const fs::path & inifile)
{
    Logs();
    LOG_INFO("Parsing the ini file <{}>", inifile);
    // extract and set the project root dir
    fs::path p(inifile);
    _projectRootDir = weakly_canonical(p).parent_path();
    TiXmlDocument doc(inifile.string());
    if(!doc.LoadFile()) {
        LOG_ERROR("{}", doc.ErrorDesc());
        LOG_ERROR("Could not parse the ini file");
        return false;
    }
    TiXmlElement * xMainNode = doc.RootElement();
    if(!xMainNode) {
        LOG_ERROR("Root element does not exist");
        return false;
    }

    if(xMainNode->ValueStr() != "JPSreport") {
        LOG_ERROR("Root element value is not 'JPSreport'.");
        return false;
    }

    // geometry
    if(xMainNode->FirstChild("geometry")) {
        fs::path pathGeo(xMainNode->FirstChildElement("geometry")->Attribute("file"));
        _geometryFileName = GetProjectRootDir() / pathGeo;
        if(!fs::exists(_geometryFileName)) {
            LOG_ERROR("Geometry File <{}> does not exist", _geometryFileName.string());
            return false;
        }
        _geometryFileName = fs::canonical(_geometryFileName);
        LOG_INFO("Geometry File is: <{}>", _geometryFileName.string());
    }

    // trajectories
    TiXmlNode * xTrajectories = xMainNode->FirstChild("trajectories");
    if(xTrajectories) {
        // add the extension point
        string fmt =
            "." + string(xmltoa(xMainNode->FirstChildElement("trajectories")->Attribute("format")));
        LOG_INFO("Format of the trajectory file is: <{}>", fmt);
        if(fmt == ".xml") {
            _fileFormat = FORMAT_XML_PLAIN;
        } else if(fmt == ".txt") {
            _fileFormat = FORMAT_PLAIN;
        } else {
            LOG_ERROR("the given trajectory format is not supported. Supply '.xml' or "
                      "'.txt' format!");
            return false;
        }

        string unit = xmltoa(xMainNode->FirstChildElement("trajectories")->Attribute("unit"), "m");
        if(unit != "m") {
            LOG_WARNING("only <m> unit is supported. Convert your units.");
            return false;
        }
        // a file descriptor was given
        for(TiXmlElement * xFile = xTrajectories->FirstChildElement("file"); xFile;
            xFile                = xFile->NextSiblingElement("file")) {
            // collect all the files given
            _trajectoriesFilename = fs::path(xFile->Attribute("name"));
            _trajectoriesFiles.push_back(_trajectoriesFilename);

            // check if the given file match the format
            if(boost::algorithm::ends_with(_trajectoriesFilename.string(), fmt)) {
                LOG_INFO("Input trajectory file is <{}>", _trajectoriesFilename.string());
            } else {
                LOG_ERROR(
                    "Wrong file extension\t<{}> for file <{}>",
                    fmt,
                    _trajectoriesFilename.string());
                return false;
            }
        }
        auto xmlpath = xTrajectories->FirstChildElement("path");
        if(xmlpath) {
            if(xmlpath->Attribute("location")) {
                _trajectoriesLocation =
                    GetProjectRootDir() / fs::path(xmlpath->Attribute("location"));
                // _trajectoriesLocation = canonical(_trajectoriesLocation);
            }
        } else {
            fs::path path_root(GetProjectRootDir());
            path_root             = canonical(path_root);
            _trajectoriesLocation = path_root.string();
        }
        LOG_INFO("Input directory for loading trajectory is <{}>", _trajectoriesLocation.string());

        // in the case no file was specified, collect all files in the specified directory
        if(_trajectoriesFiles.empty()) {
            if(exists(_trajectoriesLocation)) {
                /* print all the files and directories within directory */
                fs::path path_traj(GetTrajectoriesLocation());
                path_traj = canonical(path_traj);
                for(auto & filename :
                    boost::make_iterator_range(fs::directory_iterator(path_traj), {})) {
                    string s = filename.path().string();
                    int pos  = s.find_last_of('/');
                    pos      = pos == -1 ? int(s.find_last_of('\\')) : pos;
                    s        = pos == -1 ? s : s.substr(pos + 1);
                    if(boost::algorithm::ends_with(s, fmt)) {
                        _trajectoriesFiles.push_back(s);
                        LOG_INFO("Input trajectory file is <{}>", s);
                    }
                }
            } else {
                /* could not open directory */
                LOG_ERROR("could not open the directory <{}>", _trajectoriesLocation.string());
                return false;
            }
        }
    }

    // max CPU
    if(xMainNode->FirstChild("num_threads")) {
        TiXmlNode * numthreads = xMainNode->FirstChild("num_threads")->FirstChild();
        if(numthreads) {
#ifdef _OPENMP
            omp_set_num_threads(xmltoi(numthreads->Value(), omp_get_max_threads()));
#endif
        }
        LOG_INFO("Using <{}> threads", omp_get_max_threads());
    }

    // output directory
    _outputDir = GetProjectRootDir() / "Output";
    if(xMainNode->FirstChild("output")) {
        string tmp = xMainNode->FirstChildElement("output")->Attribute("location");
        fs::path tmpPath(tmp);
        _outputDir = tmpPath;
        if(tmp.empty()) {
            _outputDir = GetProjectRootDir() / "Output";
        }
        if(!_outputDir.is_absolute()) {
            _outputDir = _projectRootDir / _outputDir;
        }
    } else
        LOG_INFO("Default output directory");
    if(!exists(_outputDir)) {
        // does not exist yet. mkdir
        bool res = fs::create_directory(_outputDir);
        if(res == false) {
            LOG_ERROR("Could not create the directory <{}>", _outputDir.string());
            return false;
        } else
            LOG_INFO("created directory <{}>", _outputDir.string());
    }
    LOG_INFO("Output directory for results is: <{}>", _outputDir.string());

    // measurement area
    if(xMainNode->FirstChild("measurement_areas")) {
        string unit = "";
        if(xMainNode->FirstChildElement("measurement_areas")->Attribute("unit"))
            unit = xMainNode->FirstChildElement("measurement_areas")->Attribute("unit");
        if(unit != "m") {
            LOG_WARNING("only <m> unit is supported. Convert your units.");
            return false;
        }

        for(TiXmlNode * xMeasurementArea_B =
                xMainNode->FirstChild("measurement_areas")->FirstChild("area_B");
            xMeasurementArea_B;
            xMeasurementArea_B = xMeasurementArea_B->NextSibling("area_B")) {
            MeasurementArea_B * areaB = new MeasurementArea_B();
            areaB->_id                = xmltoi(xMeasurementArea_B->ToElement()->Attribute("id"));
            areaB->_type              = xMeasurementArea_B->ToElement()->Attribute("type");
            if(xMeasurementArea_B->ToElement()->Attribute("zPos")) {
                if(string(xMeasurementArea_B->ToElement()->Attribute("zPos")) != "None") {
                    areaB->_zPos = xmltof(xMeasurementArea_B->ToElement()->Attribute("zPos"));
                } else {
                    areaB->_zPos = 10000001.0;
                }
            } else {
                areaB->_zPos = 10000001.0;
            }

            polygon_2d poly;
            LOG_INFO("Measure area id  <{}> with type <{}>", areaB->_id, areaB->_type);
            int num_verteces = 0;
            for(TiXmlElement * xVertex = xMeasurementArea_B->FirstChildElement("vertex"); xVertex;
                xVertex                = xVertex->NextSiblingElement("vertex")) {
                // Note: Attributes are optional, their existence needs to be checked since xmltof
                // returns 0.0 if unknown
                if(xVertex->Attribute("x") != nullptr && xVertex->Attribute("y") != nullptr) {
                    // DEPRECATED FORMAT should be removed in future
                    double box_px = xmltof(xVertex->Attribute("x")) * M2CM;
                    double box_py = xmltof(xVertex->Attribute("y")) * M2CM;
                    boost::geometry::append(poly, boost::geometry::make<point_2d>(box_px, box_py));
                    LOG_INFO(
                        "Measure area points  <{:.3f}, {:.3f}>", box_px * CMtoM, box_py * CMtoM);
                    num_verteces++;
                } else if(
                    xVertex->Attribute("px") != nullptr && xVertex->Attribute("py") != nullptr) {
                    // NEW FORMAT
                    double box_px = xmltof(xVertex->Attribute("px")) * M2CM;
                    double box_py = xmltof(xVertex->Attribute("py")) * M2CM;
                    boost::geometry::append(poly, boost::geometry::make<point_2d>(box_px, box_py));
                    LOG_INFO(
                        "Measure area points  <{:.3f}, {:.3f}>", box_px * CMtoM, box_py * CMtoM);
                    num_verteces++;
                } else {
                    LOG_WARNING("Invalid vertex format given.");
                }
            }
            if(num_verteces < 3 && num_verteces > 0)
                LOG_WARNING(
                    "Less than 3 measure area points given ({}). At least 3 or nothing "
                    "at all!!",
                    num_verteces);

            correct(poly); // in the case the Polygone is not closed
            areaB->_poly = poly;

            TiXmlElement * xLength =
                xMeasurementArea_B->FirstChildElement("length_in_movement_direction");
            if(xLength) {
                areaB->_length = xmltof(xLength->Attribute("distance"));
                LOG_INFO("Length in movement direction {:.3f}", areaB->_length);
            }

            TiXmlElement * xLengthOrthogonal =
                xMeasurementArea_B->FirstChildElement("length_orthogonal_to_movement_direction");
            if(xLengthOrthogonal) {
                areaB->_lengthOrthogonal = xmltof(xLengthOrthogonal->Attribute("distance"));
                LOG_INFO(
                    "Length orthogonal to movement direction {:.3f}", areaB->_lengthOrthogonal);
            }
            // delta y for methods E and F (orthogonal to movement direction)

            _measurementAreasByIDs[areaB->_id] = areaB;
        }
        for(TiXmlNode * xMeasurementArea_L =
                xMainNode->FirstChild("measurement_areas")->FirstChild("area_L");
            xMeasurementArea_L;
            xMeasurementArea_L = xMeasurementArea_L->NextSibling("area_L")) {
            MeasurementArea_L * areaL = new MeasurementArea_L();
            areaL->_id                = xmltoi(xMeasurementArea_L->ToElement()->Attribute("id"));
            areaL->_type              = xMeasurementArea_L->ToElement()->Attribute("type");
            if(xMeasurementArea_L->ToElement()->Attribute("zPos")) {
                if(string(xMeasurementArea_L->ToElement()->Attribute("zPos")) != "None") {
                    areaL->_zPos = xmltof(xMeasurementArea_L->ToElement()->Attribute("zPos"));
                } else {
                    areaL->_zPos = 10000001.0;
                }
            } else {
                areaL->_zPos = 10000001.0;
            }
            LOG_INFO("Measurement area id  <{}> with type <{}>", areaL->_id, areaL->_type);

            if(xMeasurementArea_L->FirstChildElement("start")->Attribute("x") != nullptr &&
               xMeasurementArea_L->FirstChildElement("start")->Attribute("y") != nullptr) {
                areaL->_lineStartX =
                    xmltof(xMeasurementArea_L->FirstChildElement("start")->Attribute("x")) * M2CM;
                areaL->_lineStartY =
                    xmltof(xMeasurementArea_L->FirstChildElement("start")->Attribute("y")) * M2CM;
            } else if(
                xMeasurementArea_L->FirstChildElement("start")->Attribute("px") != nullptr &&
                xMeasurementArea_L->FirstChildElement("start")->Attribute("py") != nullptr) {
                // NEW FORMAT
                // Note: argument can be changed to "required" (once the deprecated format is
                // removed), existence would no longer need to be checked
                areaL->_lineStartX =
                    xmltof(xMeasurementArea_L->FirstChildElement("start")->Attribute("px")) * M2CM;
                areaL->_lineStartY =
                    xmltof(xMeasurementArea_L->FirstChildElement("start")->Attribute("py")) * M2CM;
            } else {
                LOG_ERROR("Invalid definition of measurement line start");
                exit(EXIT_FAILURE);
            }

            if(xMeasurementArea_L->FirstChildElement("end")->Attribute("x") != nullptr &&
               xMeasurementArea_L->FirstChildElement("end")->Attribute("y") != nullptr) {
                // DEPRECATED FORMAT should be removed in future
                areaL->_lineEndX =
                    xmltof(xMeasurementArea_L->FirstChildElement("end")->Attribute("x")) * M2CM;
                areaL->_lineEndY =
                    xmltof(xMeasurementArea_L->FirstChildElement("end")->Attribute("y")) * M2CM;
            } else if(
                xMeasurementArea_L->FirstChildElement("end")->Attribute("px") != nullptr &&
                xMeasurementArea_L->FirstChildElement("end")->Attribute("py") != nullptr) {
                // NEW FORMAT
                // Note: argument can be changed to "required" (once the deprecated format is
                // removed), existence would no longer need to be checked
                areaL->_lineEndX =
                    xmltof(xMeasurementArea_L->FirstChildElement("end")->Attribute("px")) * M2CM;
                areaL->_lineEndY =
                    xmltof(xMeasurementArea_L->FirstChildElement("end")->Attribute("py")) * M2CM;
            } else {
                LOG_ERROR("Invalid definition of measurement line end");
                exit(EXIT_FAILURE);
            }

            _measurementAreasByIDs[areaL->_id] = areaL;
            LOG_INFO(
                "Measurement line starts from  <{:.3f}, {:.3f}> to <{:.3f}, {:.3f}>",
                areaL->_lineStartX * CMtoM,
                areaL->_lineStartY * CMtoM,
                areaL->_lineEndX * CMtoM,
                areaL->_lineEndY * CMtoM);
        }
    }
    // instantaneous velocity
    TiXmlNode * xVelocity = xMainNode->FirstChild("velocity");
    if(xVelocity) {
        string FrameSteps = "10";
        if(xMainNode->FirstChildElement("velocity")->Attribute("frame_step")) {
            FrameSteps   = xMainNode->FirstChildElement("velocity")->Attribute("frame_step");
            _delatTVInst = atof(FrameSteps.c_str()) / 2.0;
        }
        string MovementDirection = "None";
        if(xMainNode->FirstChildElement("velocity")->Attribute("set_movement_direction")) {
            MovementDirection =
                xMainNode->FirstChildElement("velocity")->Attribute("set_movement_direction");
            if(atof(MovementDirection.c_str()) < 0 && atof(MovementDirection.c_str()) > 360 &&
               MovementDirection != "None" && MovementDirection != "SeeTraj") {
                LOG_WARNING("The movement direction should be set between 0 to 360 or None!");
                return false;
            }
        }
        if(xMainNode->FirstChildElement("velocity")->Attribute("ignore_backward_movement")) {
            if(string(xMainNode->FirstChildElement("velocity")
                          ->Attribute("ignore_backward_movement")) == "true") {
                _IgnoreBackwardMovement = true;
            } else {
                _IgnoreBackwardMovement = false;
            }
        }
        if(MovementDirection == "None") {
            _vComponent             = "B"; // both components
            _IgnoreBackwardMovement = false;
            LOG_INFO(
                "Both x and y-component of coordinates will be used to calculate instantaneous "
                "velocity over <{}> frames",
                FrameSteps);
        } else if(MovementDirection == "SeeTraj") {
            _vComponent = "F";
            LOG_INFO(
                "The component defined in the trajectory file will be used to calculate "
                "instantaneous velocity over <{}> frames",
                FrameSteps);
        } else {
            _vComponent = MovementDirection;
            LOG_INFO(
                "The instantaneous velocity in the direction of <{}>"
                " will be calculated over <{}> frames",
                MovementDirection,
                FrameSteps);
        }
    }
    // Method A
    TiXmlElement * xMethod_A = xMainNode->FirstChildElement("method_A");
    if(xMethod_A) {
        if(string(xMethod_A->Attribute("enabled")) == "true") {
            _isMethodA = true;
            LOG_INFO("Method A is selected");
            /*               _timeIntervalA =
               xmltoi(xMethod_A->FirstChildElement("frame_interval")->GetText());
                 Log->Write("INFO: \tFrame interval used for calculating flow in Method A is <%d>
               frame",_timeIntervalA);*/
            for(TiXmlElement * xMeasurementArea =
                    xMainNode->FirstChildElement("method_A")->FirstChildElement("measurement_area");
                xMeasurementArea;
                xMeasurementArea = xMeasurementArea->NextSiblingElement("measurement_area")) {
                int id = xmltoi(xMeasurementArea->Attribute("id"));

                if(_measurementAreasByIDs[id]->_type == "Line") {
                    _areaIDforMethodA.push_back(id);
                    LOG_INFO("Measurement area id <{}> will be used for analysis", id);
                } else {
                    LOG_WARNING(
                        "Measurement area id <{}> will NOT be used for analysis (Type "
                        "<{}> is not Line)",
                        id,
                        _measurementAreasByIDs[id]->_type);
                }

                if(xMeasurementArea->Attribute("frame_interval")) {
                    if(string(xMeasurementArea->Attribute("frame_interval")) != "None") {
                        _timeIntervalA.push_back(
                            xmltoi(xMeasurementArea->Attribute("frame_interval")));
                        LOG_INFO(
                            "Frame interval used for calculating flow is <{}> frame",
                            xmltoi(xMeasurementArea->Attribute("frame_interval")));
                    } else {
                        _timeIntervalA.push_back(100);
                    }
                } else {
                    _timeIntervalA.push_back(100);
                }
            }
        }
    }
    // method B
    TiXmlElement * xMethod_B = xMainNode->FirstChildElement("method_B");
    if(xMethod_B)

        if(string(xMethod_B->Attribute("enabled")) == "true") {
            _isMethodB = true;
            LOG_INFO("Method B is selected");
            for(TiXmlElement * xMeasurementArea =
                    xMainNode->FirstChildElement("method_B")->FirstChildElement("measurement_area");
                xMeasurementArea;
                xMeasurementArea = xMeasurementArea->NextSiblingElement("measurement_area")) {
                _areaIDforMethodB.push_back(xmltoi(xMeasurementArea->Attribute("id")));
                LOG_INFO(
                    "Measurement area id <{}> will be used for analysis",
                    xmltoi(xMeasurementArea->Attribute("id")));
            }
        }
    // method C
    TiXmlElement * xMethod_C = xMainNode->FirstChildElement("method_C");
    if(xMethod_C) {
        if(string(xMethod_C->Attribute("enabled")) == "true") {
            _isMethodC = true;
            LOG_INFO("Method C is selected");
            for(TiXmlElement * xMeasurementArea =
                    xMainNode->FirstChildElement("method_C")->FirstChildElement("measurement_area");
                xMeasurementArea;
                xMeasurementArea = xMeasurementArea->NextSiblingElement("measurement_area")) {
                _areaIDforMethodC.push_back(xmltoi(xMeasurementArea->Attribute("id")));
                LOG_INFO(
                    "Measurement area id <{}> will be used for analysis",
                    xmltoi(xMeasurementArea->Attribute("id")));
            }
        }
    }

    // method D
    TiXmlElement * xMethod_D = xMainNode->FirstChildElement("method_D");
    if(xMethod_D) {
        LOG_INFO("Method D is selected with following options");
        if(auto configData = ParseDIJParams(xMethod_D)) {
            _configDataD = configData.value();
            _isMethodD   = true;
        }
    }

    // method E
    TiXmlElement * xMethod_E = xMainNode->FirstChildElement("method_E");
    if(xMethod_E) {
        if(string(xMethod_E->Attribute("enabled")) == "true") {
            _isMethodE = true;
            LOG_INFO("Method E is selected");
            for(TiXmlElement * xMeasurementArea =
                    xMainNode->FirstChildElement("method_E")->FirstChildElement("measurement_area");
                xMeasurementArea;
                xMeasurementArea = xMeasurementArea->NextSiblingElement("measurement_area")) {
                int id      = xmltoi(xMeasurementArea->Attribute("id"));
                int line_id = xmltoi(xMeasurementArea->Attribute("line_id"));

                if(_measurementAreasByIDs[id]->_type == "BoundingBox" &&
                   _measurementAreasByIDs[line_id]->_type == "Line") {
                    if(IsInMeasureArea(
                           dynamic_cast<MeasurementArea_L *>(GetMeasurementArea(line_id)),
                           dynamic_cast<MeasurementArea_B *>(GetMeasurementArea(id)))) {
                        _areaIDforMethodE.push_back(id);
                        _lineIDforMethodE.push_back(line_id);
                        LOG_INFO("Measurement area id <{}> will be used for analysis", id);
                        if(xMeasurementArea->Attribute("frame_interval")) {
                            if(string(xMeasurementArea->Attribute("frame_interval")) != "None") {
                                _timeIntervalE.push_back(
                                    xmltoi(xMeasurementArea->Attribute("frame_interval")));
                                LOG_INFO(
                                    "Frame interval used for calculating density is <{}> frames",
                                    xmltoi(xMeasurementArea->Attribute("frame_interval")));
                            } else {
                                _timeIntervalE.push_back(-1);
                            }
                        } else {
                            _timeIntervalE.push_back(-1);
                        }
                    } else {
                        LOG_WARNING(
                            "Measurement area id <{}> with line id <{}> will NOT be used for "
                            "analysis: The line is not located within the measurement area.",
                            id,
                            line_id);
                    }
                } else {
                    LOG_WARNING(
                        "Measurement area id <{}> will NOT be used for analysis: Either type "
                        "of measurement area ({}) is not BoundingBox, or type of line ({}) "
                        "is not Line.",
                        id,
                        _measurementAreasByIDs[id]->_type,
                        _measurementAreasByIDs[line_id]->_type);
                }
            }
        }
    }

    // method F
    TiXmlElement * xMethod_F = xMainNode->FirstChildElement("method_F");
    if(xMethod_F) {
        if(string(xMethod_F->Attribute("enabled")) == "true") {
            _isMethodF = true;
            LOG_INFO("Method F is selected");
            for(TiXmlElement * xMeasurementArea =
                    xMainNode->FirstChildElement("method_F")->FirstChildElement("measurement_area");
                xMeasurementArea;
                xMeasurementArea = xMeasurementArea->NextSiblingElement("measurement_area")) {
                int id      = xmltoi(xMeasurementArea->Attribute("id"));
                int line_id = xmltoi(xMeasurementArea->Attribute("line_id"));

                if(_measurementAreasByIDs[id]->_type == "BoundingBox" &&
                   _measurementAreasByIDs[line_id]->_type == "Line") {
                    if(IsInMeasureArea(
                           dynamic_cast<MeasurementArea_L *>(GetMeasurementArea(line_id)),
                           dynamic_cast<MeasurementArea_B *>(GetMeasurementArea(id)))) {
                        _areaIDforMethodF.push_back(id);
                        _lineIDforMethodF.push_back(line_id);
                        LOG_INFO("Measurement area id <{}> will be used for analysis", id);
                        if(xMeasurementArea->Attribute("frame_interval")) {
                            if(string(xMeasurementArea->Attribute("frame_interval")) != "None") {
                                _timeIntervalF.push_back(
                                    xmltoi(xMeasurementArea->Attribute("frame_interval")));
                                LOG_INFO(
                                    "Frame interval used for calculating density is <{}> frames",
                                    xmltoi(xMeasurementArea->Attribute("frame_interval")));
                            } else {
                                _timeIntervalF.push_back(-1);
                            }
                        } else {
                            _timeIntervalF.push_back(-1);
                        }
                    } else {
                        LOG_WARNING(
                            "Measurement area id <{}> with line id <{}> will NOT be used for "
                            "analysis: The line is not located within the measurement area.",
                            id,
                            line_id);
                    }
                } else {
                    LOG_WARNING(
                        "Measurement area id <{}> will NOT be used for analysis: Either type "
                        "of measurement area ({}) is not BoundingBox, or type of line ({}) "
                        "is not Line.",
                        id,
                        _measurementAreasByIDs[id]->_type,
                        _measurementAreasByIDs[line_id]->_type);
                }
            }
        }
    }

    // method G
    TiXmlElement * xMethod_G = xMainNode->FirstChildElement("method_G");
    if(xMethod_G) {
        if(string(xMethod_G->Attribute("enabled")) == "true") {
            _isMethodG = true;
            LOG_INFO("Method G is selected");
            for(TiXmlElement * xMeasurementArea =
                    xMainNode->FirstChildElement("method_G")->FirstChildElement("measurement_area");
                xMeasurementArea;
                xMeasurementArea = xMeasurementArea->NextSiblingElement("measurement_area")) {
                int id = xmltoi(xMeasurementArea->Attribute("id"));

                if(_measurementAreasByIDs[id]->_type == "BoundingBox") {
                    TiXmlElement * xPoint1 = xMainNode->FirstChildElement("method_G")
                                                 ->FirstChildElement("measurement_area")
                                                 ->FirstChildElement("point_1");
                    TiXmlElement * xPoint2 = xMainNode->FirstChildElement("method_G")
                                                 ->FirstChildElement("measurement_area")
                                                 ->FirstChildElement("point_2");
                    TiXmlElement * numberPolygons = xMainNode->FirstChildElement("method_G")
                                                        ->FirstChildElement("measurement_area")
                                                        ->FirstChildElement("number_areas");

                    if(xPoint1->Attribute("x") && xPoint1->Attribute("y") &&
                       xPoint2->Attribute("x") && xPoint2->Attribute("y")) {
                        _areaIDforMethodG.push_back(id);
                        LOG_INFO("Measurement area id <{}> will be used for analysis", id);
                        if(xMeasurementArea->Attribute("frame_interval")) {
                            if(string(xMeasurementArea->Attribute("frame_interval")) != "None") {
                                _timeIntervalG.push_back(
                                    xmltoi(xMeasurementArea->Attribute("frame_interval")));
                                LOG_INFO(
                                    "Frame interval used for calculation is <{}> frames",
                                    xmltoi(xMeasurementArea->Attribute("frame_interval")));
                            } else {
                                _timeIntervalG.push_back(-1);
                            }
                        } else {
                            _timeIntervalG.push_back(-1);
                        }
                        if(xMeasurementArea->Attribute("dt")) {
                            if(string(xMeasurementArea->Attribute("dt")) != "None") {
                                _dtMethodG.push_back(xmltoi(xMeasurementArea->Attribute("dt")));
                                LOG_INFO(
                                    "Small frame interval (dt) used for calculation is <{}> frames",
                                    xmltoi(xMeasurementArea->Attribute("dt")));
                            } else {
                                _dtMethodG.push_back(4); // what is a good default value?
                            }
                        } else {
                            _dtMethodG.push_back(4); // what is a good default value?
                        }
                        if(numberPolygons->Attribute("n")) {
                            _numberPolygonsMethodG.push_back(
                                xmltoi(numberPolygons->Attribute("n")));
                        } else {
                            _numberPolygonsMethodG.push_back(10); // what is a good default value?
                        }

                        vector<point_2d> tmpPoints;
                        double x = xmltof(xPoint1->Attribute("x")) * M2CM;
                        double y = xmltof(xPoint1->Attribute("y")) * M2CM;
                        tmpPoints.push_back(boost::geometry::make<point_2d>(x, y));
                        x = xmltof(xPoint2->Attribute("x")) * M2CM;
                        y = xmltof(xPoint2->Attribute("y")) * M2CM;
                        tmpPoints.push_back(boost::geometry::make<point_2d>(x, y));
                        _pointsMethodG.push_back(tmpPoints);
                    } else {
                        LOG_WARNING(
                            "Measurement area id <{}> will NOT be used for analysis "
                            "(no side of measurement area was given)",
                            id,
                            _measurementAreasByIDs[id]->_type);
                    }
                } else {
                    LOG_WARNING(
                        "Measurement area id <{}> will NOT be used for analysis (Type "
                        "<{}> is not BoundingBox)",
                        id,
                        _measurementAreasByIDs[id]->_type);
                }
            }
        }
    }

    // method H
    TiXmlElement * xMethod_H = xMainNode->FirstChildElement("method_H");
    if(xMethod_H) {
        if(string(xMethod_H->Attribute("enabled")) == "true") {
            _isMethodH = true;
            LOG_INFO("Method H is selected");
            for(TiXmlElement * xMeasurementArea =
                    xMainNode->FirstChildElement("method_H")->FirstChildElement("measurement_area");
                xMeasurementArea;
                xMeasurementArea = xMeasurementArea->NextSiblingElement("measurement_area")) {
                int id = xmltoi(xMeasurementArea->Attribute("id"));

                if(_measurementAreasByIDs[id]->_type == "BoundingBox") {
                    _areaIDforMethodH.push_back(id);
                    LOG_INFO("Measurement area id <{}> will be used for analysis", id);
                    if(xMeasurementArea->Attribute("frame_interval")) {
                        if(string(xMeasurementArea->Attribute("frame_interval")) != "None") {
                            _timeIntervalH.push_back(
                                xmltoi(xMeasurementArea->Attribute("frame_interval")));
                            LOG_INFO(
                                "Frame interval for calculation is <{}> frames",
                                xmltoi(xMeasurementArea->Attribute("frame_interval")));
                        } else {
                            _timeIntervalH.push_back(-1);
                        }
                    } else {
                        _timeIntervalH.push_back(-1);
                    }
                } else {
                    LOG_WARNING(
                        "Measurement area id <{}> will NOT be used for analysis (Type "
                        "<{}> is not BoundingBox)",
                        id,
                        _measurementAreasByIDs[id]->_type);
                }
            }
        }
    }

    LOG_INFO("Finish parsing inifile");
    if(!(_isMethodA || _isMethodB || _isMethodC || _isMethodD || _isMethodE || _isMethodF ||
         _isMethodG || _isMethodH)) {
        LOG_WARNING("No measurement method enabled. Nothing to do.");
        exit(EXIT_SUCCESS);
    }
    return true;
}

bool ArgumentParser::ParseInputFiles(const string & inifile)
{
    if(ParseInifile(inifile)) {
        if(auto geometry = ParseGeometry(_geometryFileName); geometry.has_value()) {
            _geometry = geometry.value();
            return true;
        }
        return false;
    }

    return false;
}

std::optional<std::vector<polygon_2d>> ArgumentParser::ParseGeometry(const fs::path & geometryFile)
{
    LOG_INFO("ReadGeometry with {}.", geometryFile.string());

    std::vector<polygon_2d> geoPoly;

    TiXmlDocument docGeo(geometryFile.string());
    if(!docGeo.LoadFile()) {
        LOG_ERROR("{}", docGeo.ErrorDesc());
        LOG_ERROR("could not parse the geometry file");
        return std::nullopt;
    }

    TiXmlElement * xRootNode = docGeo.RootElement();
    if(!xRootNode) {
        LOG_ERROR("Root element does not exist");
        return std::nullopt;
    }

    if(xRootNode->ValueStr() != "geometry") {
        LOG_ERROR("Root element value is not 'geometry'.");
        return std::nullopt;
    }
    if(xRootNode->Attribute("unit"))
        if(string(xRootNode->Attribute("unit")) != "m") {
            LOG_ERROR(
                "Only the unit m (meters) is supported. \n\tYou supplied [{}]",
                xRootNode->Attribute("unit"));
            return std::nullopt;
        }

    double version = xmltof(xRootNode->Attribute("version"), -1);

    if(version != std::stod(JPS_VERSION) && version != std::stod(JPS_OLD_VERSION)) {
        LOG_ERROR(" Wrong geometry version!");
        LOG_ERROR(" Only version >= {} supported", JPS_VERSION);
        LOG_ERROR(" Please update the version of your geometry file to {}", JPS_VERSION);
        return std::nullopt;
    }

    // processing the rooms node
    TiXmlNode * xRoomsNode = xRootNode->FirstChild("rooms");
    if(!xRoomsNode) {
        LOG_ERROR("The geometry should have at least one room and one subroom");
        return std::nullopt;
    }

    for(TiXmlElement * xRoom = xRoomsNode->FirstChildElement("room"); xRoom;
        xRoom                = xRoom->NextSiblingElement("room")) {
        string room_id = xmltoa(xRoom->Attribute("id"), "-1");

        boost::geometry::model::multi_polygon<polygon_2d> room;
        // parsing the subrooms
        // processing the rooms node
        for(TiXmlElement * xSubRoom = xRoom->FirstChildElement("subroom"); xSubRoom;
            xSubRoom                = xSubRoom->NextSiblingElement("subroom")) {
            polygon_2d subroom;

            string subroom_id = xmltoa(xSubRoom->Attribute("id"), "-1");

            // looking for polygons (walls)
            for(TiXmlElement * xPolyVertices = xSubRoom->FirstChildElement("polygon");
                xPolyVertices;
                xPolyVertices = xPolyVertices->NextSiblingElement("polygon")) {
                for(TiXmlElement * xVertex = xPolyVertices->FirstChildElement("vertex");
                    xVertex && xVertex != xPolyVertices->LastChild("vertex");
                    xVertex = xVertex->NextSiblingElement("vertex")) {
                    double x1 = xmltof(xVertex->Attribute("px"));
                    double y1 = xmltof(xVertex->Attribute("py"));
                    double x2 = xmltof(xVertex->NextSiblingElement("vertex")->Attribute("px"));
                    double y2 = xmltof(xVertex->NextSiblingElement("vertex")->Attribute("py"));
                    append(subroom, make<point_2d>(x1 * M2CM, y1 * M2CM));
                    append(subroom, make<point_2d>(x2 * M2CM, y2 * M2CM));
                }
            }
            correct(subroom);

            std::vector<polygon_2d> obstacles;
            // looking for obstacles
            for(TiXmlElement * xObstacle = xSubRoom->FirstChildElement("obstacle"); xObstacle;
                xObstacle                = xObstacle->NextSiblingElement("obstacle")) {
                // looking for polygons (walls)
                polygon_2d obstacle;

                for(TiXmlElement * xPolyVertices = xObstacle->FirstChildElement("polygon");
                    xPolyVertices;
                    xPolyVertices = xPolyVertices->NextSiblingElement("polygon")) {
                    for(TiXmlElement * xVertex = xPolyVertices->FirstChildElement("vertex");
                        xVertex && xVertex != xPolyVertices->LastChild("vertex");
                        xVertex = xVertex->NextSiblingElement("vertex")) {
                        double x1 = xmltof(xVertex->Attribute("px"));
                        double y1 = xmltof(xVertex->Attribute("py"));
                        double x2 = xmltof(xVertex->NextSiblingElement("vertex")->Attribute("px"));
                        double y2 = xmltof(xVertex->NextSiblingElement("vertex")->Attribute("py"));

                        append(obstacle, make<point_2d>(x1 * M2CM, y1 * M2CM));
                        append(obstacle, make<point_2d>(x2 * M2CM, y2 * M2CM));
                    }
                }
                correct(obstacle);
            }

            // Add obstacles as holes in subroom polygon
            int k = 1;
            for(auto && obstacle : obstacles) {
                subroom.inners().resize(k++);
                subroom.inners().back();
                model::ring<point_2d> & inner = subroom.inners().back();
                for(auto && tmp_point : obstacle.outer()) {
                    append(inner, make<point_2d>(tmp_point.x() * M2CM, tmp_point.y() * M2CM));
                }
                correct(subroom);
            }
            unique(subroom);
            geoPoly.emplace_back(subroom);
        }
    }
    return std::optional<std::vector<polygon_2d>>{geoPoly};
}

std::optional<ConfigData_D> ArgumentParser::ParseDIJParams(TiXmlElement * xMethod)
{
    if(string(xMethod->Attribute("enabled")) == "false") {
        LOG_INFO("Method is disabled");
        return nullopt;
    }

    ConfigData_D configData;
    for(TiXmlElement * xMeasurementArea = xMethod->FirstChildElement("measurement_area");
        xMeasurementArea;
        xMeasurementArea = xMeasurementArea->NextSiblingElement("measurement_area")) {
        configData.areaIDs.push_back(xmltoi(xMeasurementArea->Attribute("id")));
        LOG_INFO(
            "Measurement area id <{}> will be used for analysis",
            xmltoi(xMeasurementArea->Attribute("id")));

        if(xMeasurementArea->Attribute("start_frame") &&
           string(xMeasurementArea->Attribute("start_frame")) != "None") {
            configData.startFrames.push_back(xmltoi(xMeasurementArea->Attribute("start_frame")));
            LOG_INFO(
                "the analysis starts from frame <{}>",
                xmltoi(xMeasurementArea->Attribute("start_frame")));
        } else {
            configData.startFrames.push_back(-1);
        }
        if(xMeasurementArea->Attribute("stop_frame") &&
           string(xMeasurementArea->Attribute("stop_frame")) != "None") {
            configData.stopFrames.push_back(xmltoi(xMeasurementArea->Attribute("stop_frame")));
            LOG_INFO(
                "the analysis stops from frame <{}>",
                xmltoi(xMeasurementArea->Attribute("stop_frame")));
        } else {
            configData.stopFrames.push_back(-1);
        }

        if(xMeasurementArea->Attribute("local_IFD") &&
           string(xMeasurementArea->Attribute("local_IFD")) == "true") {
            configData.calcLocalIFD.push_back(true);
            LOG_INFO("Local individual FD will be output");
        } else {
            configData.calcLocalIFD.push_back(false);
        }
    }
    if(xMethod->FirstChildElement("one_dimensional") &&
       string(xMethod->FirstChildElement("one_dimensional")->Attribute("enabled")) == "true") {
        configData.isOneDimensional = true;
        LOG_INFO("The data will be analyzed with one dimensional way!!");
    }

    if(xMethod->FirstChildElement("cut_by_circle") &&
       string(xMethod->FirstChildElement("cut_by_circle")->Attribute("enabled")) == "true") {
        configData.cutByCircle = true;
        configData.cutRadius =
            xmltof(xMethod->FirstChildElement("cut_by_circle")->Attribute("radius")) * M2CM;
        configData.circleEdges =
            xmltoi(xMethod->FirstChildElement("cut_by_circle")->Attribute("edges"));
        LOG_INFO(
            "Each Voronoi cell will be cut by a circle with the radius of <{}> m",
            configData.cutRadius * CMtoM);
        LOG_INFO(
            "The circle is discretized to a polygon with <{}> edges!!", configData.circleEdges);
    }

    if(xMethod->FirstChildElement("steadyState")) {
        _steadyStart = xmltof(xMethod->FirstChildElement("steadyState")->Attribute("start"));
        _steadyEnd   = xmltof(xMethod->FirstChildElement("steadyState")->Attribute("end"));
        LOG_INFO("the steady state is from <{}> to <{}> frames", _steadyStart, _steadyEnd);
    }

    if(xMethod->FirstChildElement("profiles") &&
       string(xMethod->FirstChildElement("profiles")->Attribute("enabled")) == "true") {
        TiXmlElement * xProfile = xMethod->FirstChildElement("profiles");
        configData.getProfile   = true;
        configData.gridSizeX    = xmltof(xProfile->Attribute("grid_size_x")) * M2CM;
        configData.gridSizeY    = xmltof(xProfile->Attribute("grid_size_y")) * M2CM;
        LOG_INFO("Profiles will be calculated");
        LOG_INFO(
            "The discretized grid size in x, y direction is: <{}> by <{}> m^2",
            configData.gridSizeX * CMtoM,
            configData.gridSizeY * CMtoM);

        // read in start and stop frame
        if(xProfile->Attribute("start_frame") &&
           string(xProfile->Attribute("start_frame")) != "None") {
            configData.startFrames.push_back(xmltoi(xProfile->Attribute("start_frame")));
            LOG_INFO(
                "the profile analysis starts from frame <{}>",
                xmltoi(xProfile->Attribute("start_frame")));
        } else {
            configData.startFrames.push_back(-1);
        }
        if(xProfile->Attribute("stop_frame") &&
           string(xProfile->Attribute("stop_frame")) != "None") {
            configData.stopFrames.push_back(xmltoi(xProfile->Attribute("stop_frame")));
            LOG_INFO(
                "the profile analysis stops from frame <{}>",
                xmltoi(xProfile->Attribute("stop_frame")));
        } else {
            configData.stopFrames.push_back(-1);
        }

        // TODO: restructuring needed. profiles option should be independet of MA so that start and
        // stop frame can be set here.
        // create MA with polygon points outside the geometry
        MeasurementArea_B * areaB = new MeasurementArea_B();
        areaB->_id                = -2;
        areaB->_type              = "Bounding Box";
        polygon_2d poly;
        areaB->_poly                       = poly;
        areaB->_zPos                       = 10000001.0; // TODO: why do we need to do this??
        _measurementAreasByIDs[areaB->_id] = areaB;

        // set config parameters
        configData.areaIDs.push_back(areaB->_id);
        configData.calcLocalIFD.push_back(false);
    }

    if(xMethod->FirstChildElement("use_blind_points") &&
       string(xMethod->FirstChildElement("use_blind_points")->Attribute("enabled")) == "false") {
        configData.useBlindPoints = false;
        LOG_INFO("Use of blind points disabled");
    }

    if(xMethod->FirstChildElement("vel_calculation") &&
       string(xMethod->FirstChildElement("vel_calculation")->Attribute("type")) == "Arithmetic") {
        /** Arithmetic velocity calculation is chosen.
         * Calculates the velocity of pedestrians based as arithmetic mean of their instantaneous
         *velocity. Method is independent of size of Voronoi cells.
         **/
        configData.velocityCalcFunc = [](const polygon_list & polygons,
                                         const vector<double> & individualVelocity,
                                         const polygon_2d &
                                         /*measurementArea*/) -> double {
            double arithmeticVelocity = 0;
            int pedsInMeasurementArea = polygons.size();
            double velocitySum =
                std::accumulate(individualVelocity.begin(), individualVelocity.end(), 0.0);

            if(pedsInMeasurementArea != 0) {
                arithmeticVelocity = velocitySum / (1.0 * pedsInMeasurementArea);
            }
            return arithmeticVelocity;
        };
        configData.velocityType = "Arithmetic";
        LOG_INFO("Arithmetic velocity calculation is used.");
    } else {
        LOG_INFO("Default Voronoi velocity calculation is used.");
    }

    if(xMethod->FirstChildElement("global_IFD") &&
       string(xMethod->FirstChildElement("global_IFD")->Attribute("enabled")) == "true") {
        auto xGlobal = xMethod->FirstChildElement("global_IFD");
        LOG_INFO(
            "Global IFD data will be calculated. Bounding box is created as measurement area.");

        // create MA with polygon points outside the geometry
        MeasurementArea_B * areaB = new MeasurementArea_B();
        areaB->_id                = -1;
        areaB->_type              = "Bounding Box";
        polygon_2d poly;
        areaB->_poly                       = poly;
        areaB->_zPos                       = 10000001.0; // TODO: why do we need to do this??
        _measurementAreasByIDs[areaB->_id] = areaB;

        // set config parameters
        configData.areaIDs.push_back(areaB->_id);
        configData.calcLocalIFD.push_back(false);

        // read in start and stop frame
        if(xGlobal->Attribute("start_frame") &&
           string(xGlobal->Attribute("start_frame")) != "None") {
            configData.startFrames.push_back(xmltoi(xGlobal->Attribute("start_frame")));
            LOG_INFO(
                "the global IFD analysis starts from frame <{}>",
                xmltoi(xGlobal->Attribute("start_frame")));
        } else {
            configData.startFrames.push_back(-1);
        }
        if(xGlobal->Attribute("stop_frame") && string(xGlobal->Attribute("stop_frame")) != "None") {
            configData.stopFrames.push_back(xmltoi(xGlobal->Attribute("stop_frame")));
            LOG_INFO(
                "the global IFD analysis stops from frame <{}>",
                xmltoi(xGlobal->Attribute("stop_frame")));
        } else {
            configData.stopFrames.push_back(-1);
        }
    }

    return configData;
}

const fs::path & ArgumentParser::GetGeometryFilename() const
{
    return _geometryFileName;
}

const FileFormat & ArgumentParser::GetFileFormat() const
{
    return _fileFormat;
}

const fs::path & ArgumentParser::GetTrajectoriesLocation() const
{
    return _trajectoriesLocation;
}

const fs::path & ArgumentParser::GetOutputLocation() const
{
    return _outputDir;
}

const fs::path & ArgumentParser::GetTrajectoriesFilename() const
{
    return _trajectoriesFilename;
}

std::string ArgumentParser::GetVComponent() const
{
    return _vComponent;
}

bool ArgumentParser::GetIgnoreBackwardMovement() const
{
    return _IgnoreBackwardMovement;
}

int ArgumentParser::GetDelatT_Vins() const
{
    return _delatTVInst;
}


bool ArgumentParser::GetIsMethodA() const
{
    return _isMethodA;
}

vector<int> ArgumentParser::GetTimeIntervalA() const
{
    return _timeIntervalA;
}

vector<int> ArgumentParser::GetTimeIntervalE() const
{
    return _timeIntervalE;
}

vector<int> ArgumentParser::GetTimeIntervalF() const
{
    return _timeIntervalF;
}

vector<int> ArgumentParser::GetTimeIntervalG() const
{
    return _timeIntervalG;
}

vector<int> ArgumentParser::GetTimeIntervalH() const
{
    return _timeIntervalH;
}

vector<int> ArgumentParser::GetDtMethodG() const
{
    return _dtMethodG;
}

bool ArgumentParser::GetIsMethodB() const
{
    return _isMethodB;
}

bool ArgumentParser::GetIsMethodC() const
{
    return _isMethodC;
}

bool ArgumentParser::GetIsMethodD() const
{
    return _isMethodD;
}

bool ArgumentParser::GetIsMethodE() const
{
    return _isMethodE;
}

bool ArgumentParser::GetIsMethodF() const
{
    return _isMethodF;
}

bool ArgumentParser::GetIsMethodG() const
{
    return _isMethodG;
}

bool ArgumentParser::GetIsMethodH() const
{
    return _isMethodH;
}

double ArgumentParser::GetSteadyStart() const
{
    return _steadyStart;
}

double ArgumentParser::GetSteadyEnd() const
{
    return _steadyEnd;
}

vector<int> ArgumentParser::GetAreaIDforMethodA() const
{
    return _areaIDforMethodA;
}

vector<int> ArgumentParser::GetAreaIDforMethodB() const
{
    return _areaIDforMethodB;
}

vector<int> ArgumentParser::GetAreaIDforMethodC() const
{
    return _areaIDforMethodC;
}

vector<int> ArgumentParser::GetAreaIDforMethodE() const
{
    return _areaIDforMethodE;
}

vector<int> ArgumentParser::GetAreaIDforMethodF() const
{
    return _areaIDforMethodF;
}

vector<int> ArgumentParser::GetAreaIDforMethodG() const
{
    return _areaIDforMethodG;
}

vector<int> ArgumentParser::GetAreaIDforMethodH() const
{
    return _areaIDforMethodH;
}

vector<int> ArgumentParser::GetLineIDforMethodE() const
{
    return _lineIDforMethodE;
}

vector<int> ArgumentParser::GetLineIDforMethodF() const
{
    return _lineIDforMethodF;
}

std::vector<vector<point_2d>> ArgumentParser::GetPointsMethodG() const
{
    return _pointsMethodG;
}

std::vector<int> ArgumentParser::GetNumPolyMethodG() const
{
    return _numberPolygonsMethodG;
}

MeasurementArea * ArgumentParser::GetMeasurementArea(int id)
{
    if(_measurementAreasByIDs.count(id) == 0) {
        LOG_ERROR("Measurement id [%d] not found.", id);
        exit(EXIT_FAILURE);
    }
    return _measurementAreasByIDs[id];
}

const std::vector<polygon_2d> & ArgumentParser::GetGeometry() const
{
    return _geometry;
}

bool ArgumentParser::IsInMeasureArea(MeasurementArea_L * line, MeasurementArea_B * area)
{
    double lx1 = line->_lineStartX;
    double ly1 = line->_lineStartY;
    double lx2 = line->_lineEndX;
    double ly2 = line->_lineEndY;

    point_2d Line_pt0(lx1, ly1);
    point_2d Line_pt1(lx2, ly2);

    return covered_by(Line_pt0, area->_poly) && covered_by(Line_pt1, area->_poly) &&
           line->_zPos == area->_zPos;
}
