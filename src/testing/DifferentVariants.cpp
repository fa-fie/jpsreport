#include "../Analysis.h"
#include "../general/Logger.h"
#include "DifferentVariants.h"

#include <fstream>
#include <iostream>

using std::ofstream;
using std::string;
using std::vector;

DifferentVariants::DifferentVariants(MeasurementArea_B * area)
{
    _realV           = 1; 
    // set this to the real velocity to caluclate the error later
    // if no real velocity is given, use std::numeric_limits<double>::quiet_NaN()
    // in that case, the average velocity is calculated instead of the error
    _fps             = 16;
    _areaForTesting  = area;
    _measureAreaId   = boost::lexical_cast<string>(_areaForTesting->_id);
    _dx              = area->_length;
}

DifferentVariants::~DifferentVariants() = default;

void DifferentVariants::RunTests(const PedData & peddata) 
{
    _outputLocation = peddata.GetOutputLocation();
    _peds_t         = peddata.GetPedIDsByFrameNr();
    _xCor           = peddata.GetXCor();
    _yCor           = peddata.GetYCor();
    _minFrame       = peddata.GetMinFrame();
    _fps            = peddata.GetFps();
    _firstFrame     = peddata.GetFirstFrame();

    if(_dx < 0) {
        LOG_WARNING("The measurement area length in movement direction for testing variants is not assigned!");
        exit(EXIT_FAILURE);
    }

    /*
    Explanation for the different variants:
    1   last frame before entrance (i.e. last frame not in MA) -> exit equivalently (not including on line)
    2   last frame before entrance (i.e. last frame not in MA) -> exit equivalently (including on line)
    3   first frame after entrance (i.e. first frame fully in MA) -> exit equivalently (not including on line)
    4   first frame after entrance (i.e. first frame fully in MA) -> exit equivalently (including on line)
    5   combination of 1 and 2 -> first frame after entrance (not on line) -> last frame before exit (not on line) 
    6   combination of 1 and 2 -> first frame after entrance (or on line) -> last frame before exit (or on line)
    7   frame with smallest distance to entrance or exit (if equal -> last frame)
    8   frame with smallest distance to entrance or exit (if equal -> first frame)
    */

    std::ofstream errorFile;
    std::ofstream avgFile;
    if(!isnan(_realV)) {
        // if data from real experiments is used, _realV might not be a certain value
        errorFile = GetFile("error_" + std::to_string((int) _fps) + "_fps", "error");
        errorFile << "#framerate:\t" << _fps << "\n#real velocity (m/s):\t" << _realV
                  << "\n\n#variant\terror (deviation from real velocity)\n";
    } else {
        avgFile = GetFile("avg_" + std::to_string((int) _fps) + "_fps", "average_v");
        avgFile << "#framerate:\t" << _fps << "\n\n#variant\taverage velocity\n";
    }

    int numFrames     = peddata.GetNumFrames();
    int numPeds       = peddata.GetNumPeds();
    string foldername = std::to_string((int) _fps) + "_fps";
    for(int variant = 1; variant < 9; variant++) {
        vector<vector<int>> TinTout = GetTinTout(numFrames, _areaForTesting->_poly, numPeds, variant);
        vector<int> tIn             = TinTout[0];
        vector<int> tOut            = TinTout[1];
        std::ofstream file          = GetFile("variant_" + std::to_string(variant), foldername);
        file << "#framerate:\t" << _fps;
        file << "\n#length in movement direction (m):\t" << _dx;
        file << "\n\n#person index\tvelocity\tdifference in frames\tentry frame\texit frame\n";
        double sumVelocity = 0;
        int pedsInArea     = 0;
        for(int ped = 0; ped < numPeds; ped++) {
            if (tOut[ped] > 0) {
                // if pedestrian crossed the measurement area
                double velocity = _dx / ((tOut[ped] - tIn[ped]) / _fps);
                file << ped << "\t" << velocity << "\t" << (tOut[ped] - tIn[ped]) << "\t"
                     << tIn[ped] + peddata.GetMinFrame() << "\t"
                     << tOut[ped] + peddata.GetMinFrame() << "\n";
                pedsInArea++;
                sumVelocity += velocity;
            }
        }
        double avgVelocity = sumVelocity / pedsInArea;
        file.close();
        if(!isnan(_realV)) {
            double deviation   = abs(avgVelocity - _realV);
            errorFile << variant << "\t" << deviation << "\n";
        } else {
            avgFile << variant << "\t" << avgVelocity << "\n";
        }
    }
    errorFile.close();
}

vector<vector<int>> DifferentVariants::GetTinTout(
    int numFrames,
    polygon_2d polygon,
    int numPeds,
    int variant)
{
    vector<bool> IsinMeasurezone(numPeds, false);
    vector<int> tIn(numPeds, 0);
    vector<int> tOut(numPeds, 0);

    segment entrance;
    segment exit;
    for(int frameNr = 0; frameNr < numFrames; frameNr++) {
        vector<int> ids = _peds_t[frameNr];
        for(int ID : ids) {
            int x = _xCor(ID, frameNr);
            int y = _yCor(ID, frameNr);
            int nextX, nextY;
            if((frameNr + 1) < numFrames) {
                nextX = _xCor(ID, frameNr + 1);
                nextY = _yCor(ID, frameNr + 1);
            }

            // switch different variants
            switch(variant) {
                case 1:
                    if((frameNr + 1) < numFrames && !(nextX == 0 && nextY == 0)) {
                        if(covered_by(make<point_2d>((nextX), (nextY)), polygon) &&
                           !(IsinMeasurezone[ID])) {
                            tIn[ID]             = frameNr;
                            IsinMeasurezone[ID] = true;
                        } else if(
                            (!within(make<point_2d>((nextX), (nextY)), polygon)) &&
                            IsinMeasurezone[ID]) {
                            tOut[ID]            = frameNr;
                            IsinMeasurezone[ID] = false;
                        }
                    }
                    break;
                case 2:
                    if((frameNr + 1) < numFrames && !(nextX == 0 && nextY == 0)) {
                        if(within(make<point_2d>((nextX), (nextY)), polygon) &&
                           !(IsinMeasurezone[ID])) {
                            tIn[ID]             = frameNr;
                            IsinMeasurezone[ID] = true;
                        } else if(
                            (!covered_by(make<point_2d>((nextX), (nextY)), polygon)) &&
                            IsinMeasurezone[ID]) {
                            tOut[ID]            = frameNr;
                            IsinMeasurezone[ID] = false;
                        }
                    }
                    break;
                case 3:
                    if(within(make<point_2d>((x), (y)), polygon) && !(IsinMeasurezone[ID]) &&
                       !(nextX == 0 && nextY == 0)) {
                        tIn[ID]             = frameNr;
                        IsinMeasurezone[ID] = true;
                    } else if(
                        (!covered_by(make<point_2d>((x), (y)), polygon)) &&
                        IsinMeasurezone[ID]) {
                        tOut[ID]            = frameNr;
                        IsinMeasurezone[ID] = false;
                    }
                    break;
                case 4:
                    if(covered_by(make<point_2d>((x), (y)), polygon) && !(IsinMeasurezone[ID]) &&
                       !(nextX == 0 && nextY == 0)) {
                        tIn[ID]             = frameNr;
                        IsinMeasurezone[ID] = true;
                    } else if(
                        (!within(make<point_2d>((x), (y)), polygon)) &&
                        IsinMeasurezone[ID]) {
                        tOut[ID]            = frameNr;
                        IsinMeasurezone[ID] = false;
                    }
                    break;
                case 5:
                    if((frameNr + 1) < numFrames && !(nextX == 0 && nextY == 0)) {
                        if(within(make<point_2d>((x), (y)), polygon) &&
                           !(IsinMeasurezone[ID])) {
                            tIn[ID]             = frameNr;
                            IsinMeasurezone[ID] = true;
                        } else if(
                            (!within(make<point_2d>((nextX), (nextY)), polygon)) &&
                            IsinMeasurezone[ID]) {
                            tOut[ID]            = frameNr;
                            IsinMeasurezone[ID] = false;
                        }
                    }
                    break;
                case 6:
                    if((frameNr + 1) < numFrames && !(nextX == 0 && nextY == 0)) {
                        if(covered_by(make<point_2d>((x), (y)), polygon) && !(IsinMeasurezone[ID])) {
                            tIn[ID]             = frameNr;
                            IsinMeasurezone[ID] = true;
                        } else if(
                            (!covered_by(make<point_2d>((nextX), (nextY)), polygon)) &&
                            IsinMeasurezone[ID]) {
                            tOut[ID]            = frameNr;
                            IsinMeasurezone[ID] = false;
                        }
                    }
                    break;
                case 7:
                    if((frameNr + 1) < numFrames) {
                        if(!covered_by(make<point_2d>(x, y), polygon) &&
                           !IsinMeasurezone[ID] &&
                           covered_by(make<point_2d>(nextX, nextY), polygon) && !(nextX == 0 && nextY == 0)) {
                            // the area is entered between these frames
                            for(int p = 0; p < 4; p++) {
                                segment edge0(polygon.outer()[p], polygon.outer()[p + 1]);
                                segment edge1(make<point_2d>(x, y), make<point_2d>(nextX, nextY));
                                if(intersects(edge0, edge1)) {
                                    entrance = edge0;
                                    if(p < 2) {
                                        segment edge2(
                                            polygon.outer()[p + 2], polygon.outer()[p + 3]);
                                        exit = edge2;
                                    } else if (p == 2) {
                                        segment edge2(polygon.outer()[4], polygon.outer()[1]);
                                        exit = edge2;
                                    } else if (p == 3) {
                                        segment edge2(polygon.outer()[1], polygon.outer()[2]);
                                        exit = edge2;
                                    }
                                }
                            }
                            double dist1 = boost::geometry::distance(make<point_2d>(x, y), entrance);
                            double dist2 =
                                boost::geometry::distance(make<point_2d>(nextX, nextY), entrance);
                            if(dist1 < dist2) {
                                tIn[ID] = frameNr;
                            } else {
                                tIn[ID] = frameNr + 1;
                            }
                            IsinMeasurezone[ID] = true;
                        } else if(
                            covered_by(make<point_2d>((x), (y)), polygon) &&
                            IsinMeasurezone[ID] &&
                            !covered_by(make<point_2d>((nextX), (nextY)), polygon)) {
                            // the area is exited between these frames
                            double dist1 = boost::geometry::distance(make<point_2d>(x, y), exit);
                            double dist2 =
                                boost::geometry::distance(make<point_2d>(nextX, nextY), exit);
                            if(dist1 < dist2) {
                                tOut[ID] = frameNr;
                            } else {
                                tOut[ID] = frameNr + 1;
                            }
                            IsinMeasurezone[ID] = false;
                        }
                    }
                    break;
                case 8:
                    if((frameNr + 1) < numFrames) {
                        if(!covered_by(make<point_2d>(x, y), polygon) && !IsinMeasurezone[ID] &&
                           covered_by(make<point_2d>(nextX, nextY), polygon) &&
                           !(nextX == 0 && nextY == 0)) {
                            // the area is entered between these frames
                            for(int p = 0; p < 4; p++) {
                                segment edge0(polygon.outer()[p], polygon.outer()[p + 1]);
                                segment edge1(make<point_2d>(x, y), make<point_2d>(nextX, nextY));
                                if(intersects(edge0, edge1)) {
                                    entrance = edge0;
                                    entrance = edge0;
                                    if(p < 2) {
                                        segment edge2(
                                            polygon.outer()[p + 2], polygon.outer()[p + 3]);
                                        exit = edge2;
                                    } else if(p == 2) {
                                        segment edge2(polygon.outer()[4], polygon.outer()[1]);
                                        exit = edge2;
                                    } else if(p == 3) {
                                        segment edge2(polygon.outer()[1], polygon.outer()[2]);
                                        exit = edge2;
                                    }
                                }
                            }
                            double dist1 =
                                boost::geometry::distance(make<point_2d>(x, y), entrance);
                            double dist2 =
                                boost::geometry::distance(make<point_2d>(nextX, nextY), entrance);
                            if(dist1 <= dist2) {
                                tIn[ID] = frameNr;
                            } else {
                                tIn[ID] = frameNr + 1;
                            }
                            IsinMeasurezone[ID] = true;
                        } else if(
                            covered_by(make<point_2d>((x), (y)), polygon) && IsinMeasurezone[ID] &&
                            !covered_by(make<point_2d>((nextX), (nextY)), polygon)) {
                            // the area is exited between these frames
                            double dist1 = boost::geometry::distance(make<point_2d>(x, y), exit);
                            double dist2 =
                                boost::geometry::distance(make<point_2d>(nextX, nextY), exit);
                            if(dist1 <= dist2) {
                                tOut[ID] = frameNr;
                            } else {
                                tOut[ID] = frameNr + 1;
                            }
                            IsinMeasurezone[ID] = false;
                        }
                    }
                    break;
            }
        }
    }

    vector<vector<int>> output;
    output.push_back(tIn);
    output.push_back(tOut);
    return output;
}

std::ofstream DifferentVariants::GetFile(string fname, string foldername)
{
    fs::path tmp(fname + ".dat");
    tmp               = _outputLocation / "Testing_Variants" / foldername / tmp.string();
    string filename   = tmp.string();
    fs::path filepath = fs::path(filename.c_str()).parent_path();
    if(fs::is_directory(filepath) == false) {
        if(fs::create_directories(filepath) == false && fs::is_directory(filepath) == false) {
            LOG_ERROR("cannot create the directory <{}>", filepath.string());
            exit(EXIT_FAILURE);
        }
        LOG_INFO("create the directory <{}>", filepath.string());
    }
    std::ofstream file(tmp.string());
    return file;
}
