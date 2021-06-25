#ifndef DIFFERENTVARIANTS_H
#define DIFFERENTVARIANTS_H

#include "../general/Macros.h"
#include "../methods/MeasurementArea.h"
#include "../methods/PedData.h"
#include "tinyxml.h"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

typedef boost::geometry::model::segment<boost::geometry::model::d2::point_xy<double>> segment;

class DifferentVariants
{
public:
    DifferentVariants(MeasurementArea_B * area);
    virtual ~DifferentVariants();
    void RunTests(const PedData & peddata);

private:
    fs::path _outputLocation;
    std::map<int, std::vector<int>> _peds_t;
    ub::matrix<double> _xCor;
    ub::matrix<double> _yCor;
    std::vector<int> _firstFrame;
    int _minFrame;
    float _fps;
    float _realV;

    std::string _measureAreaId;
    MeasurementArea_B * _areaForTesting;
    double _dx;

    std::vector<std::vector<int>> GetTinTout(int numFrames, polygon_2d polygon, int numPeds, int variant);
    std::ofstream GetFile(std::string fname, std::string foldername);
};

#endif /* DIFFERENTVARIANTS_H */
