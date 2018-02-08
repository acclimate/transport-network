/*
  Copyright (C) 2017 Sven Willner <sven.willner@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published
  by the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <gdal_alg.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <netcdf>
#include <stdexcept>
#include <vector>
#include "nvector.h"
#include "version.h"
#include "csv-parser.h"
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
#include "tqdm/tqdm.h"
#endif
using namespace netCDF;
using namespace netCDF::exceptions;

#define TYPE_COUNT  3
#define TYPE_REGION 0
#define TYPE_PORT   1
#define TYPE_SEA    2



inline double distance(const OGRPoint& p1, const OGRPoint& p2) {
    const auto R = 6371;
    const auto PI = 3.14159265;
    const auto latsin = std::sin((p1.getY() - p2.getY()) * PI / 360);
    const auto lonsin = std::sin((p1.getX() - p2.getX()) * PI / 360);
    const auto a = latsin * latsin + cos(p1.getY() * PI / 180) * cos(p2.getY() * PI / 180) * lonsin * lonsin;
    return 2 * R * std::atan2(std::sqrt(a), std::sqrt(1 - a));
}

static inline void rtrim(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int c) {
                                                   return !std::isspace(c);
                                               }).base(), s.end());
}

template <typename T>
bool inside_bool(const T& object,const std::vector<T>& list) {
    auto res = std::find(list.begin(), list.end(), object);
    return res != list.end();
}

template <typename T>
std::size_t inside_number(const T& object,const std::vector<T>& list) {
    auto res = std::find(list.begin(), list.end(), object);
    if (res == list.end()) {
        throw std::runtime_error("fuck");
    } else {
        return res - list.begin();
    }
}


int main(int argc, char* argv[]) {
    
#ifndef DEBUG
    try {
#endif
        std::vector<const char*> type_map;
        type_map.push_back("region");
        type_map.push_back("port");
        type_map.push_back("sea");

        std::cout << "begin" << std::endl;
        std::vector<std::string> ind;
        std::vector<std::string> ind_sea;
        std::vector<double>     lat;
        std::vector<double>     lon;
        std::vector<unsigned char> type;
        std::vector<const char*> ind_c;
        
        std::ifstream infile_sea("/home/kikula/Documents/p/Acclimate/Projekte/Infrastructure/Shipping_Routes/water_objects.csv");
        csv::Parser parser(infile_sea);
        
        do {
            const auto c = parser.read<std::string, double, double,std::string>();
            if (std::get<3>(c) == "sea" || std::get<3>(c) == "port") {             
                ind.push_back(std::get<0>(c));
                lat.push_back(std::get<1>(c));
                lon.push_back(std::get<2>(c));
                if ( std::get<3>(c) == "sea") {
                    type.push_back(2);
                }
                if ( std::get<3>(c) == "port") {
                    type.push_back(1);
                }
            }
            ind_sea.push_back(std::get<0>(c));
        } while (parser.next_row());
        
        std::vector<std::string> adm_string;
        adm_string.push_back("adm0");
        adm_string.push_back("adm1");
        
        std::vector<std::string> coding_string;
        adm_string.push_back("ISO");
        adm_string.push_back("HASC-1");
        
        
        
        std::vector<OGRFeature*> features(0);
        std::vector<OGRGeometry*> geometries(0);
        std::vector<OGRPoint> centroids(0);
        std::vector<std::string> ids(0);
        std::vector<bool> connected(0,false);
        
        
        std::ofstream file("graph.dot");
        file << "graph {\n";
        GDALAllRegister();
        const std::string shapefilename = argv[1];
        auto infile = static_cast<GDALDataset*>(GDALOpenEx(shapefilename.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr));
        for (int part = 0; part < 2; part++) {
            
            //~ const std::string layername = argv[2];
            std::string layername = adm_string[part];
            if (!infile) {
                throw std::runtime_error("could not open shape file");
            }
            OGRLayer* inlayer = infile->GetLayerByName(layername.c_str());
            inlayer->ResetReading();
            if (!inlayer) {
                throw std::runtime_error("could not read layer from shape file");
            }
            const std::size_t length = inlayer->GetFeatureCount(); // ACHTUNG DA SIZE WICHTIG IST; MUSS ES WOANDERS ANGEPASST WERDEN
            std::size_t starter = ids.size();
            
            features.resize(features.size() + length);
            geometries.resize(geometries.size() + length);
            centroids.resize(centroids.size() + length);
            ids.resize(ids.size() + length);
            connected.resize(connected.size() + length);
        {
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
            tqdm::Params p;
            p.desc = "Input " + std::to_string(part + 1);
            p.ascii = "";
            p.f = stdout;
            tqdm::RangeTqdm<std::size_t> it{tqdm::RangeIterator<std::size_t>(length), tqdm::RangeIterator<std::size_t>(length, length), p};
#endif
            inlayer->ResetReading();
#pragma omp parallel for default(shared) schedule(guided)
            for (std::size_t i = starter; i < ids.size(); ++i) {
                OGRFeature* feature;
#pragma omp critical(feature)
                { feature = inlayer->GetFeature(i - starter  + 1); }
                features[i] = feature;
                //~ OGRGeometry* geometry = feature->GetGeometryRef();  //->SimplifyPreserveTopology(0.001);
                OGRGeometry* geometry = feature->GetGeometryRef()->SimplifyPreserveTopology(0.001);
                geometry->Centroid(&centroids[i]);
                geometries[i] = geometry;  //->ConvexHull();
                std::string id = feature->GetFieldAsString("ISO");
                if (id.empty()) {
                    id = std::to_string(i);
                } else {
                    rtrim(id);
                    std::replace(id.begin(), id.end(), '.', '_');
                }
#pragma omp critical(foutput)
                {
                    if (id.find("-") == std::string::npos) {
                        file << "    // " << id << " [label=\"" << id << "\"];\n" << std::flush;
                    }
                }
                ids[i] = id;
                ind.push_back(id);
                type.push_back(0);
                lat.push_back(centroids[i].getY());
                lon.push_back(centroids[i].getX());
                // std::cout << id << " " << centroids[i].getX() << " " << centroids[i].getY() << std::endl;
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
                ++it;
#endif
            }
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
            it.close();
#endif
        }
        }
        const std::size_t size = ids.size();
        
        nvector<unsigned char,2> matrix(0,ind.size(),ind.size());
        {
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
            tqdm::Params p;
            p.desc = "Sea Matrix";
            p.ascii = "";
            p.f = stdout;
            const auto total = ind_sea.size();
            tqdm::RangeTqdm<std::size_t> it{tqdm::RangeIterator<std::size_t>(total), tqdm::RangeIterator<std::size_t>(total, total), p};
#endif
        
        std::ifstream infile_matrix("/home/kikula/Documents/p/Acclimate/Projekte/Infrastructure/Shipping_Routes/matrix.csv");
        csv::Parser parser2(infile_matrix);
        int k = 0;
        do {
            const auto row = inside_number(ind_sea[k],ind);
            int p = 0;
            do {
                const auto col = inside_number(ind_sea[p],ind);
                matrix(row,col) = parser2.read<unsigned char>();
                if (matrix(row,col)) {
#pragma omp critical(foutput)
                    { file << "    " << ind[row] << " -- " << ind[col] << ";\n" << std::flush; }
                }
                p++;
            } while (parser2.next_col());
            k++;
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
            ++it;
#endif           
        } while (parser2.next_row());
        
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
            it.close();
#endif
        }    

        return 0 ;
        {
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
            tqdm::Params p;
            p.desc = "Ports";
            p.ascii = "";
            p.f = stdout;
            const auto total = ids.size();
            tqdm::RangeTqdm<std::size_t> it{tqdm::RangeIterator<std::size_t>(total), tqdm::RangeIterator<std::size_t>(total, total), p};
#endif
            const auto PORT = TYPE_PORT;
#pragma omp parallel for default(shared) schedule(dynamic)
        for (std::size_t country = 0; country < ids.size(); country++) {
            OGRGeometry* geometry1 = geometries[country];
            const auto i = inside_number(ids[country], ind);
            for (std::size_t port = 0; port  < ind.size() ; port++) {
                if (type[port] == PORT) {
                    OGRPoint geometry2 = OGRPoint(lon[port],lat[port]);
                    bool inside = 
                                //(ids[country] == "DEU") && 
                                geometry1->Contains(&geometry2);
                    if (inside) {
#pragma omp critical(foutput)
                    { file << "    " << ind[i] << " -- " << ind[port] << ";\n" << std::flush; }
                        matrix(i,port) = 1;
                        matrix(port,i) = 1;
                    }
                }
            }
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
            ++it;
#endif
        }
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
        it.close();
#endif
        }

        {
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
            tqdm::Params p;
            p.desc = "Network";
            p.ascii = "";
            p.f = stdout;
            const auto total = size * (size - 1) / 2;
            tqdm::RangeTqdm<std::size_t> it{tqdm::RangeIterator<std::size_t>(total), tqdm::RangeIterator<std::size_t>(total, total), p};
#endif

            std::vector<std::pair<std::size_t, std::size_t>> pairs;
            pairs.reserve(total);
            for (std::size_t i = 0; i < size; ++i) {
                for (std::size_t j = i + 1; j < size; ++j) {
                    pairs.emplace_back(std::make_pair(i, j));
                }
            }
            std::random_shuffle(std::begin(pairs), std::end(pairs));
        
#pragma omp parallel for default(shared) schedule(dynamic)
            for (std::size_t k = 0; k < pairs.size(); ++k) {
                const auto p = pairs[k];
                const auto i = p.first;
                const auto j = p.second;
                OGRGeometry* geometry1 = geometries[i];
                OGRGeometry* geometry2 = geometries[j];
                // OGRFeature* feature1 = features[i];
                // OGRFeature* feature2 = features[j];
                const std::string& id1 = ids[i];  // feature1->GetFieldAsString("ISO");
                const std::string& id2 = ids[j];  // feature2->GetFieldAsString("ISO");
                bool touches =     
                            //~ (id1 == "DEU" && id2 == "FRA") &&
                            geometry1->Intersects(geometry2);// || geometry1->Touches(geometry2));
                if (touches) {
                    connected[i] = true;
                    connected[j] = true;
#pragma omp critical(foutput)
                    { file << "    " << id1 << " -- " << id2 << ";\n" << std::flush; }
                    matrix(inside_number(ids[i],ind),inside_number(ids[j],ind)) = 1;
                    matrix(inside_number(ids[j],ind),inside_number(ids[i],ind)) = 1;
                }
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
#pragma omp critical(output)
                { ++it; }
#endif
            }
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
            it.close();
#endif
        }
        for (std::size_t i = 0; i < size; ++i) {
            if (!connected[i]) {
                file << "    // " << ids[i] << " unconnected\n";
            }
        }
        
        for (std::size_t i = 0; i < ind.size(); i++) {
            ind_c.push_back(ind[i].c_str());
        }
        file << "}\n";
        file.close();
        
        try {
            NcFile test("routes_test.nc", NcFile::replace,NcFile::nc4);

            NcDim indDim = test.addDim("index", ind.size());
            NcDim typDim = test.addDim("typeindex", TYPE_COUNT);

            NcVar indVar = test.addVar("index",  ncString, indDim);
            NcVar mapVar = test.addVar("typeindex",  ncString, typDim);
            NcVar typVar = test.addVar("type",  ncByte, indDim);
            NcVar latVar = test.addVar("latitude", ncDouble, indDim);
            NcVar lonVar = test.addVar("longitude",  ncDouble, indDim);
            NcVar matVar = test.addVar("connections", ncByte, {indDim,indDim}) ;

            indVar.putAtt("units", "");
            mapVar.putAtt("units", "");
            typVar.putAtt("units", "");
            latVar.putAtt("units", "degrees_north");
            lonVar.putAtt("units", "degrees_east");
            matVar.putAtt("units", "");
            
            indVar.putVar(&ind_c[0]);
            mapVar.putVar(&type_map[0]);
            typVar.putVar(&type[0]);
            latVar.putVar(&lat[0]);
            lonVar.putVar(&lon[0]);
            matVar.putVar(&matrix.data()[0]);
            std::cout << "for the win" << std::endl;
        

        }
        catch(NcException& e)
        {
            std::cout << e.what() <<  std::endl;
        }

        {
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
            tqdm::Params p;
            p.desc = "Cleaning up";
            p.ascii = "";
            p.f = stdout;
            tqdm::RangeTqdm<std::size_t> it{tqdm::RangeIterator<std::size_t>(size), tqdm::RangeIterator<std::size_t>(size, size), p};
#endif
            for (std::size_t i = 0; i < size; ++i) {
                OGRFeature::DestroyFeature(features[i]);
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
                ++it;
#endif
            }
            features.clear();
            geometries.clear();
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
            it.close();
#endif
        }
        
        GDALClose(infile);
        
#ifndef DEBUG
    } catch (std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return 255;
    }
#endif
    

    return 0;
}
