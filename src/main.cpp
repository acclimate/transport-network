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

#include <csv-parser.h>
#include <gdal_alg.h>
#include <gdal_priv.h>
#include <nvector.h>
#include <ogrsf_frmts.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <netcdf>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include "settingsnode.h"
#include "version.h"
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
#include <tqdm/tqdm.h>
#endif
using namespace netCDF;
using namespace netCDF::exceptions;

#define TYPE_REGION 0
#define TYPE_PORT 1
#define TYPE_SEA 2
#define TYPE_COUNT 3

class ProgressBar {
  protected:
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
    std::unique_ptr<tqdm::RangeTqdm<std::size_t>> it;
#endif

  public:
    ProgressBar(std::string description, std::size_t length) {
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
        tqdm::Params p;
        p.desc = description;
        p.ascii = "";
        p.f = stdout;
        it.reset(new tqdm::RangeTqdm<std::size_t>{tqdm::RangeIterator<std::size_t>(length), tqdm::RangeIterator<std::size_t>(length, length), p});
#endif
    }
    inline void tick() {
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
#pragma omp critical(output)
        { ++(*it); }
#endif
    }
    ~ProgressBar() {
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
        it->close();
#endif
    }
};

int main(int argc, char* argv[]) {
#ifndef DEBUG
    try {
#endif
        std::ifstream settings_file(argv[1]);
        if (!settings_file) {
            throw std::runtime_error("cannot open settings file");
        }
        settings::SettingsNode settings(settings_file);
        const settings::SettingsNode& input = settings["input"];
        if (!settings["output"].has("graphdot") && !settings["output"].has("netcdf")) {
            throw std::runtime_error("no output defined");
        }

        // Read ISO3->ISO2 mapping
        std::unordered_map<std::string, std::string> iso_map;
        {
            std::ifstream infile(input["iso3_to_iso2"].as<std::string>());
            if (!infile) {
                throw std::runtime_error("could not open iso mapping file");
            }
            csv::Parser parser(infile);
            do {
                const auto c = parser.read<std::string, std::string>();
                iso_map.emplace(std::get<0>(c), std::get<1>(c));
            } while (parser.next_row());
        }

        // Read regions to ignore
        std::set<std::string> ignores;
        for (const auto& parameters : settings["ignore"].as_sequence()) {
            ignores.insert(parameters.template as<std::string>());
        }

        // For all objects:
        std::vector<std::string> ids;
        std::vector<double> lat;
        std::vector<double> lon;
        std::vector<unsigned char> type;

        // Only for regions:
        std::vector<OGRFeature*> features;
        std::vector<OGRGeometry*> geometries;
        std::vector<OGRPoint> centroids;

        // Read regions
        GDALDataset* infile;
        {
            GDALAllRegister();
            const std::string shapefilename = input["gdal"].as<std::string>();
            infile = static_cast<GDALDataset*>(GDALOpenEx(shapefilename.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr));
            if (!infile) {
                throw std::runtime_error("could not open shape file");
            }
            for (int part = 0; part < input["level"].as<int>(); ++part) {
                std::string layername = "adm" + std::to_string(part);
                OGRLayer* inlayer = infile->GetLayerByName(layername.c_str());
                if (!inlayer) {
                    throw std::runtime_error("could not read layer from shape file");
                }
                inlayer->ResetReading();
                const std::size_t length = inlayer->GetFeatureCount();
                std::string encoding;
                if (layername == "adm0") {
                    encoding = "ISO";
                } else if (layername == "adm1") {
                    encoding = "HASC_1";
                }
                ProgressBar progress("Input " + layername, length);
                inlayer->ResetReading();
#pragma omp parallel for default(shared) schedule(guided)
                for (std::size_t i = 0; i < length; ++i) {
                    OGRFeature* feature;
#pragma omp critical(feature)
                    { feature = inlayer->GetFeature(i + 1); }
                    std::string id = feature->GetFieldAsString(encoding.c_str());
                    if (ignores.find(id) != ignores.end()) {
                        OGRFeature::DestroyFeature(feature);
                    } else {
                        if (encoding == "HASC_1") {
                            if (id.length() != 5 || id[2] != '.') {
                                const std::string& iso_quant = iso_map.at(feature->GetFieldAsString("ISO"));
                                int n = feature->GetFieldAsInteger("ID_1");
                                id = iso_quant + "." + (n < 10 ? "0" : "") + std::to_string(n);
                            }
                        }
                        OGRGeometry* geometry = feature->GetGeometryRef()->SimplifyPreserveTopology(input["resolution"].as<double>());
                        OGRPoint centroid;
                        geometry->Centroid(&centroid);
#pragma omp critical(vectors)
                        {
                            features.push_back(feature);
                            geometries.push_back(geometry);
                            ids.push_back(id);
                            type.push_back(TYPE_REGION);
                            lat.push_back(centroid.getY());
                            lon.push_back(centroid.getX());
                        }
                    }
                    progress.tick();
                }
            }
        }

        const std::size_t regions_count = ids.size();

        std::vector<std::size_t> sea_object_indices;

        // Read sea objects
        {
            std::ifstream infile(input["sea_routes"].as<std::string>());
            if (!infile) {
                throw std::runtime_error("could not open sea routes file");
            }
            csv::Parser parser(infile);
            std::size_t row = 0;
            do {
                const auto c = parser.read<std::string, double, double, std::string>();
                const auto id = std::get<0>(c);
                if (std::get<3>(c) == "sea" || std::get<3>(c) == "port") {
                    ids.push_back(id);
                    sea_object_indices.push_back(row + regions_count);
                    lat.push_back(std::get<1>(c));
                    lon.push_back(std::get<2>(c));
                    if (std::get<3>(c) == "sea") {
                        type.push_back(TYPE_SEA);
                    } else { // std::get<3>(c) == "port"
                        type.push_back(1);
                    }
                } else if (std::get<3>(c) == "region") {
                    const auto index = std::find(std::begin(ids), std::end(ids), id);
                    if (index == std::end(ids)) {
                        throw std::runtime_error("unknown region " + id);
                    }
                    sea_object_indices.push_back(index - std::begin(ids));
                } else {
                    throw std::runtime_error("unknown object type " + std::get<3>(c));
                }
                ++row;
            } while (parser.next_row());
        }

        const std::size_t sea_objects_count = sea_object_indices.size();

        nvector<unsigned char, 2> matrix(0, ids.size(), ids.size());
        std::vector<bool> connected(ids.size(), false);

        // Read sea matrix
        if (sea_objects_count > 0) {
            ProgressBar progress("Sea Matrix", sea_objects_count);
            std::ifstream infile(input["sea_matrix"].as<std::string>());
            if (!infile) {
                throw std::runtime_error("could not open sea matrix file");
            }
            csv::Parser parser(infile);
            std::size_t row = 0;
            do {
                std::size_t col = 0;
                do {
                    const auto value = parser.read<unsigned char>();
                    const auto i = sea_object_indices[row];
                    const auto j = sea_object_indices[col];
                    matrix(i, j) = value;
                    matrix(j, i) = value;
                    if (value) {
                        connected[i] = true;
                        connected[j] = true;
                    }
                    ++col;
                } while (parser.next_col() && col < sea_objects_count);
                ++row;
                progress.tick();
            } while (parser.next_row() && row < sea_objects_count);
        }

        {
            ProgressBar progress("Ports", ids.size());
#pragma omp parallel for default(shared) schedule(dynamic)
            for (std::size_t region = 0; region < regions_count; region++) {
                OGRGeometry* geometry1 = geometries[region];
                for (std::size_t port = regions_count; port < ids.size(); port++) {
                    if (type[port] == TYPE_PORT) {
                        OGRPoint geometry2 = OGRPoint(lon[port], lat[port]);
                        bool inside = geometry1->Contains(&geometry2);
                        if (inside) {
                            connected[region] = true;
                            connected[port] = true;
                            matrix(region, port) = 1;
                            matrix(port, region) = 1;
                        }
                    }
                }
                progress.tick();
            }
        }

        {
            const auto total = regions_count * (regions_count - 1) / 2;
            ProgressBar progress("Region network", total);

            std::vector<std::pair<std::size_t, std::size_t>> pairs;
            pairs.reserve(total);
            for (std::size_t i = 0; i < regions_count; ++i) {
                for (std::size_t j = i + 1; j < regions_count; ++j) {
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
                bool touches = geometry1->Intersects(geometry2);
                if (touches) {
                    const std::string& id1 = ids[i];
                    const std::string& id2 = ids[j];
                    connected[i] = true;
                    connected[j] = true;
                    matrix(i, j) = 1;
                    matrix(j, i) = 1;
                }
                progress.tick();
            }
        }

        if (settings["output"].has("graphdot")) {
            std::ofstream file(settings["output"]["graphdot"].as<std::string>());
            file << "graph {\n";
            for(const auto& id : ids) {
                file << "    // " << id << " [label=\"" << id << "\"];\n";
            }
            for (std::size_t i = 0; i < ids.size(); ++i) {
                for (std::size_t j = 0; j < ids.size(); ++j) {
                    file << "    " << ids[i] << " -- " << ids[j] << ";\n";
                }
            }
            for (std::size_t i = 0; i < ids.size(); ++i) {
                if (!connected[i]) {
                    file << "    // " << ids[i] << " unconnected\n";
                }
            }
            file << "}\n";
            file.close();
        }

        if (settings["output"].has("netcdf")) {
            NcFile file(settings["output"]["netcdf"].as<std::string>(), NcFile::replace, NcFile::nc4);

            NcDim indDim = file.addDim("index", ids.size());
            NcDim typDim = file.addDim("typeindex", TYPE_COUNT);

            NcVar indVar = file.addVar("index", ncString, indDim);
            NcVar mapVar = file.addVar("typeindex", ncString, typDim);
            NcVar typVar = file.addVar("type", ncByte, indDim);
            NcVar latVar = file.addVar("latitude", ncDouble, indDim);
            NcVar lonVar = file.addVar("longitude", ncDouble, indDim);
            NcVar matVar = file.addVar("connections", ncByte, {indDim, indDim});

            latVar.putAtt("units", "degrees_north");
            lonVar.putAtt("units", "degrees_east");

            std::vector<const char*> type_map = {"region", "port", "sea"};
            mapVar.putVar(&type_map[0]);

            std::vector<const char*> ids_c(ids.size());
            for (std::size_t i = 0; i < ids.size(); ++i) {
                ids_c[i] = ids[i].c_str();
            }
            indVar.putVar(&ids_c[0]);

            typVar.putVar(&type[0]);
            latVar.putVar(&lat[0]);
            lonVar.putVar(&lon[0]);
            matVar.putVar(&matrix(0, 0));

            file.close();
        }

        {
            ProgressBar progress("Cleaning up", regions_count);
            for (std::size_t i = 0; i < regions_count; ++i) {
                OGRFeature::DestroyFeature(features[i]);
                progress.tick();
            }
            features.clear();
            geometries.clear();
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
