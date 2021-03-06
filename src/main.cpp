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
        //p.f = stdout;
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

enum Type { REGION = 0, PORT = 1, SEA = 2 };
std::vector<const char*> type_map = {"region", "port", "sea"};

class Location {
  public:
    std::string id;
    std::string name;
    double lat;
    double lon;
    Type type;
    OGRFeature* feature;
    OGRGeometry* geometry;
    std::size_t level;
    Location(std::string id_p,
             std::string name_p,
             double lat_p,
             double lon_p,
             Type type_p,
             OGRFeature* feature_p = nullptr,
             OGRGeometry* geometry_p = nullptr,
             std::size_t level_p = 0)
        : id(id_p), name(name_p), lat(lat_p), lon(lon_p), type(type_p), feature(feature_p), geometry(geometry_p), level(level_p) {}
};

static void print_usage(const char* program_name) {
    std::cerr << "Transport preprocessing for the acclimate model\n"
                 "Version: " ACCLIMATE_TRANSPORT_VERSION
                 "\n\n"
                 "Authors: Sven Willner <sven.willner@pik-potsdam.de>\n"
                 "         Kilian Kuhla <kilian.kuhla@pik-potsdam.de>\n"
                 "\n"
                 "Usage:   "
              << program_name
              << " (<option> | <settingsfile>)\n"
                 "Options:\n"
                 "  -h, --help     Print this help text\n"
                 "  -v, --version  Print version"
              << std::endl;
}

int main(int argc, char* argv[]) {
#ifndef DEBUG
    try {
#endif
        if (argc != 2) {
            print_usage(argv[0]);
            return 1;
        }
        const std::string arg = argv[1];
        if (arg.length() > 1 && arg[0] == '-') {
            if (arg == "--version" || arg == "-v") {
                std::cout << ACCLIMATE_TRANSPORT_VERSION << std::endl;
                return 0;
            } else if (arg == "--help" || arg == "-h") {
                print_usage(argv[0]);
                return 0;
            } else {
                print_usage(argv[0]);
                return 1;
            }
        }
        settings::SettingsNode settings;
        if (arg == "-") {
            std::cin >> std::noskipws;
            settings = settings::SettingsNode(std::cin);
        } else {
            std::ifstream settings_file(arg);
            if (!settings_file) {
                throw std::runtime_error("cannot open settings file");
            }
            settings = settings::SettingsNode(settings_file);
        }

        const settings::SettingsNode& input = settings["input"];
        if (!settings["output"].has("graphdot") && !settings["output"].has("netcdf") && !settings["output"].has("yaml")) {
            throw std::runtime_error("no output defined");
        }

        const settings::SettingsNode iso_map = settings["iso_mapping"];
        const auto resolution = input["resolution"].as<double>();
        const auto port_threshold = input["port_threshold"].as<double>();

        std::vector<Location> locations;
        std::unordered_map<std::string, std::size_t> indices;

        // Read regions to ignore
        std::set<std::string> ignores;
        for (const auto& parameters : settings["ignore"].as_sequence()) {
            ignores.insert(parameters.template as<std::string>());
        }

        // Read regions
        GDALDataset* infile;
        {
            GDALAllRegister();
            const std::string shapefilename = input["gdal"].as<std::string>();
            infile = static_cast<GDALDataset*>(GDALOpenEx(shapefilename.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr));
            if (!infile) {
                throw std::runtime_error("could not open shape file");
            }
            for (std::size_t level = 0; level < input["level"].as<std::size_t>(); ++level) {
                std::string layername = "adm" + std::to_string(level);
                OGRLayer* inlayer = infile->GetLayerByName(layername.c_str());
                if (!inlayer) {
                    throw std::runtime_error("could not read layer from shape file");
                }
                inlayer->ResetReading();
                const std::size_t length = inlayer->GetFeatureCount();
                std::string encoding;
                std::string name_field;
                if (layername == "adm0") {
                    encoding = "ISO";
                    name_field = "NAME_ENGLISH";
                } else if (layername == "adm1") {
                    encoding = "HASC_1";
                    name_field = "NAME_1";
                }
                ProgressBar progress("Input " + layername, length);
                inlayer->ResetReading();
#pragma omp parallel for default(shared) schedule(dynamic)
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
                                const std::string& iso_quant = iso_map[feature->GetFieldAsString("ISO")].as<std::string>();
                                int n = feature->GetFieldAsInteger("ID_1");
                                id = iso_quant + "." + (n < 10 ? "0" : "") + std::to_string(n);
                            }
                        }
                        const std::string& rawname = feature->GetFieldAsString(name_field.c_str());
                        std::string name;
                        const auto p = rawname.find_first_of("\n\r");
                        if (p != std::string::npos) {
                            name = rawname.substr(0, p);
                        } else {
                            name = rawname;
                        }
                        OGRGeometry* geometry;
                        if (resolution > 0) {
                            geometry = feature->GetGeometryRef()->SimplifyPreserveTopology(resolution);
                        } else {
                            geometry = feature->GetGeometryRef();
                        }
                        OGRPoint centroid;
                        geometry->Centroid(&centroid);
#pragma omp critical(locations)
                        {
                            indices.emplace(id, locations.size());
                            locations.push_back(Location{id, name, centroid.getY(), centroid.getX(), Type::REGION, feature, geometry, level});
                        }
                    }
                    progress.tick();
                }
            }
        }

        const auto shapefile_regions_count = locations.size();

        // Read ports
        for (const auto& p : settings["ports"].as_map()) {
            indices.emplace(p.first, locations.size());
            locations.push_back({p.first, p.second["name"].as<std::string>(), p.second["lat"].as<double>(), p.second["lon"].as<double>(), Type::PORT});
        }

        const auto ports_count = locations.size() - shapefile_regions_count;

        // Read seas
        for (const auto& s : settings["seas"].as_map()) {
            indices.emplace(s.first, locations.size());
            locations.push_back({s.first, s.second["name"].as<std::string>(), s.second["lat"].as<double>(), s.second["lon"].as<double>(), Type::SEA});
        }

        // Read additional regions
        if (settings.has("additional_regions")) {
            for (const auto& p : settings["additional_regions"].as_map()) {
                indices.emplace(p.first, locations.size());
                locations.push_back({p.first, p.second["name"].as<std::string>(), p.second["lat"].as<double>(), p.second["lon"].as<double>(), Type::REGION});
            }
        }

        // All location have been read, create connections
        std::vector<unsigned char> matrix(locations.size() * locations.size(), 0);

        // Read additional connections
        if (settings.has("additional_connections")) {
            for (const auto& location : settings["additional_connections"].as_map()) {
                const auto i = indices.find(location.first);
                if (i == std::end(indices)) {
                    std::cout << "Warning: Could not find " + location.first << std::endl;
                } else {
                    for (const auto& other : location.second.as_sequence()) {
                        const std::string id = other.as<std::string>();
                        const auto j = indices.find(id);
                        if (j == std::end(indices)) {
                            std::cout << "Warning: Could not find " + id << std::endl;
                        } else {
                            matrix[i->second * locations.size() + j->second] = 1;
                            matrix[j->second * locations.size() + i->second] = 1;
                        }
                    }
                }
            }
        }

        // Read connections of seas and ports
        for (const auto& type : {"seas", "ports"}) {
            for (const auto& location : settings[type].as_map()) {
                const auto i = indices.find(location.first);
                if (location.second.has("connections")) {
                    for (const auto& other : location.second["connections"].as_sequence()) {
                        const std::string id = other.as<std::string>();
                        const auto j = indices.find(id);
                        if (j == std::end(indices)) {
                            std::cout << "Warning: Could not find " + id << std::endl;
                        } else {
                            matrix[i->second * locations.size() + j->second] = 1;
                            matrix[j->second * locations.size() + i->second] = 1;
                        }
                    }
                }
            }
        }

        // Find regions of ports
        {
            ProgressBar progress("Ports", ports_count);
#pragma omp parallel for default(shared) schedule(dynamic)
            for (std::size_t port = shapefile_regions_count; port < shapefile_regions_count + ports_count; ++port) {
                const auto& port_location = locations[port];
                OGRPoint null_point = OGRPoint(0, 0);
                OGRSpatialReference target_ref;
                target_ref.SetAE(port_location.lat, port_location.lon, 0, 0);
                OGRPoint port_point = OGRPoint(port_location.lon, port_location.lat);
                bool region_found = false;
                for (std::size_t region = 0; region < shapefile_regions_count; ++region) {
                    if (locations[region].geometry->Distance(&port_point) < 90) {
                        OGRGeometry* geometry;
                        OGRCoordinateTransformation* transformation;
                        geometry = locations[region].geometry->clone();
                        transformation = OGRCreateCoordinateTransformation(geometry->getSpatialReference(), &target_ref);
                        geometry->transform(transformation);
                        const auto dist = geometry->Distance(&null_point) / 1000;
                        OGRGeometryFactory::destroyGeometry(geometry);
                        if (dist < port_threshold) {
                            region_found = true;
                            matrix[region * locations.size() + port] = 1;
                            matrix[port * locations.size() + region] = 1;
                        }
                    }
                }
                if (!region_found) {
#pragma omp critical(output)
                    { std::cout << "Warning: no region found for port " << port_location.id << " (" << port_location.name << ")" << std::endl; }
                }
                progress.tick();
            }
        }

        // Find connected regions
        {
            const auto total = shapefile_regions_count * (shapefile_regions_count - 1) / 2;
            ProgressBar progress("Region network", total);

            std::vector<std::pair<std::size_t, std::size_t>> pairs;
            pairs.reserve(total);
            for (std::size_t i = 0; i < shapefile_regions_count; ++i) {
                for (std::size_t j = i + 1; j < shapefile_regions_count; ++j) {
                    pairs.emplace_back(std::make_pair(i, j));
                }
            }
            std::random_shuffle(std::begin(pairs), std::end(pairs));

#pragma omp parallel for default(shared) schedule(dynamic)
            for (std::size_t k = 0; k < pairs.size(); ++k) {
                const auto p = pairs[k];
                const auto i = p.first;
                const auto j = p.second;
                bool touches = locations[i].geometry->Intersects(locations[j].geometry);
                if (touches) {
                    matrix[i * locations.size() + j] = 1;
                    matrix[j * locations.size() + i] = 1;
                }
                progress.tick();
            }
        }

        // Show subnets
        {
            std::vector<int> subnet(locations.size(), -1);
            int subnet_id = 0;
            while(true) {
                bool subnet_found = false;
                std::size_t size = 0;
                std::size_t level;
                for (std::size_t i = 0; i < locations.size(); ++i) {
                    if (subnet[i] < 0 && locations[i].type == Type::REGION) {
                        std::cout << "Subnet:" << std::endl;
                        subnet_found = true;
                        level = locations[i].level;
                        subnet[i] = subnet_id;
                        std::cout << "  " << locations[i].id << " " << locations[i].name << std::endl;
                        ++size;
                        ++subnet_id;
                        break;
                    }
                }
                if (!subnet_found) {
                    break;
                }

                while(true) {
                    bool unconnected_found = false;
                    for (std::size_t i = 0; i < locations.size(); ++i) {
                        if (subnet[i] == subnet_id - 1) {
                            for (std::size_t j = 0; j < locations.size(); ++j) {
                                if (matrix[i * locations.size() + j] && subnet[j] < 0 && (locations[j].type != Type::REGION || level == locations[j].level)) {
                                    unconnected_found = true;
                                    subnet[j] = subnet_id - 1;
                                    std::cout << "  " << locations[j].id << " " << locations[j].name << std::endl;
                                    ++size;
                                }
                            }
                        }
                    }
                    if (!unconnected_found) {
                        std::cout << "  (size: " << size << ")" << std::endl;
                        break;
                    }
                }

                for (std::size_t i = 0; i < locations.size(); ++i) {
                    if (locations[i].type != Type::REGION) {
                        subnet[i] = -1;
                    }
                }
            }
        }

        // YAML output
        if (settings["output"].has("yaml")) {
            std::ofstream file(settings["output"]["yaml"].as<std::string>());
            for (std::size_t i = 0; i < locations.size(); ++i) {
                const auto& location = locations[i];
                file << location.id << ":\n  name: " << location.name << "\n  lat: " << location.lat << "\n  lon: " << location.lon
                     << "\n  type: " << type_map[location.type] << "\n  connections:";
                bool connected = false;
                for (std::size_t j = 0; j < locations.size(); ++j) {
                    if (matrix[i * locations.size() + j]) {
                        file << "\n    - " << locations[j].id;
                        connected = true;
                    }
                }
                if (!connected) {
                    file << " []";
                }
                file << "\n";
            }
            file.close();
        }

        // graphdot output
        if (settings["output"].has("graphdot")) {
            std::ofstream file(settings["output"]["graphdot"].as<std::string>());
            file << "graph {\n";
            for (const auto& location : locations) {
                file << "    // " << location.id << " [label=\"" << location.name << "\"];\n";
            }
            for (std::size_t i = 0; i < locations.size(); ++i) {
                for (std::size_t j = 0; j < locations.size(); ++j) {
                    if (matrix[i * locations.size() + j]) {
                        file << "    " << locations[i].id << " -- " << locations[j].id << ";\n";
                    }
                }
            }
            file << "}\n";
            file.close();
        }

        // NetCDF output
        if (settings["output"].has("netcdf")) {
            netCDF::NcFile file(settings["output"]["netcdf"].as<std::string>(), netCDF::NcFile::replace, netCDF::NcFile::nc4);

            auto index_dim = file.addDim("index", locations.size());

            {
                auto type_dim = file.addDim("typeindex", type_map.size());
                auto var = file.addVar("typeindex", netCDF::ncString, type_dim);
                var.putVar(&type_map[0]);
            }

            {
                auto var = file.addVar("index", netCDF::ncString, index_dim);
                std::vector<const char*> c;
                c.reserve(locations.size());
                for (const auto& location : locations) {
                    c.push_back(location.id.c_str());
                }
                var.putVar(&c[0]);
            }

            {
                auto var = file.addVar("name", netCDF::ncString, index_dim);
                std::vector<const char*> c;
                c.reserve(locations.size());
                for (const auto& location : locations) {
                    c.push_back(location.name.c_str());
                }
                var.putVar(&c[0]);
            }

            {
                auto var = file.addVar("type", netCDF::ncByte, index_dim);
                std::vector<unsigned char> c;
                c.reserve(locations.size());
                for (const auto& location : locations) {
                    c.push_back(location.type);
                }
                var.putVar(&c[0]);
            }

            {
                auto var = file.addVar("latitude", netCDF::ncDouble, index_dim);
                var.putAtt("units", "degrees_north");
                std::vector<double> c;
                c.reserve(locations.size());
                for (const auto& location : locations) {
                    c.push_back(location.lat);
                }
                var.putVar(&c[0]);
            }

            {
                auto var = file.addVar("longitude", netCDF::ncDouble, index_dim);
                var.putAtt("units", "degrees_east");
                std::vector<double> c;
                c.reserve(locations.size());
                for (const auto& location : locations) {
                    c.push_back(location.lon);
                }
                var.putVar(&c[0]);
            }

            {
                auto var = file.addVar("connections", netCDF::ncByte, {index_dim, index_dim});
                var.putVar(&matrix[0]);
            }

            file.close();
        }

        // Clean up before closing the file
        for (auto& location : locations) {
            if (location.feature) {
                OGRFeature::DestroyFeature(location.feature);
            }
        }
        locations.clear();

        GDALClose(infile);

        std::cout << "Done" << std::endl;

#ifndef DEBUG
    } catch (std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return 255;
    }
#endif

    return 0;
}
