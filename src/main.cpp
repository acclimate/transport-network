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
#include <stdexcept>
#include <vector>
#include "version.h"
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
#include "tqdm/tqdm.h"
#endif

inline double distance(const OGRPoint& p1, const OGRPoint& p2) {
    const auto R = 6371;
    const auto PI = 3.14159265;
    const auto latsin = std::sin((p1.getY() - p2.getY()) * PI / 360);
    const auto lonsin = std::sin((p1.getX() - p2.getX()) * PI / 360);
    const auto a = latsin * latsin + cos(p1.getY() * PI / 180) * cos(p2.getY() * PI / 180) * lonsin * lonsin;
    return 2 * R * std::atan2(std::sqrt(a), std::sqrt(1 - a));
}

static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int c) {
                                                   return !std::isspace(c);
                                               }).base(), s.end());
}

int main(int argc, char* argv[]) {
#ifndef DEBUG
    try {
#endif
        GDALAllRegister();

        const std::string shapefilename = argv[1];
        const std::string layername = argv[2];

        auto infile = static_cast<GDALDataset*>(GDALOpenEx(shapefilename.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr));
        if (!infile) {
            throw std::runtime_error("could not open shape file");
        }
        OGRLayer* inlayer = infile->GetLayerByName(layername.c_str());
        if (!inlayer) {
            throw std::runtime_error("could not read layer from shape file");
        }

        const std::size_t size = inlayer->GetFeatureCount();
        std::vector<OGRFeature*> features(size);
        std::vector<OGRGeometry*> geometries(size);
        std::vector<OGRPoint> centroids(size);
        std::vector<std::string> ids(size);
        std::vector<bool> connected(size, false);

        std::ofstream file("graph.dot");
        file << "graph {\n";
        {
#ifdef ACCLIMATE_TRANSPORT_WITH_TQDM
            tqdm::Params p;
            p.desc = "Input";
            p.ascii = "";
            p.f = stdout;
            tqdm::RangeTqdm<std::size_t> it{tqdm::RangeIterator<std::size_t>(size), tqdm::RangeIterator<std::size_t>(size, size), p};
#endif

            inlayer->ResetReading();
#pragma omp parallel for default(shared) schedule(guided)
            for (std::size_t i = 0; i < size; ++i) {
                OGRFeature* feature;
#pragma omp critical(feature)
                { feature = inlayer->GetFeature(i + 1); }
                features[i] = feature;
                OGRGeometry* geometry = feature->GetGeometryRef();  //->SimplifyPreserveTopology(0.001);
                geometry->Centroid(&centroids[i]);
                geometries[i] = geometry;  //->ConvexHull();
                std::string id = feature->GetFieldAsString(argv[3]);
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
                // std::cout << id << " " << centroids[i].getX() << " " << centroids[i].getY() << std::endl;
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
                    //(distance(centroids[i], centroids[j]) < 4000) &&
                        geometry1->Intersects(geometry2);// || geometry1->Touches(geometry2));
                if (touches) {
                    connected[i] = true;
                    connected[j] = true;
#pragma omp critical(foutput)
                    { file << "    " << id1 << " -- " << id2 << ";\n" << std::flush; }
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

        file << "}\n";
        file.close();

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
