/*
 * LSST Data Management System
 * Copyright 2013 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

/// \file
/// \brief Print the layout of partitions for the specified
///        configuration.

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <cmath>

#include "boost/program_options.hpp"
#include "boost/shared_ptr.hpp"

#include "lsst/partition/Chunker.h"

#include "lsst/partition/Geometry.h"

namespace po    = boost::program_options;
namespace part  = lsst::partition;

namespace {

    typedef std::map<int32_t, std::string> Chunk2WorkerMap;

    Chunk2WorkerMap parseChunk2WorkerMap(std::string const & filename) {

        Chunk2WorkerMap result;

        std::ifstream infile {filename,  std::ifstream::in};
        for (std::string line; std::getline(infile, line);) {

            std::stringstream is(line);
            int32_t chunk;
            std::string node;
            is >> chunk;
            is >> node;
            
            result[chunk] = node;
        }
        return result;
    }

    void dumpHeader() {

        std::cout
            << "\n"
            << "        |        RA [degree]          |      DECL [degree]          |              |            \n"
            << "     id |--------------+--------------+--------------+--------------|    Area [sr] |     Worker \n"
            << "        |          Min |          Max |          Min |          Max |              |            \n"
            << " -------+--------------+--------------+--------------+--------------+--------------+------------\n";
    }
    void dumpFooter() {

        std::cout
            << "\n";
    }

    void dumpRow(
        int32_t            const   chunkId,
        part::SphericalBox const & box,
        Chunk2WorkerMap    const & chunk2worker) {

        const double area = box.area();
        if (!(std::isnormal(area) && area > 1e-7)) return;


        std::cout << "  ";

        std::cout.width(5);
        std::cout << chunkId;
        std::cout << " | ";

        std::cout.width(12);
        std::cout.precision(3);
        std::cout << std::fixed << box.getLonMin();
        std::cout << " | ";

        std::cout.width(12);
        std::cout.precision(3);
        std::cout << std::fixed << box.getLonMax();
        std::cout << " | ";

        std::cout.width(12);
        std::cout.precision(3);
        std::cout << std::fixed << box.getLatMin();
        std::cout << " | ";

        std::cout.width(12);
        std::cout.precision(3);
        std::cout << std::fixed << box.getLatMax();
        std::cout << " | ";

        std::cout.width(12);
        std::cout.precision(6);
        std::cout << std::fixed << box.area();
        std::cout << " | ";

        Chunk2WorkerMap::const_iterator const & itr = chunk2worker.find(chunkId);
        const std::string worker{itr == chunk2worker.end() ? "" : itr->second};

        std::cout.width(10);
        std::cout << worker;

        std::cout << "\n";
    }
}

static char const * help =
    "The tool will report a layout of partitions for the specified\n"
    "configuration of stripes and overlaps.\n";

int main(int argc, char const * const * argv) {

    try {
        po::options_description desc("\\_______________ Layout", 80);
        desc.add_options()

        ("help,h",    help)
        ("verbose,v", "Produce verbose output.")

        ("part.num-stripes",     po::value<int>()   ->default_value(85),    "Chunk file name prefix.")
        ("part.num-sub-stripes", po::value<int>()   ->default_value(12),    "The number of sub-stripes to divide each stripe into.")
        ("part.overlap",         po::value<double>()->default_value(0.01),  "Chunk/sub-chunk overlap radius (deg).")

        ("chunk2worker", po::value<std::string>(), "Chunk-to-worker map.")

        ("chunk", po::value<std::vector<int32_t>>(), "Chunk identifier.");

        po::positional_options_description pos_descr;
        pos_descr.add("chunk", -1);

        po::variables_map vm;
        po::store(
            po::command_line_parser(argc, argv).options(desc).positional(pos_descr).run(),
            vm);
        po::notify(vm);

        const bool verbose = vm.count("verbose") > 0;
        
        const int    numStripes             = vm["part.num-stripes"]    .as<int>();
        const int    numSubStripesPerStripe = vm["part.num-sub-stripes"].as<int>();
        const double overlap                = vm["part.overlap"]        .as<double>();

        if (verbose) {
            std::cout
                << "\n"
                << "** Configuration **\n"
                << "\n"
                << "  part.num-stripes:     " << numStripes << "\n"
                << "  part.num-sub-stripes: " << numSubStripesPerStripe << "\n"
                << "  part.overlap:         " << overlap << "\n"
                << std::endl;
        }

        std::vector<int32_t> chunks =
            vm.count("chunk") ? vm["chunk"].as<std::vector<int32_t>>() : std::vector<int32_t>();

        if (chunks.empty())
            for (int32_t chunkId = 0 ; chunkId < 20000; ++chunkId)
                chunks.push_back(chunkId);
                
        ::Chunk2WorkerMap chunk2worker;
        if (vm.count("chunk2worker")) chunk2worker =
            ::parseChunk2WorkerMap(vm["chunk2worker"].as<std::string>());

        if (verbose) {
            std::cout << "  chunk2worker size: " << chunk2worker.size() << "\n";
        }

        part::Chunker chunker{overlap, numStripes, numSubStripesPerStripe};

        ::dumpHeader();

        for (int32_t chunkId : chunks) {
            try {
                ::dumpRow(chunkId, chunker.getChunkBounds(chunkId), chunk2worker);
            } catch (std::runtime_error const & ex) {
                // src/Geometry.cc: throw std::runtime_error("Spherical box longitude angle max < min.");
                ;
            }
        }

        ::dumpFooter();

    } catch (std::exception const & ex) {
        std::cerr << ex.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
