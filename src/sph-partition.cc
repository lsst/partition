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
/// \brief The partitioner for tables which have a single
///        partitioning position.

#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"
#include "boost/shared_ptr.hpp"

#include "lsst/partition/Chunker.h"
#include "lsst/partition/ChunkReducer.h"
#include "lsst/partition/CmdLineUtils.h"
#include "lsst/partition/Csv.h"

namespace fs = boost::filesystem;
namespace po = boost::program_options;


namespace lsst {
namespace partition {

class Worker : public ChunkReducer {
public:
    Worker(po::variables_map const & vm);

    /// Compute all partitioning locations of each input
    /// record and store an output record per-location.
    void map(char const * const begin, char const * const end, Silo & silo);

    static void defineOptions(po::options_description & opts);

private:
    csv::Editor _editor;
    std::pair<int,int> _pos;
    int _chunkIdField;
    int _subChunkIdField;
    Chunker _chunker;
    std::vector<ChunkLocation> _locations;
};

Worker::Worker(po::variables_map const & vm) :
    ChunkReducer(vm),
    _editor(vm),
    _pos(-1, -1),
    _chunkIdField(-1),
    _subChunkIdField(-1),
    _chunker(vm)
{
    // Map field names of interest to field indexes.
    if (vm.count("part.pos") == 0) {
        throw std::runtime_error("The --part.pos option was not specified.");
    }
    FieldNameResolver fields(_editor);
    std::string s = vm["part.pos"].as<std::string>();
    std::pair<std::string, std::string> p = parseFieldNamePair("part.pos", s);
    _pos.first = fields.resolve("part.pos", s, p.first);
    _pos.second = fields.resolve("part.pos", s, p.second);
    if (vm.count("part.chunk") != 0) {
        s = vm["part.chunk"].as<std::string>();
        _chunkIdField = fields.resolve("part.chunk", s);
    }
    s = vm["part.sub-chunk"].as<std::string>();
    _subChunkIdField = fields.resolve("part.sub-chunk", s);
}

void Worker::map(char const * const begin,
                 char const * const end,
                 Worker::Silo & silo)
{
    typedef std::vector<ChunkLocation>::const_iterator LocIter;
    std::pair<double, double> sc;
    char const * cur = begin;
    while (cur < end) {
        cur = _editor.readRecord(cur, end);
        sc.first = _editor.get<double>(_pos.first);
        sc.second = _editor.get<double>(_pos.second);
        // Locate partitioning position and output a record for each location.
        _locations.clear();
        _chunker.locate(sc, -1, _locations);
        assert(!_locations.empty());
        for (LocIter i = _locations.begin(), e = _locations.end(); i != e; ++i) {
            _editor.set(_chunkIdField, i->chunkId);
            _editor.set(_subChunkIdField, i->subChunkId);
            silo.add(*i, _editor);
        }
    }
}

void Worker::defineOptions(po::options_description & opts) {
    po::options_description part("\\_______________ Partitioning", 80);
    part.add_options()
        ("part.prefix",
         po::value<std::string>()->default_value("chunk"),
         "Chunk file name prefix.")
        ("part.chunk",
         po::value<std::string>(),
         "Optional chunk ID output field name. This field name is appended "
         "to the output field name list if it isn't already included.")
        ("part.sub-chunk",
         po::value<std::string>()->default_value("subChunkId"),
         "Sub-chunk ID output field name. This field field name is appended "
         "to the output field name list if it isn't already included.")
        ("part.pos",
         po::value<std::string>(),
         "The partitioning longitude and latitude angle field names, "
         "separated by a comma.");
    Chunker::defineOptions(part);
    opts.add(part);
    defineOutputOptions(opts);
    csv::Editor::defineOptions(opts);
    defineInputOptions(opts);
}

typedef Job<Worker> PartitionJob;

}} // namespace lsst::partition


static char const * help =
    "The spherical partitioner partitions one or more input CSV files in\n"
    "preparation for loading into database worker nodes. This boils down to\n"
    "assigning each input position to locations in a 2-level subdivision\n"
    "scheme, where a location consists of a chunk and sub-chunk ID, and\n"
    "then bucket-sorting input records into output files by chunk ID.\n"
    "Chunk files can then be distributed to worker nodes for loading.\n"
    "\n"
    "A partitioned data-set can be built-up incrementally by running the\n"
    "partitioner with disjoint input file sets and the same output directory.\n"
    "Beware - the output CSV format, partitioning parameters, and worker\n"
    "node count MUST be identical between runs. Additionally, only one\n"
    "partitioner process should write to a given output directory at a\n"
    "time. If any of these conditions are not met, then the resulting\n"
    "chunk files will be corrupt and/or useless.\n";

int main(int argc, char const * const * argv) {
    namespace part = lsst::partition;

    try {
        po::options_description options;
        part::PartitionJob::defineOptions(options);
        po::variables_map vm;
        part::parseCommandLine(vm, options, argc, argv, help);
        part::ensureOutputFieldExists(vm, "part.chunk");
        part::ensureOutputFieldExists(vm, "part.sub-chunk");
        part::makeOutputDirectory(vm, true);
        part::PartitionJob job(vm);
        boost::shared_ptr<part::ChunkIndex> index =
            job.run(part::makeInputLines(vm));
        if (!index->empty()) {
            fs::path d(vm["out.dir"].as<std::string>());
            fs::path f = vm["part.prefix"].as<std::string>() + "_index.bin";
            index->write(d / f, false);
        }
        if (vm.count("verbose") != 0) {
            index->write(std::cout, 0);
            std::cout << std::endl;
        } else {
            std::cout << *index << std::endl;
        }
    } catch (std::exception const & ex) {
        std::cerr << ex.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

// FIXME(smm): The partitioner should store essential parameters so that
//             it can detect whether the same ones are used by incremental
//             additions to a partitioned data-set.

