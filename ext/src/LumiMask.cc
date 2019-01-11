#include "LumiMask.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <iostream>

LumiMask LumiMask::fromJSON(const std::string& file, LumiMask::Run firstRun, LumiMask::Run lastRun)
{
  const bool noRunFilter = ( firstRun == 0 ) && ( lastRun == 0 );
  boost::property_tree::ptree ptree;
  boost::property_tree::read_json(file, ptree);

  std::vector<LumiMask::LumiBlockRange> accept;
  for ( const auto& runEntry : ptree ) {
    const LumiMask::Run run = std::stoul(runEntry.first);
    if ( noRunFilter || ( ( firstRun <= run ) && ( run <= lastRun ) ) ) {
      for ( const auto& lrEntry : runEntry.second ) {
        const auto lrNd = lrEntry.second;
        if ( lrNd.size() != 2 ) {
          std::cout << "ERROR: format should be run : list-of-lists, where each sublist has exactly two entries" << std::endl;
        } else {
          const LumiMask::LumiBlock firstLumi = std::stoul(lrNd.begin()->second.data());
          const LumiMask::LumiBlock lastLumi  = std::stoul((++lrNd.begin())->second.data());
          // std::cout << run << " : "<< firstLumi << "-" << lastLumi << std::endl;
          accept.emplace_back(run, firstLumi, lastLumi);
        }
      }
    // } else {
    //   std::cout << "Excluding run " << run << " (out of range " << firstRun << "-" << lastRun << ")" << std::endl;
    }
  }
  return LumiMask(accept);
}
